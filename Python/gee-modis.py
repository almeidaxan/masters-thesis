# ------------------------------------------------------------------------------------
# Copyright 2016 Alexandre E. Almeida, Nathalia M. C. dos Santos. All rights reserved.
# ------------------------------------------------------------------------------------
# Description:
# 	- Download Modis Terra & Aqua images and clip them using a shapefile
#	- Output is a .zip containing .tif images for every band
# ------------------------------------------------------------------------------------

import argparse
import ee
import os
import shapefile
import socket
import sys
import time
import wget
import zipfile
from calendar import monthrange
from datetime import datetime

def valid_date(s):
	try:
		return datetime.strptime(s, '%Y-%m-%d')
	except ValueError:
		msg = "Not a valid date: '{0}'.".format(s)
		raise argparse.ArgumentTypeError(msg)

# argparse
parser = argparse.ArgumentParser(description='Downloads TERRA/AQUA MODIS images.')
parser.add_argument('shape', help='shapefile of the region of interest (do not include the extension); shape projection must be longlat')
parser.add_argument('product', help='choose a MODIS product to be downloaded', choices=['vi', 'cloud'])
parser.add_argument('periodStart', help='period start date (format YYYY-MM-DD)', type=valid_date)
parser.add_argument('periodEnd', help='period end date (format YYYY-MM-DD)', type=valid_date)
args = parser.parse_args()

shape = args.shape
product = args.product
per_start = args.periodStart
per_end = args.periodEnd

if(per_start < valid_date('2000-02-01')):
	per_start = valid_date('2000-02-01')

# earth engine init
ee.Initialize()

# chdir to shape folder and reads input shape
os.chdir('shape')
sf = shapefile.Reader(shape)
shapes = sf.shapes()
shapespoints = shapes[0].points
os.chdir('..'+os.sep)

# provided shape cannot be a MultiPolygon
if(len(shapes) > 1):
	sys.exit('The provided shape cannot be a MultiPolygon type. Please join all Polygons into a single Polygon.')

# creating a polygon using the input shape
tmp = list()
for i in range(0, len(shapespoints)-1):
	tmp.append(list(shapespoints[i]))
bnd = ee.Geometry.Polygon(tmp)

# if the provided shapefile is too complex, use the bounding box instead
# OBS: still not sure why this error happens, the only conclusion so far is
#	that the socket.error is called when the shapefile is too complex e.g.
#	when a lot of small Polygons have been joined into a single Polygon
try:
	bnd.getInfo()
except socket.error:
	bnd = ee.Geometry.Rectangle(shapes[0].bbox)
	tmp = bnd.coordinates().getInfo()

# chdir to raster folder
if(not(os.path.exists('raster'))):
	os.mkdir('raster')
os.chdir('raster')

# creating shape subfolder inside raster
if(not(os.path.exists(shape))):
	os.mkdir(shape)
if(not(os.path.exists(shape+os.sep+'SRTMGL1_003'))):
	os.mkdir(shape+os.sep+'SRTMGL1_003')
os.chdir(shape)

# function to calculate VIs using reflectance bands & reprojects to WGS 84
def getVi(image):
	# scale existing bands with scaling factor
	evi = image.select(['EVI']).multiply(0.0001)
	ndvi = image.select(['NDVI']).multiply(0.0001)
	blue = image.select(['sur_refl_b03']).multiply(0.0001)
	vzen = image.select(['ViewZenith']).multiply(0.01)
	# calculating evi2
	evi2 = image.expression('2.5*(NIR*0.0001 - R*0.0001) / (NIR*0.0001 + 2.4*R*0.0001 + 1)', {'R': image.select(['sur_refl_b01']), 'NIR': image.select(['sur_refl_b02'])}).clamp(-1,1)
	blueMask = blue.lte(0.1);
	vzenMask = blue.lte(32.5);
	image = (
		image
		.addBands(evi.select([0],['evi']))
		.addBands(evi2.select([0],['evi2']))
		.addBands(ndvi.select([0],['ndvi']))
		.reproject('EPSG:4326', scale=250)
		.updateMask(blueMask)
		.updateMask(vzenMask)
	)
	return(image)

def mod09Cloud(image):
	return image.select(['state_1km']).expression('((b(0)/1024)%2)').reproject('EPSG:4326', scale=1000)

def mod35Cloud(image):
	return image.select(['state_1km']).expression('(b(0)%4)==1|(b(0)%4)==2').reproject('EPSG:4326', scale=1000)

# {vi} - MOD13 Terra & Aqua 250m datasets
if(product == 'vi'):
	for sat in ['MOD13Q1', 'MYD13Q1']:
		# create a subfolder with the satellite dataset name
		if(not(os.path.exists(sat))):
			os.mkdir(sat)
		os.chdir(sat)

		# define satellite images to work with
		imgCol = ee.ImageCollection('MODIS/'+sat).filterBounds(bnd).filterDate(per_start, per_end).map(getVi)

		# quantity of images in imgCol
		imgColLen = imgCol.size().getInfo()

		# if no image is available, warns the user and exits
		if(imgColLen == 0):
			print('\nNo images from '+sat+' match the criteria for the selected period.')
			os.chdir('..'+os.sep)
			continue

		# creates a list with all images from imageCollection
		imgList = imgCol.toList(imgColLen)

		# creates a list with the name of all images in imgList
		imgListNames = []
		for i in range(0, imgColLen):
			imgListNames += [ee.Image(imgList.get(i)).get('system:index')]
		imgListNames = ee.List(imgListNames)
		imgListNames = imgListNames.getInfo()

		# start the download task
		dlTimes = []
		for i in range(0, imgColLen):
			# count individual download times for each file
			start_time = time.time()
			# verifiy if the file already exists
			arq = imgListNames[i]+'.zip'
			cond = os.path.exists(arq)
			if(cond):
				# if yes, prints a message and skip
				print(str(i+1).zfill(len(str(imgColLen)))+' of '+str(imgColLen)+') '+arq+' - Already downloaded!')
				continue
			else:
				# if not, downloads the file
				print(str(i+1).zfill(len(str(imgColLen)))+' of '+str(imgColLen)+') '+arq+' - Downloading...')
				# tries to connect to the server in 'tryMax' times in a row
				# if it fails, gives up and warns user
				tryNo = 1
				tryMax = 3
				while(not(cond) and tryNo<=tryMax):
					try:
						path = ee.Image(imgList.get(i)).getDownloadUrl({
							'scale': 250,
							'bands': [{'id': 'evi'}, {'id': 'evi2'}, {'id': 'ndvi'}, {'id': 'blue'}, {'id': 'vzen'}],
							'region': str(tmp)
						})
						ignore = wget.download(url=path, bar=None, out=arq)
						cond = os.path.exists(arq)
						if(cond):
							try:
								zipfile.ZipFile(arq)
							except zipfile.BadZipfile:
								cond = False
								os.remove(arq)
					except ee.ee_exception.EEException:
						tryNo += 1
						sys.stdout.write('\033[F\033[K')
						print(str(i+1).zfill(len(str(imgColLen)))+' of '+str(imgColLen)+') '+arq+' - Connection to the server timed out. Trying to download again... (Try '+str(tryNo)+' of '+str(tryMax)+')')
				elapsed_time = time.time() - start_time
				dlTimes.append(elapsed_time)
				sys.stdout.write('\033[F\033[K')
				if(cond):
					print(str(i+1).zfill(len(str(imgColLen)))+' of '+str(imgColLen)+') '+arq+' - Download successful in '+str(int(elapsed_time))+' seconds.')
				else:
					print(str(i+1).zfill(len(str(imgColLen)))+' of '+str(imgColLen)+') '+arq+' - Download failed. Please retry again later. ('+str(int(elapsed_time))+' seconds elapsed)')

		# if at least one csv has been downloaded, print time statistics
		if(dlTimes):
			print('\nTIME STATISTICS:')
			print('Total elapsed time: '+str(int(sum(dlTimes)))+' seconds')
			print('Average time per file: '+str(int(sum(dlTimes)/len(dlTimes)))+' seconds')
			print('Fastest/slowest time: '+str(int(min(dlTimes)))+' seconds'+' / '+str(int(max(dlTimes)))+' seconds')

		os.chdir('..'+os.sep)

	# also download SRTM Digital Elevation 30m data for the region
	os.chdir('SRTMGL1_003')
	arq = shape+'_elevation'+'.zip'
	cond = os.path.exists(arq)
	if(cond):
		print('\nElevation) '+arq+' - Already downloaded!')
	else:
		print('\nElevation) '+arq+' - Downloading elevation data...')
		tryNo = 1
		tryMax = 3
		while(not(cond) and tryNo<=tryMax):
			try:
				path = ee.Image('USGS/SRTMGL1_003').clip(bnd).getDownloadUrl({
					'scale': 30,
					'region': str(tmp)
				})
				ignore = wget.download(url=path, bar=None, out=arq)
				cond = os.path.exists(arq)
			except ee.ee_exception.EEException:
				tryNo += 1
				sys.stdout.write('\033[F\033[K')
				print('Elevation) '+arq+' - Connection to the server timed out. Trying to download again... (Try '+str(tryNo)+' of '+str(tryMax)+')')
		sys.stdout.write('\033[F\033[K')
		if(cond):
			print('Elevation) '+arq+' - Download successful.')
		else:
			print('Elevation) '+arq+' - Download failed. Please retry again later.')

# {cloud} - MOD09GA Terra dataset
if(product == 'cloud'):
	sat = 'MOD09GA'

	# create a subfolder with the satellite dataset name
	if(not(os.path.exists(sat))):
		os.mkdir(sat)
	os.chdir(sat)

	# start the download task
	dlTimes = []
	i = 0
	imgLen = (len(range(per_start.year, per_end.year+1))-2)*12 + (12 - per_start.month + 1) + per_end.month
	for year in range(per_start.year, per_end.year+1):
		m1 = 1
		m2 = 12
		if(year == per_start.year):
			m1 = per_start.month
		if(year == per_end.year):
			m2 = per_end.month
		for month in range(m1, m2+1):
			fileName = str(year)+'-'+str(month).zfill(2)

			date_start = valid_date(str(year)+'-'+str(month)+'-01')
			date_end = valid_date(str(year)+'-'+str(month)+'-'+str(monthrange(year, month)[1]))

			mod09 = ee.ImageCollection('MODIS/'+sat).filterDate(date_start, date_end).filterBounds(bnd).map(mod09Cloud)
			mod35 = ee.ImageCollection('MODIS/'+sat).filterDate(date_start, date_end).filterBounds(bnd).map(mod35Cloud)

			# define mean reducer and map of constants
			MEAN = ee.call('Reducer.mean')
			c100 = ee.Image(100) # to multiply by 100

			# calculate mean cloudiness (%), rename, and convert to integer
			mod09a = mod09.reduce(MEAN).select([0], ['MOD09']).multiply(c100).int8()
			mod35a = mod35.reduce(MEAN).select([0], ['MOD35']).multiply(c100).int8()

			# stack both cloud layers bands on an image
			img = mod09a.addBands(mod35a)

			# count individual download times for each file
			start_time = time.time()
			# verify if the file already exists
			arq = fileName+'.zip'
			cond = os.path.exists(arq)
			if(cond):
				# if yes, prints a message and skip
				print(str(i+1).zfill(len(str(imgLen)))+' of '+str(imgLen)+') '+arq+' - Already downloaded!')
				i = i + 1
				continue
			else:
				# if not, downloads the file
				print(str(i+1).zfill(len(str(imgLen)))+' of '+str(imgLen)+') '+arq+' - Downloading...')
				# tries to connect to the server in 'tryMax' times in a row
				# if it fails, gives up and warns user
				tryNo = 1
				tryMax = 3
				while(not(cond) and tryNo<=tryMax):
					try:
						path = ee.Image(img).getDownloadUrl({
							'scale': 1000,
							'bands': [{'id': 'MOD09'}, {'id': 'MOD35'}],
							'region': str(tmp)
						})
						ignore = wget.download(url=path, bar=None, out=arq)
						cond = os.path.exists(arq)
						if(cond):
							try:
								zipfile.ZipFile(arq)
							except zipfile.BadZipfile:
								cond = False
								os.remove(arq)
					except ee.ee_exception.EEException:
						tryNo += 1
						sys.stdout.write('\033[F\033[K')
						print(str(i+1).zfill(len(str(imgLen)))+' of '+str(imgLen)+') '+arq+' - Connection to the server timed out. Trying to download again... (Try '+str(tryNo)+' of '+str(tryMax)+')')
				elapsed_time = time.time() - start_time
				dlTimes.append(elapsed_time)
				sys.stdout.write('\033[F\033[K')
				if(cond):
					print(str(i+1).zfill(len(str(imgLen)))+' of '+str(imgLen)+') '+arq+' - Download successful in '+str(int(elapsed_time))+' seconds.')
				else:
					print(str(i+1).zfill(len(str(imgLen)))+' of '+str(imgLen)+') '+arq+' - Download failed. Please retry again later. ('+str(int(elapsed_time))+' seconds elapsed)')
			i = i + 1

	# print time statistics
	if(dlTimes):
		print('\nTIME STATISTICS:')
		print('Total elapsed time: '+str(int(sum(dlTimes)))+' seconds')
		print('Average time per file: '+str(int(sum(dlTimes)/len(dlTimes)))+' seconds')
		print('Fastest/slowest time: '+str(int(min(dlTimes)))+' seconds'+' / '+str(int(max(dlTimes)))+' seconds')
