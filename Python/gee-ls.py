# ------------------------------------------------------------------------------------
# Copyright 2016 Alexandre E. Almeida, Nathalia M. C. dos Santos. All rights reserved.
# ------------------------------------------------------------------------------------
# Description:
# 	- Download Landsat images and clip them using a shapefile
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
from datetime import datetime

def valid_date(s):
	try:
		return datetime.strptime(s, "%Y-%m-%d")
	except ValueError:
		msg = "Not a valid date: '{0}'.".format(s)
		raise argparse.ArgumentTypeError(msg)

# argparse
satsSR = ['LT4_SR', 'LT5_SR', 'LE7_SR', 'LC8_SR']
satsTOA = ['LT4_L1T_TOA_FMASK', 'LT5_L1T_TOA_FMASK', 'LE7_L1T_TOA_FMASK', 'LC8_L1T_TOA_FMASK']
satsNames = ['4', '5', '7', '8']
satsProd = ['SR', 'TOA']
parser = argparse.ArgumentParser(description='Downloads Landsat images.')
parser.add_argument('shape', help='shapefile of the region of interest (do not include the extension); shape projection must be longlat')
parser.add_argument('satellite', help='Which Landsat # to use', choices=satsNames)
parser.add_argument('satprod', help='Which Landsat product to use', choices=satsProd)
parser.add_argument('periodStart', help='period start date (format YYYY-MM-DD)', type=valid_date)
parser.add_argument('periodEnd', help='period end date (format YYYY-MM-DD)', type=valid_date)
args = parser.parse_args()

shape = args.shape
per_start = args.periodStart
per_end = args.periodEnd
if(args.satprod == 'SR'):
	sat = satsSR[satsNames.index(args.satellite)]
if(args.satprod == 'TOA'):
	sat = satsTOA[satsNames.index(args.satellite)]

# earth engine init
ee.Initialize()

# chdir to shape folder and reads input shape
os.chdir('shape')
sf = shapefile.Reader(shape)
shapes = sf.shapes()
shapespoints = shapes[0].points
os.chdir('..'+os.sep)

# provided shapefile cannot be a MultiPolygon
if(len(shapes) > 1):
	sys.exit('The provided shape cannot be a MultiPolygon type. Please join all Polygons into a single Polygon.')

# creating a polygon from input shape
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

# chdir to raster folder
if(not(os.path.exists('raster'))):
	os.mkdir('raster')
os.chdir('raster')

# creating subfolders inside raster
if(not(os.path.exists(shape))):
	os.mkdir(shape)
if(not(os.path.exists(shape+os.sep+sat))):
	os.mkdir(shape+os.sep+sat)
if(not(os.path.exists(shape+os.sep+'SRTMGL1_003'))):
	os.mkdir(shape+os.sep+'SRTMGL1_003')
os.chdir(shape+os.sep+sat)

# function to calculate and add the VI band & cloud cover band (LS 4,5,7)
def getVi(image):
	if(args.satprod == 'SR'):
		cfmask = image.select(['cfmask']).multiply(1)
		image2 = image.select(['cfmask'],['ignore'])
		# scale existing bands with scaling factor
		scalingFactor = 0.0001
	if(args.satprod == 'TOA'):
		cfmask = image.select(['fmask']).multiply(1)
		image2 = image.select(['fmask'],['ignore'])
		# don't use scaling factor for TOA
		scalingFactor = 1
	b1 = image.select(['B1']).multiply(scalingFactor)
	b2 = image.select(['B2']).multiply(scalingFactor)
	b3 = image.select(['B3']).multiply(scalingFactor)
	b4 = image.select(['B4']).multiply(scalingFactor)
	b5 = image.select(['B5']).multiply(scalingFactor)
	b7 = image.select(['B7']).multiply(scalingFactor)
	evi = image.expression('2.5 * (NIR - R) / (NIR + 6.0*R - 7.5*B + 1)', {'R': b3, 'NIR': b4, 'B': b1}).clamp(-1,1)
	evi2 = image.expression('2.5 * (NIR - R) / (NIR + 2.4*R + 1)', {'R': b3, 'NIR': b4}).clamp(-1,1)
	ndvi = image.normalizedDifference(['B4','B3'])
	image2 = (
		image2
		.addBands(b1.select([0],['B1']))
		.addBands(b2.select([0],['B2']))
		.addBands(b3.select([0],['B3']))
		.addBands(b4.select([0],['B4']))
		.addBands(b5.select([0],['B5']))
		.addBands(b7.select([0],['B7']))
		.addBands(cfmask.select([0],['cfmask']))
		.addBands(evi.select([0],['evi']))
		.addBands(evi2.select([0],['evi2']))
		.addBands(ndvi.select([0],['ndvi']))
	)
	return(image2)

# function to calculate and add the VI band & cloud cover band (LS 8)
def getVi8(image):
	if(args.satprod == 'SR'):
		cfmask = image.select(['cfmask']).multiply(1)
		image2 = image.select(['cfmask'],['ignore'])
		# scale existing bands with scaling factor
		scalingFactor = 0.0001
	if(args.satprod == 'TOA'):
		cfmask = image.select(['fmask']).multiply(1)
		image2 = image.select(['fmask'],['ignore'])
		# don't use scaling factor for TOA
		scalingFactor = 1
	b1 = image.select(['B1']).multiply(scalingFactor)
	b2 = image.select(['B2']).multiply(scalingFactor)
	b3 = image.select(['B3']).multiply(scalingFactor)
	b4 = image.select(['B4']).multiply(scalingFactor)
	b5 = image.select(['B5']).multiply(scalingFactor)
	b6 = image.select(['B6']).multiply(scalingFactor)
	b7 = image.select(['B7']).multiply(scalingFactor)
	evi = image.expression('2.5 * (NIR - R) / (NIR + 6.0*R - 7.5*B + 1)', {'R': b4, 'NIR': b5, 'B': b2}).clamp(-1,1)
	evi2 = image.expression('2.5 * (NIR - R) / (NIR + 2.4*R + 1)', {'R': b4, 'NIR': b5}).clamp(-1,1)
	ndvi = image.normalizedDifference(['B5','B4'])
	image2 = (
		image2
		.addBands(b1.select([0],['B1']))
		.addBands(b2.select([0],['B2']))
		.addBands(b3.select([0],['B3']))
		.addBands(b4.select([0],['B4']))
		.addBands(b5.select([0],['B5']))
		.addBands(b5.select([0],['B6']))
		.addBands(b7.select([0],['B7']))
		.addBands(cfmask.select([0],['cfmask']))
		.addBands(evi.select([0],['evi']))
		.addBands(evi2.select([0],['evi2']))
		.addBands(ndvi.select([0],['ndvi']))
	)
	return(image2)

# define image collection
if(args.satellite == '8'):
	imgCol = (
		ee.ImageCollection('LANDSAT/'+sat)
		.filterBounds(bnd)
		.filterDate(per_start, per_end)
		.map(getVi8)
	)
else:
	imgCol = (
		ee.ImageCollection('LANDSAT/'+sat)
		.filterBounds(bnd)
		.filterDate(per_start, per_end)
		.map(getVi)
	)

# quantity of images in imgCol
imgColLen = imgCol.size().getInfo()

# if no image is available, warns the user and exits
if(imgColLen == 0):
	print('\nNo images from '+sat+' match the criteria for the selected period.\n')
	os.chdir('..'+os.sep)
	sys.exit()

# creates a list with all images from imgCol
imgList = imgCol.toList(imgColLen)

# creates a list with the name of all images in imgList
imgListNames = []
for i in range(0, imgColLen):
	imgListNames += [ee.Image(imgList.get(i)).get('system:index')]
imgListNames = ee.List(imgListNames)
imgListNames = imgListNames.getInfo()

# creates a list with the name of all bands of the imgCol
imgBands = ee.Image(imgList.get(1)).bandNames().getInfo()
imgBandsListDict = list()
for i, val in enumerate(imgBands):
	if(str(val) != 'ignore'):
		imgBandsListDict.append({'id': str(val)})

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
				path = ee.Image(imgList.get(i)).clip(bnd).getDownloadUrl({
					'scale': 30,
					'bands': imgBandsListDict,
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

# also download SRTM Digital Elevation 30m data for the region
os.chdir('..'+os.sep)
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
