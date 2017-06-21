# ------------------------------------------------------------------------------------
# Copyright 2016 Alexandre E. Almeida, Nathalia M. C. dos Santos, all rights reserved.
# ------------------------------------------------------------------------------------
# Description:
# 	- Download Landsat pixels time series (for the whole period available)
#	- Output is a .csv
# ------------------------------------------------------------------------------------

import argparse
import csv
import ee
import os
import signal
import ssl
import sys
import thread
import threading
import time
from datetime import datetime
from operator import mul

# argparse
satsSR = ['LT4_SR', 'LT5_SR', 'LE7_SR', 'LC8_SR']
satsTOA = ['LT4_L1T_TOA_FMASK', 'LT5_L1T_TOA_FMASK', 'LE7_L1T_TOA_FMASK', 'LC8_L1T_TOA_FMASK']
satsNames = ['4', '5', '7', '8']
satsProd = ['SR', 'TOA']
parser = argparse.ArgumentParser(description='Downloads Landsat time series for individual pixels, covering the whole period available.')
parser.add_argument('csvName', help='csv file name (do not include the extension)')
parser.add_argument('satellite', help='Which Landsat # to use', choices=satsNames)
parser.add_argument('satprod', help='Which Landsat product to use', choices=satsProd)
args = parser.parse_args()

csvName = args.csvName
if(args.satprod == 'SR'):
	sat = satsSR[satsNames.index(args.satellite)]
if(args.satprod == 'TOA'):
	sat = satsTOA[satsNames.index(args.satellite)]

# earth engine API init
ee.Initialize()

# convert Julian day to YYYY-MM-DD date
def julianDayToDate(s):
	return datetime.strptime(s, '%Y%j').strftime('%Y_%m_%d')

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

# function to extract the time series of the pixel {vi}
def getViPx(img):
	val = img.sample(region=bounds, scale=30, numPixels=1).first()
	return(val)

# defines csv that contains the plots coordinates
folName = csvName
csvName += '.csv'

# create a folder inside csv folder with the csvName name
os.chdir('csv')
if not(os.path.exists(folName)):
	os.mkdir(folName)

# define image collection
imgCol = ee.ImageCollection('LANDSAT/'+sat)

# create a subfolder with the satellite dataset name
os.chdir(folName)
if not(os.path.exists(sat)):
	os.mkdir(sat)
os.chdir('..'+os.sep)

# calculate number of csv rows (containing coords)
# to iterate, excluding header
with open(csvName, 'rb') as csvfile:
	maxRows = len(csvfile.readlines())-1

# calculate longest file name length in the csv
maxLen = 0
with open(csvName, 'rb') as csvfile:
	csvreader = csv.reader(csvfile)
	csvreader.next()
	for row in csvreader:
		if len(row[0]) > maxLen:
			maxLen = len(row[0])

# start the download task
print('')
dlTimes = []
with open(csvName, 'rb') as csvfile:
	csvreader = csv.reader(csvfile)
	# skips the header row
	csvreader.next()
	# calculate longest file name
	for row in csvreader:
		plotName = row[0]
		arq = folName+os.sep+sat+os.sep+plotName+'.csv'
		arqPrint = folName+os.sep+sat+os.sep+(plotName+'.csv').ljust(maxLen+4)
		cond = os.path.exists(arq)
		# checks if csv has already been downloaded
		if cond: # if cond and os.stat(arq).st_size>20:
			# if csv has already been downloaded, warns the user and skip
			print(str(csvreader.line_num-1).zfill(len(str(maxRows)))+' of '+str(maxRows)+') '+arqPrint+' - File already downloaded.')
		else:
			# if not, attempts downloads the csv
			start_time = time.time()
			print(str(csvreader.line_num-1).zfill(len(str(maxRows)))+' of '+str(maxRows)+') '+arqPrint+' - Downloading...')
			latCen = float(row[1])
			longCen = float(row[2])
			bounds = ee.Geometry.Point([longCen, latCen])
			if sat == 'LC8_SR':
				values = imgCol.filterBounds(bounds).map(getVi8).map(getViPx,True).filterMetadata('cfmask', 'equals', 0)
			else:
				values = imgCol.filterBounds(bounds).map(getVi).map(getViPx,True).filterMetadata('cfmask', 'equals', 0)
			# tries to connect to the server in 'tryMax' times in a row
			# if it fails, gives up and warns user
			tryNo = 1
			tryMax = 5
			suc = False
			while tryNo <= tryMax:
				timer = threading.Timer(90, os.kill, args=[os.getpid(), signal.SIGINT])
				timer.start()
				try:
					values = values.getInfo()
					timer.cancel()
					valuesAux = values['features']
					values = []
					for i in range(len(valuesAux)):
						values += [julianDayToDate(str(valuesAux[i]['id'])[9:16])] + \
						[round(valuesAux[i]['properties']['B1'],4)] + \
						[round(valuesAux[i]['properties']['B2'],4)] + \
						[round(valuesAux[i]['properties']['B3'],4)] + \
						[round(valuesAux[i]['properties']['B4'],4)] + \
						[round(valuesAux[i]['properties']['B5'],4)]
						if sat == 'LC8_SR':
							values += [round(valuesAux[i]['properties']['B6'],4)]
						values += [round(valuesAux[i]['properties']['B7'],4)] + \
						[round(valuesAux[i]['properties']['ndvi'],4)] + \
						[round(valuesAux[i]['properties']['evi'],4)] + \
						[round(valuesAux[i]['properties']['evi2'],4)]
					suc = True
					break
				except (KeyboardInterrupt, ee.ee_exception.EEException, ssl.SSLEOFError):
					timer.cancel()
					tryNo += 1
					sys.stdout.write('\033[F\033[K')
					print(str(csvreader.line_num-1).zfill(len(str(maxRows)))+' of '+str(maxRows)+') '+arqPrint+' - Connection to the server timed out. Trying to download again... (Try '+str(tryNo)+' of '+str(tryMax)+')')
			if suc:
				with open(arq, 'wb') as csvfile:
					csvwriter = csv.writer(csvfile, delimiter=',')
					if sat == 'LC8_SR':
						csvwriter.writerow(['date', 'B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'ndvi', 'evi', 'evi2'])
					else:
						csvwriter.writerow(['date', 'B1', 'B2', 'B3', 'B4', 'B5', 'B7', 'ndvi', 'evi', 'evi2'])
					# number of columns in output csv
					if len(values) > 0:
						numCol = len(values)/len(valuesAux)
						for i in range(1, len(values)/numCol + 1):
							j = (i-1) * numCol
							# write a csv line for each timestamp
							bands = [values[1+j]] + \
									[values[2+j]] + \
									[values[3+j]] + \
									[values[4+j]] + \
									[values[5+j]] + \
									[values[6+j]]
							if sat == 'LC8_SR':
								bands += [values[7+j]]
							# only run if no saturated bands are present
							if all(b < 2 for b in bands):
								if sat == 'LC8_SR':
									csvwriter.writerow(
										[values[0+j]] + \
										bands + \
										[values[8+j]] + \
										[values[9+j]] + \
										[values[10+j]]
									)
								else:
									csvwriter.writerow(
										[values[0+j]] + \
										bands + \
										[values[7+j]] + \
										[values[8+j]] + \
										[values[9+j]]
									)
			elapsed_time = time.time() - start_time
			dlTimes.append(elapsed_time)
			sys.stdout.write('\033[F\033[K')
			if suc:
				print(str(csvreader.line_num-1).zfill(len(str(maxRows)))+' of '+str(maxRows)+') '+arqPrint+' - Download successful in '+str(int(elapsed_time))+' seconds.')
			else:
				print(str(csvreader.line_num-1).zfill(len(str(maxRows)))+' of '+str(maxRows)+') '+arqPrint+' - Download failed. Please retry again later. ('+str(int(elapsed_time))+' seconds elapsed)')

# if at least one csv has been downloaded, print time statistics
if dlTimes:
	print('\nTIME STATISTICS:')
	print('Total elapsed time: '+str(int(sum(dlTimes)))+' seconds')
	print('Average time per file: '+str(int(sum(dlTimes)/len(dlTimes)))+' seconds')
	print('Fastest/slowest time: '+str(int(min(dlTimes)))+' seconds'+' / '+str(int(max(dlTimes)))+' seconds\n')
