# ------------------------------------------------------------------------------------
# Copyright 2016 Alexandre E. Almeida, Nathalia M. C. dos Santos, all rights reserved.
# ------------------------------------------------------------------------------------
# Description:
# 	- Download Modis Terra&Aqua pixels time series (for the whole period available)
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

# argparse
parser = argparse.ArgumentParser(description='Downloads TERRA/AQUA MODIS time series for individual pixels, covering the whole period available.')
parser.add_argument('csvName', help='csv file name (do not include the extension)')
parser.add_argument('product', help='choose a MODIS product to be downloaded', choices=['vi', 'bands', 'fire'])
args = parser.parse_args()

csvName = args.csvName
product = args.product

# earth engine init
ee.Initialize()

# function to calculate VIs using reflectance bands and reproject to WGS 84
def getVi(img):
	EVI2 = img.expression('2.5*(NIR*0.0001 - R*0.0001) / (NIR*0.0001 + 2.4*R*0.0001 + 1)', {'R': img.select(['sur_refl_b01']), 'NIR': img.select(['sur_refl_b02'])}).clamp(-1,1)
	EVI2 = EVI2.select([0],['EVI2'])
	img = img.addBands(EVI2).reproject('EPSG:4326', scale=250)
	return(img)

# function to extract the time series of the pixel {bands, fire}
def getPx(img):
	img = img.reproject('EPSG:4326', scale=500)
	val = img.sample(region=bounds, scale=500, numPixels=1).first()
	return(val)

# function to extract the time series of the pixel {vi}
def getViPx(img):
	val = img.sample(region=bounds, scale=250, numPixels=1).first()
	return(val)

# defines csv that contains the plots coordinates
folName = csvName
csvName = csvName+'.csv'

# create a folder inside csv folder with the csvName name
os.chdir('csv')
if(not(os.path.exists(folName))):
	os.mkdir(folName)

print('')
# {vi} - MOD13 Terra & Aqua 250m datasets
if(product == 'vi'):
	for sat in ['MOD13Q1', 'MYD13Q1']:
		imgCol = ee.ImageCollection('MODIS/'+sat)

		# create a subfolder with the satellite dataset name
		os.chdir(folName)
		if(not(os.path.exists(sat))):
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
				if(len(row[0])>maxLen):
					maxLen = len(row[0])

		# start the download task
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
				# checks if csv has already been downloaded and if it is not an empty file (20 bytes)
				if(cond and os.stat(arq).st_size>20):
					# if csv has already been downloaded, warns the user and skip
					print(str(csvreader.line_num-1).zfill(len(str(maxRows)))+' of '+str(maxRows)+') '+arqPrint+' - File already downloaded.')
				else:
					# if not, downloads the csv
					start_time = time.time()
					print(str(csvreader.line_num-1).zfill(len(str(maxRows)))+' of '+str(maxRows)+') '+arqPrint+' - Downloading...')
					latCen = float(row[1])
					longCen = float(row[2])
					bounds = ee.Geometry.Point([longCen, latCen])
					values = imgCol.filterBounds(bounds).map(getVi).map(getViPx,True).filterMetadata('ViewZenith', 'less_than', 3250).filterMetadata('sur_refl_b03', 'less_than', 1000)
					# tries to connect to the server in 'tryMax' times in a row
					# if it fails, gives up and warns user
					tryNo = 1
					tryMax = 5
					suc = False
					while(tryNo<=tryMax):
						timer = threading.Timer(90, os.kill, args=[os.getpid(), signal.SIGINT])
						timer.start()
						try:
							values = values.getInfo()
							timer.cancel()
							valuesAux = values['features']
							values = []
							for i in range(len(valuesAux)):
								values = values + \
								[str(valuesAux[i]['id'])[12:]] + \
								[round(valuesAux[i]['properties']['NDVI']*0.0001,4)] + \
								[round(valuesAux[i]['properties']['EVI']*0.0001,4)] + \
								[round(valuesAux[i]['properties']['EVI2'],4)]
							suc = True
							break
						except ee.ee_exception.EEException, ssl.SSLEOFError:
							timer.cancel()
							tryNo += 1
							sys.stdout.write('\033[F\033[K')
							print(str(csvreader.line_num-1).zfill(len(str(maxRows)))+' of '+str(maxRows)+') '+arqPrint+' - Connection to the server timed out. Trying to download again... (Try '+str(tryNo)+' of '+str(tryMax)+')')
					if(suc):
						with open(arq, 'wb') as csvfile:
							csvwriter = csv.writer(csvfile, delimiter=',')
							csvwriter.writerow(['date', 'ndvi', 'evi', 'evi2'])
							# number of columns in output csv
							if(len(values) > 0):
								numCol = len(values)/len(valuesAux)
								for i in range(1, len(values)/numCol + 1):
									j = (i-1) * numCol
									# write a csv line for each timestamp
									csvwriter.writerow(
										[values[0+j]] + \
										[values[1+j]] + \
										[values[2+j]] + \
										[values[3+j]]
									)
					elapsed_time = time.time() - start_time
					dlTimes.append(elapsed_time)
					sys.stdout.write('\033[F\033[K')
					if(suc):
						print(str(csvreader.line_num-1).zfill(len(str(maxRows)))+' of '+str(maxRows)+') '+arqPrint+' - Download successful in '+str(int(elapsed_time))+' seconds.')
					else:
						print(str(csvreader.line_num-1).zfill(len(str(maxRows)))+' of '+str(maxRows)+') '+arqPrint+' - Download failed. Please retry again later. ('+str(int(elapsed_time))+' seconds elapsed)')

		# if at least one csv has been downloaded, print time statistics
		if(dlTimes):
			print('\nTIME STATISTICS:')
			print('Total elapsed time: '+str(int(sum(dlTimes)))+' seconds')
			print('Average time per file: '+str(int(sum(dlTimes)/len(dlTimes)))+' seconds')
			print('Fastest/slowest time: '+str(int(min(dlTimes)))+' seconds'+' / '+str(int(max(dlTimes)))+' seconds\n')

# {bands} - MOD09 Terra & Aqua 500m datasets
if(product == 'bands'):
	for sat in ['MOD09A1', 'MYD09A1']:
		imgCol = ee.ImageCollection('MODIS/'+sat)

		# create a subfolder with the satellite dataset name
		os.chdir(folName)
		if(not(os.path.exists(sat))):
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
				if(len(row[0])>maxLen):
					maxLen = len(row[0])

		# start the download task
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
				# checks if csv has already been downloaded and if it is not an empty file (20 bytes)
				if(cond and os.stat(arq).st_size>25):
					# if csv has already been downloaded, warns the user and skip
					print(str(csvreader.line_num-1).zfill(len(str(maxRows)))+' of '+str(maxRows)+') '+arqPrint+' - File already downloaded.')
				else:
					# if not, downloads the csv
					start_time = time.time()
					print(str(csvreader.line_num-1).zfill(len(str(maxRows)))+' of '+str(maxRows)+') '+arqPrint+' - Downloading...')
					latCen = float(row[1])
					longCen = float(row[2])
					bounds = ee.Geometry.Point([longCen, latCen])
					values = imgCol.filterBounds(bounds).map(getViPx,True).filterMetadata('ViewZenith', 'less_than', 3250).filterMetadata('sur_refl_b03', 'less_than', 1000)
					# tries to connect to the server in 'tryMax' times in a row
					# if it fails, gives up and warns user
					tryNo = 1
					tryMax = 5
					suc = False
					while(tryNo<=tryMax):
						timer = threading.Timer(90, os.kill, args=[os.getpid(), signal.SIGINT])
						timer.start()
						try:
							values = values.getInfo()
							timer.cancel()
							valuesAux = values['features']
							values = []
							for i in range(len(valuesAux)):
								values = values + \
								[str(valuesAux[i]['id'])[12:]] + \
								[round(valuesAux[i]['properties']['sur_refl_b01']*0.0001,4)] + \
								[round(valuesAux[i]['properties']['sur_refl_b02']*0.0001,4)] + \
								[round(valuesAux[i]['properties']['sur_refl_b03']*0.0001,4)] + \
								[round(valuesAux[i]['properties']['sur_refl_b04']*0.0001,4)] + \
								[round(valuesAux[i]['properties']['sur_refl_b05']*0.0001,4)] + \
								[round(valuesAux[i]['properties']['sur_refl_b06']*0.0001,4)] + \
								[round(valuesAux[i]['properties']['sur_refl_b07']*0.0001,4)]
							suc = True
							break
						except ee.ee_exception.EEException, ssl.SSLEOFError:
							timer.cancel()
							tryNo += 1
							sys.stdout.write('\033[F\033[K')
							print(str(csvreader.line_num-1).zfill(len(str(maxRows)))+' of '+str(maxRows)+') '+arqPrint+' - Connection to the server timed out. Trying to download again... (Try '+str(tryNo)+' of '+str(tryMax)+')')
					if(suc):
						with open(arq, 'wb') as csvfile:
							csvwriter = csv.writer(csvfile, delimiter=',')
							csvwriter.writerow(['date', 'b1', 'b2', 'b3', 'b4', 'b5', 'b6', 'b7'])
							# number of columns in output csv
							if(len(values) > 0):
								numCol = len(values)/len(valuesAux)
								for i in range(1, len(values)/numCol + 1):
									j = (i-1) * numCol
									# write a csv line for each timestamp
									csvwriter.writerow(
										[values[0+j]] + \
										[values[1+j]] + \
										[values[2+j]] + \
										[values[3+j]] + \
										[values[4+j]] + \
										[values[5+j]] + \
										[values[6+j]] + \
										[values[7+j]]
									)
					elapsed_time = time.time() - start_time
					dlTimes.append(elapsed_time)
					sys.stdout.write('\033[F\033[K')
					if(suc):
						print(str(csvreader.line_num-1).zfill(len(str(maxRows)))+' of '+str(maxRows)+') '+arqPrint+' - Download successful in '+str(int(elapsed_time))+' seconds.')
					else:
						print(str(csvreader.line_num-1).zfill(len(str(maxRows)))+' of '+str(maxRows)+') '+arqPrint+' - Download failed. Please retry again later. ('+str(int(elapsed_time))+' seconds elapsed)')

		# if at least one csv has been downloaded, print time statistics
		if(dlTimes):
			print('\nTIME STATISTICS:')
			print('Total elapsed time: '+str(int(sum(dlTimes)))+' seconds')
			print('Average time per file: '+str(int(sum(dlTimes)/len(dlTimes)))+' seconds')
			print('Fastest/slowest time: '+str(int(min(dlTimes)))+' seconds'+' / '+str(int(max(dlTimes)))+' seconds\n')

# {fire} - MCD45 Composite 500m dataset
if(product == 'fire'):
	for sat in ['MCD45A1']:
		imgCol = ee.ImageCollection('MODIS/051/'+sat)

		# create a subfolder with the satellite dataset name
		os.chdir(folName)
		if(not(os.path.exists(sat))):
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
				if(len(row[0])>maxLen):
					maxLen = len(row[0])

		# start the download task
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
				if(cond):
					# if csv has already been downloaded, warns the user and skip
					print(str(csvreader.line_num-1).zfill(len(str(maxRows)))+' of '+str(maxRows)+') '+arqPrint+' - File already downloaded.')
				else:
					# if not, downloads the csv
					start_time = time.time()
					print(str(csvreader.line_num-1).zfill(len(str(maxRows)))+' of '+str(maxRows)+') '+arqPrint+' - Downloading...')
					latCen = float(row[1])
					longCen = float(row[2])
					#latCen = -12.3241
					#longCen = -50.7355
					bounds = ee.Geometry.Point([longCen, latCen])
					values = imgCol.filterDate('2000-01-01', time.strftime('%Y-%m-%d')).map(getPx,True)
					# tries to connect to the server in 'tryMax' times in a row
					# if it fails, gives up and warns user
					tryNo = 1
					tryMax = 5
					suc = False
					while(tryNo<=tryMax):
						timer = threading.Timer(90, os.kill, args=[os.getpid(), signal.SIGINT])
						timer.start()
						try:
							values = values.getInfo()
							timer.cancel()
							valuesAux = values['features']
							values = []
							for i in range(len(valuesAux)):
								values = values + \
								[str(valuesAux[i]['id'])] + \
								[valuesAux[i]['properties']['burndate']] + \
								[valuesAux[i]['properties']['ba_qa']]
							suc = True
							break
						except ee.ee_exception.EEException, ssl.SSLEOFError:
							timer.cancel()
							tryNo += 1
							sys.stdout.write('\033[F\033[K')
							print(str(csvreader.line_num-1).zfill(len(str(maxRows)))+' of '+str(maxRows)+') '+arqPrint+' - Connection to the server timed out. Trying to download again... (Try '+str(tryNo)+' of '+str(tryMax)+')')
					if(suc):
						with open(arq, 'wb') as csvfile:
							csvwriter = csv.writer(csvfile, delimiter=',')
							csvwriter.writerow(['date', 'burndate', 'qa'])
							# number of columns in output csv
							if(len(values) > 0):
								numCol = len(values)/len(valuesAux)
								for i in range(1, len(values)/numCol + 1):
									j = (i-1) * numCol
									# write a csv line for each timestamp
									csvwriter.writerow(
										[values[0+j]] + \
										[values[1+j]] + \
										[values[2+j]]
									)
					elapsed_time = time.time() - start_time
					dlTimes.append(elapsed_time)
					sys.stdout.write('\033[F\033[K')
					if(suc):
						print(str(csvreader.line_num-1).zfill(len(str(maxRows)))+' of '+str(maxRows)+') '+arqPrint+' - Download successful in '+str(int(elapsed_time))+' seconds.')
					else:
						print(str(csvreader.line_num-1).zfill(len(str(maxRows)))+' of '+str(maxRows)+') '+arqPrint+' - Download failed. Please retry again later. ('+str(int(elapsed_time))+' seconds elapsed)')

		# if at least one csv has been downloaded, print time statistics
		if(dlTimes):
			print('\nTIME STATISTICS:')
			print('Total elapsed time: '+str(int(sum(dlTimes)))+' seconds')
			print('Average time per file: '+str(int(sum(dlTimes)/len(dlTimes)))+' seconds')
			print('Fastest/slowest time: '+str(int(min(dlTimes)))+' seconds'+' / '+str(int(max(dlTimes)))+' seconds\n')
