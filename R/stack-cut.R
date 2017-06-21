# -------------------------------------------------------------------------(#HEADER)
# Script:
# 		stack-cut
# Author:
# 		Alexandre Esteves Almeida (2016)
# Description:
# 		- cut RasterStack objects according to a shapefile
# -------------------------------------------------------------------------(#IMPDEF)

# load packages, install if necessary
packs <- c("bfastSpatial",
			  "raster")

# source import-define script to load required libraries and define paths
if(Sys.info()[["sysname"]] == "Windows") {
	# windows OS
	source("D:/Documents/Dropbox/Unicamp/Projeto de Mestrado/Programas/R/import-define.R")
} else {
	# unix-based OSs
	if(Sys.getenv("RSTUDIO_USER_IDENTITY") == "almeida") {
		source("/Users/almeida/Dropbox/Unicamp/Projeto de Mestrado/Programas/R/import-define.R")
	}
	if(Sys.getenv("RSTUDIO_USER_IDENTITY") == "Menini") {
		source("/Users/almeida/Dropbox/Unicamp/Pos/Programas/R/import-define.R")
	}
}

# define file name
f_area <- "fray-jorge2"
f_band <- "evi2"
f_sat <- "LT5_L1T_TOA"
f_name <- paste0(f_area, "_", f_band, ".grd")

# -------------------------------------------------------------------------(#ZIPEXT)

# filter only zip files from wd
setwd(paste0(p_raster, f_area, "/", f_sat, "/"))
files <- dir(pattern = ".zip", full.names = T)

# output folder
dir.create("extracted/", showWarnings = F)
outPath <- paste0("./extracted")

# progress bar to follow up loop progress
print("Extracting files...")
pb <- txtProgressBar(min = 0,
							max = length(files),
							style = 3)

# loop that extracts files from zips
for (i in 1:length(files)) {
	# list the content of zip
	auxList <- unzip(files[i], list = T, overwrite = F)$Name

	# list tif files indexes inside the zips
	toUnzip <- grep(".tif", auxList)

	# extracts til files from zips
	unzip(files[i],
			files = auxList[toUnzip],
			exdir = outPath)

	# update progress bar
	setTxtProgressBar(pb, i)

	rm(auxList, toUnzip)
}

# closes unzip progress bar
close(pb)

# list all .tif files inside the extracted folder
files <- dir(path = paste0(p_raster, f_area, "/", f_sat, "/", "extracted/"),
				 pattern = paste0(f_band, ".tif"))

# reorder the files by date, if there is more than one pathrow
filesDates <- substr(files, 10, 16)
files <- files[order(filesDates)]

# create rasters for all files
rasters <- lapply(paste0(p_raster, f_area, "/", f_sat, "/", "extracted/", files),
						FUN = raster)

# if any date has more than one image, merge them
dates <- as.Date(unlist(lapply(rasters, function(x) getSceneinfo(names(x))$date)))
tableDates <- table(dates)
indpen <- 0
if(sum(tableDates) > 1) {
	whichDates <- names(which(tableDates > 1))
	for(i in 1:length(whichDates)) {
		ind <- grep(whichDates[i], dates) - indpen
		newRaster <- rasters[[ind[1]]]
		newName <- names(newRaster)
		for(j in 2:length(ind)) {
			newRaster <- mosaic(newRaster, rasters[[ind[j]]], fun = max)
			rasters[[ind[j]]] <- NULL
			indpen <- indpen + 1
		}
		rasters[[ind[1]]] <- newRaster
		names(rasters[[ind[1]]]) <- newName
	}
}

# defining the stack of all the rasters
s <- raster::stack(rasters)

# reads shapefile and transforms its crs match rasters crs
shp <- shapefile(paste0(p_shape, f_area, ".shp"))
shp <- spTransform(shp, crs(s[[1]]))

# plot the cutted rasters
plot(crop(s[[1]], shp))

# delete temporary directory with unzipped images
unlink(paste0(p_raster, f_area, "/", f_sat, "/", "extracted"),
		 recursive = T)