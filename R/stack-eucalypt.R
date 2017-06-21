require(velox)
require(raster)
require(rgdal)

scriptName <- "bfast-raster"

# load packages, install if necessary
packs <- c("bfastSpatial",
		   "lubridate",
		   "raster",
		   "rgeos",
		   "robust",
		   "leaflet")

# source import-define script to load required libraries and define paths
if(Sys.info()[["sysname"]] == "Windows") {
	# windows OS
	source("D:/Documents/Dropbox/Unicamp/Projeto de Mestrado/Programas/R/import-define.R")
} else {
	# unix-based OSs
	if(Sys.getenv("USER") == "almeida") {
		source("/Users/almeida/Dropbox/Unicamp/Projeto de Mestrado/Programas/R/import-define.R")
	}
	if(Sys.getenv("USER") == "Menini") {
		source("/Users/almeida/Dropbox/Unicamp/Pos/Programas/R/import-define.R")
	}
}

# defining file names to be loaded
f_area <- "eucalipto"
f_dir <- dir(paste0(p_raster, f_area))
f_sat <- which(startsWith(f_dir, "M")) # "M" for MODIS

for(k in 1:length(f_sat)) {
	# filter only zip files from wd
	setwd(paste0(p_raster, f_area, "/", f_dir[f_sat[k]]))
	files <- dir(pattern = ".zip", full.names = T)

	# remove .tmp files if they exist
	tmp <- dir(pattern = ".tmp", full.names = T)
	tmp <- file.remove(tmp)

	# output folder
	dir.create("../extracted-modis/", showWarnings = F)
	outPath <- paste0("../extracted-modis")

	# progress bar to follow up loop progress
	print(paste("Extracting", f_dir[f_sat[k]],"files..."))
	pb <- txtProgressBar(min = 0,
						 max = length(files),
						 style = 3)

	# loop that extracts files from zips
	for (i in 1:length(files)) {
		# list the content of zip
		auxList <- unzip(files[i],
						 list = T,
						 overwrite = F)$Name

		# list .tif files indexes inside the zips
		toUnzip <- grep(".tif", auxList)

		# extracts .tif files from zips
		if(!prod(auxList[toUnzip] %in% dir(outPath))) {
			unzip(files[i],
				  files = auxList[toUnzip],
				  exdir = outPath)
		}

		# rename unzipped files to match zip name
		zipName <- strsplit(strsplit(files[i], "/")[[1]][2], "\\.")[[1]][1]
		pat <- strsplit(auxList[toUnzip[1]], "\\.")[[1]][1]
		files_pat <- dir(path = outPath, pattern = pat, full.names = T)
		files_new <- gsub(pat, zipName, files_pat)
		file.rename(files_pat, files_new)

		# update progress bar
		setTxtProgressBar(pb, i)
	}
	close(pb)
}

# list all .tif files inside the extracted folder
files <- dir(path = paste0(p_raster, f_area, "/extracted-modis/"),
			 pattern = paste0("ndvi", ".tif"))

# reorder Terra/Aqua files by date
filesDates <- substr(files, 13, 22)
files <- files[order(filesDates)]

# create rasters for all files
rasters <- lapply(paste0(p_raster, f_area, "/extracted-modis/", files),
				  FUN = raster)

# defining the stack of all the rasters
s <- raster::stack(rasters)

z <- velox(s)

# saving stack
writeRaster(
	z,
	filename = paste0(p_down, "eucalipto-full.tif"),
	options = "INTERLEAVE=BAND",
	overwrite = TRUE
)

z$write(paste0(p_shape, "eucalipto-full"))
beep()

# reads shapefile and transforms its crs match rasters crs
shp <- shapefile(paste0(p_shape, f_area, ".shp"))
shp <- spTransform(shp, crs(s[[1]]))

# NOT RUN
# takes ~60min to run
system.time({
	z <- mask(s, shp)
})
beep()
raster::plot(mask(s[[1]], shp))

# saving stack
writeRaster(
	z,
	filename = paste0(p_down, "eucalipto-recorte.tif"),
	options = "INTERLEAVE=BAND",
	overwrite = TRUE
)

# saving the shapefile
writeOGR(shp,
		 dsn = paste0(p_down, "eucalipto"),
		 layer = "eucalipto",
		 driver = "ESRI Shapefile")

# delete temporary directory with unzipped images
unlink(paste0(p_raster, f_area, "/", "extracted-modis"),
	   recursive = T)