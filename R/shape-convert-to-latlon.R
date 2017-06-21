# -------------------------------------------------------------------------(#HEADER)
# Script:
# 		shape-convert-to-latlon
# Author:
# 		Alexandre Esteves Almeida (2016)
# Description:
# 		- Batch convert shapes in any projection to latlon

# -------------------------------------------------------------------------(#IMPDEF)

# load packages, install if necessary
packs <- c("raster",
			  "rgdal")

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

# -------------------------------------------------------------------------(#MAIN)

setwd(p_shape)
files <- grep(".shp", dir())
filesnames <- substr(dir()[files], 1, nchar(dir()[files])-4)

for(i in 1:length(files)) {
	# read
	shp <- shapefile(paste0(filesnames[i], ".shp"))
	if(is.na(crs(shp))) crs(shp) <- proj_utm
	shp <- spTransform(shp, CRS(proj_ll(shp)))

	# save
	shapefile(shp,
				 filename = paste0(filesnames[i], "_ll"),
				 overwrite = T)
}