# -------------------------------------------------------------------------(#HEADER)
# Script:
# 		shape-create
# Author:
# 		Alexandre Esteves Almeida (2016)
# Description:
# 		- Create a custom square shape centered in an input latlon point

# -------------------------------------------------------------------------(#IMPDEF)

# load packages, install if necessary
packs <- c("optparse",
		   "raster",
		   "rgdal")

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

# -------------------------------------------------------------------------(#MAIN)

# # define options to run this script
# option_list = list(
# 	make_option(c("--lat"), type="double", default=NULL, metavar="DOUBLE",
# 					help="Decimal latitude of rectangle shape center."),
# 	make_option(c("--lon"), type="double", default=NULL, metavar="DOUBLE",
# 					help="Decimal longitude of rectangle shape center."),
# 	make_option(c("--wid"), type="integer", default=NULL,
# 					metavar="INTEGER", help="Rectangle shape width, in meters."),
# 	make_option(c("--hei"), type="integer", default=NULL,
# 					metavar="INTEGER", help="Rectangle shape height, in meters."),
# 	make_option(c("--out"), type="character", default=NULL,
# 					metavar="STRING", help="Output shape filename.")
# )
#
# opt_parser = OptionParser(option_list = option_list)
# opt = parse_args(opt_parser)

opt <- list()

opt$lat <- -4.0827 # -3.4249 -4.0827
opt$lon <- -41.7083 # -72.8465 -41.7083
opt$wid <- opt$hei <- 15000
opt$out <- "savfor-B"

# all options must be provided, otherwise gives an error message
if (is.null(opt$lat) | is.null(opt$lon) | is.null(opt$wid) | is.null(opt$hei) | is.null(opt$out)){
	print_help(opt_parser)
	stop("Please specify all options: 'lat', 'lon', 'wid', 'hei' and 'out'.\n", call.=FALSE)
}

# define shape center point in latlon and convert to utm
center <- SpatialPoints(coords = data.frame(opt$lon, opt$lat),
								proj4string = CRS(proj_ll))
proj_utm_def <- proj_utm(center)
center <- spTransform(center, CRS(proj_utm_def))

# keep track of coordintes converted to utm
xm <- extent(center)[1] # lon
ym <- extent(center)[3] # lat

# create polygon centered in center
xypoly <- data.frame(x = c(xm-opt$wid/2, xm+opt$wid/2, xm+opt$wid/2, xm-opt$wid/2),
					 y = c(ym-opt$hei/2, ym-opt$hei/2, ym+opt$hei/2, ym+opt$hei/2))

# finally create the desired shape and reconvert to latlon proj
outShape <- SpatialPolygons(list(Polygons(list(Polygon(xypoly)),1)),
									 proj4string = CRS(proj_utm_def))
outShape <- spTransform(outShape, CRS(proj_ll))
outShape <- SpatialPolygonsDataFrame(outShape, data.frame("ID"=0))

# write shape to file
shapefile(outShape,
		  filename = paste0(p_shape, opt$out, ".shp"),
		  overwrite = T)
