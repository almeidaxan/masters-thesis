# -------------------------------------------------------------------------
# Author:
# 		Alexandre Esteves Almeida (2016)
# ------------------------------------------------------------- IMPDEF ----

scriptName <- "bfast-csv"

# load packages, install if necessary
packs <- c("bfastSpatial",
		   "raster",
		   "lubridate",
		   "ncdf4")

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
