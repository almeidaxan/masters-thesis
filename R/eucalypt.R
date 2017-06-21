# ------------------------------------------------------------- IMPDEF ----

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

# ----------------------------------------------------------- SAMPLING ----

ras <- raster(paste0(p_raster, "eucalipto/MOD13Q1/MOD13Q1_005_2012_08_28/d69e24f856c9dbfeb279a459ad134250.ndvi.tif"))
shpArea <- paste0(p_shape, "emprise_classification2014.shp") %>% shapefile()
shpArea <- spTransform(shpArea, proj_ll)
shpEuc <- paste0(p_shape, "eucalipto.shp") %>% shapefile()
shpOther <- gDifference(shpArea, shpEuc)

# estimando a area dos poligonos
areaEuc <- sapply(shpEuc@polygons[[1]]@Polygons, function(x) {
	x@area
})

# escolhendo os 150 poligonos com maior area para evitar false positives
top <- order(areaEuc, decreasing = T)[1:150]
shpEucTop <- shpEuc
shpEucTop@polygons[[1]]@Polygons <- shpEucTop@polygons[[1]]@Polygons[top]

# mask raster from the region with Euc/Other shapes
suppressWarnings({
	rasEuc <- mask(ras, shpEucTop, updatevalue = -1, updateNA = -1)
	rasOther <- mask(ras, shpOther, updatevalue = -1, updateNA = -1)
})
pointsEuc <- which(values(rasEuc) != -1, arr.ind = T)
pointsOther <- which(values(rasOther) != -1, arr.ind = T)

# sampling Euc/Other pixels
set.seed(1)
smpEuc <- sample(pointsEuc, 250)
smpOther <- sample(pointsOther, 1000)

# defining file names to be loaded
f_area <- "eucalipto-guerric"
f_band <- "ndvi" # "ndvi", "evi", "evi2"
f_dir <- dir(paste0(p_raster, f_area))
f_sat <- which(startsWith(f_dir, "M")) # "L" for Landsat; "M" for MODIS

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
			 pattern = paste0(f_band, ".tif"))

# reorder Terra/Aqua files by date
filesDates <- substr(files, 13, 22)
files <- files[order(filesDates)]

# create rasters for all files
rasters <- lapply(paste0(p_raster, f_area, "/extracted-modis/", files),
				  FUN = raster)

# defining the stack of all the rasters
s <- raster::stack(rasters)

# if s has no time variable, create one
if(is.null(getZ(s))) {
	s <- setZ(s, as.Date(substr(names(s), 13, 22), format = "%Y_%m_%d"))
}

# create and extract time series from all the points (or a sample)
d <- zooExtract(x = s, c(smpEuc, smpOther))
d <- as.matrix(d)

# reduce date to months
dtime <- strftime(getZ(s), "%Y-%m")

# aggregate data by month
l <- length(unique(dtime))
d2 <- matrix(nrow = l, ncol = ncol(d))
for(i in 1:l) {
	w <- which(dtime == unique(dtime)[i])
	if(length(w) == 1) d2[i,] <- d[w,]
	else d2[i,] <- colMedians(d[w,], na.rm = T)
}
d2[which(is.nan(d2))] <- NA
rm(d)

# create the full year-month empty data structure
mYear <- as.numeric(substr(min(getZ(s)), 1, 4))
MYear <- as.numeric(substr(max(getZ(s)), 1, 4))
mMonth <- as.numeric(substr(min(getZ(s)), 6, 7))
MMonth <- as.numeric(substr(max(getZ(s)), 6, 7))
fullYear <- rep(mYear:MYear, each = 12)
fullMonth <- rep(1:12, times = length(mYear:MYear))
if(mMonth>1) {
	fullYear <- fullYear[-(1:(mMonth-1))]
	fullMonth <- fullMonth[-(1:(mMonth-1))]
}
fullYear <- fullYear[-((length(fullYear)-12+MMonth+1):length(fullYear))]
fullMonth <- fullMonth[-((length(fullMonth)-12+MMonth+1):length(fullMonth))]
dfull <- data.frame(time = paste0(fullYear,"-",sprintf("%02d", fullMonth)))
timediff <- setdiff(as.character(dfull$time), unique(dtime))

# merging the aggregated data with full year-month
d3 <- rbind(d2, matrix(nrow = length(timediff), ncol = ncol(d2)))
d3 <- data.frame(time = c(unique(dtime), timediff), d3)
d3 <- d3[order(d3$time),]
d3$time <- NULL
colnames(d3) <- NULL
rm(d2)

dts1 <- ts(
	d3,
	start = c(mYear, mMonth),
	end = c(MYear, MMonth),
	frequency = 12
)
rm(d3)

# filling NA gaps with linear interpolation
dts <- na.approx(dts1, rule = 2)
dts[which(is.na(dts))] <- 0
rm(dts1)

# -------------------------------------------------------------- BFAST ----

# define number of cores for parallel programming, based on OS
# detectCores()-1 cores will be used
argList <- list()
if(Sys.info()[["sysname"]] == "Windows") {
	# windows OS
	cl <- makeCluster(detectCores()-1)
	registerDoSNOW(cl)

	# progressbar (only doSNOW package is supported)
	pb <- txtProgressBar(min = 0,
						 max = ncol(dts),
						 style = 3)
	progressBar <- function(n) setTxtProgressBar(pb, n)
	argList$.options.snow <- list(progress = progressBar)
} else {
	# unix-based OSs
	registerDoMC(cores = detectCores()-1)
}

# minimum segment size (months) used to calculate h parameter in bfast
hMon <- 48
h <- hMon/nrow(dts)

# run bfast for all other pixels
system.time(bfastOut <- foreach(i = 1:ncol(dts), .packages = c("bfast")) %dopar% {
	# only runs bfast for points inside shapefile (outside = 0)
	if(sum(dts[,i]) != 0) {
		resBfast <- bfast(Yt = dts[,i],
					 max.iter = 1,
					 h = h)

		# select which outputs from bfast to keep
		out <- list()

		## 1) original time series (Yt)
		out$Yt <- resBfast$Yt
		## 2) trend component (Tt)
		out$Tt <- resBfast$output[[1]]$Tt
		## 3) seasonal component (St)
		out$St <- resBfast$output[[1]]$St
		## 4) trend breakpoints (bp.T)
		if(!resBfast$nobp$Vt) {
			out$bp.T <- resBfast$output[[1]]$bp.Vt$breakpoints
		} else {
			out$bp.T <- NA
		}
		## 5) seasonal breakpoints (bp.S)
		if(!resBfast$nobp$Wt) {
			out$bp.S <- resBfast$output[[1]]$bp.Wt$breakpoints
		} else {
			out$bp.S <- NA
		}
		## 6) trend bp magnitudes matrix (mag.T)
		# Columns described below:
		# [,1] -> value before bp
		# [,2] -> value after bp
		# [,3] -> magnitude of the bp
		out$mag.T <- resBfast$Mags
		return(out)
	} else {
		out <- list(Yt = rep(NA, nrow(dts)),
					Tt = rep(NA, nrow(dts)),
					St = rep(NA, nrow(dts)),
					bp.T = NULL,
					bp.S = NULL,
					mag.T = NULL)
	}

	# update progress bar
	if(Sys.info()[["sysname"]] == "Windows") {
		setTxtProgressBar(pb, i)
	}

	return(out)
})

if(Sys.info()[["sysname"]] == "Windows") {
	close(pb)
	stopCluster(cl)
}

beep()

# ------------------------------------------------------------- SAVING ----

# create rdata dir
dir.create(paste0(p_rdata, scriptName), showWarnings = F)

# save rdata with environment variables
save(file = paste0(p_rdata, scriptName, "/", f_area, "-modis_", f_band, "1000.RData"),
	 list = c("bfastOut", "s", "smpEuc", "smpOther"),
	 envir = .GlobalEnv)

# delete temporary directory with unzipped images
unlink(paste0(p_raster, f_area, "/", "extracted-modis/"),
	   recursive = T)

# --------------------------------------------------------------- ETC. ----

# selecting appropriate h parameter for bfast
for(j in 256:260) {
	dir.create(paste0("bfast/", j), showWarnings = F)
	for(hMon in seq(1, 90, by = 1)) {
		h <- hMon/nrow(dts)

		resBfast <- bfast(Yt = dts[,j],
						  max.iter = 1,
						  season = "harmonic",
						  h = h)

		png(paste0(p_down, "/bfast/", j, "/bfout_", hMon, ".png"), height = 500, width = 700, res = 100)
			plot(dts[,j], ylab = "NDVI", main = paste("hMon = ", hMon))
			lines(resBfast$output[[1]]$Tt + resBfast$output[[1]]$St, col = "blue")
			if(!is.na(resBfast$output[[1]]$ci.Vt[1])) {
				z <- resBfast$output[[1]]$ci.Vt$confint
				if(sum(z < 0))
					z[which(z < 0)] <- 1
				if(sum(z > length(dts[,j])))
					z[which(z > length(dts[,j]))] <- length(dts[,j])
				bps <- time(dts[,j])[z] %>% matrix(nrow = nrow(z), ncol = 3)
				abline(v = bps[,2], col = "red", lty = 2)
			}
		dev.off()
	}
}
