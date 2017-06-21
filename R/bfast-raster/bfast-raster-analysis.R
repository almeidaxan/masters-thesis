# --------------------------------------------------------- PROCESSING ----

# defining file names to be loaded
f_area <- "fray-jorge-park"
f_band <- "evi2" # "ndvi", "evi", "evi2"
f_sat <- "ls" # "ls" for Landsat, "modis" for MODIS

# loading .Rdata with bfast outputs
load(paste0(p_rdata, scriptName, "/", f_area, "-", f_sat, "_", f_band, ".Rdata"))

# loading the shape of the region
shp <- shapefile(paste0(p_shape, f_area, ".shp"))

# ################################################## PLACEHOLDER TEMPORARIO
# DESCRICAO: pixels que tem bp.S (somente alguns espalhados) causam
# inconsistencias nas analises. SOLUCAO TEMPORARIA: pegar informacao dos
# pixels da vizinhanca enquanto não arrumo isso
z1 <- sapply(bfastOut, function(x) {
	if(length(x$bp.S)==0) {
		NA
	} else {
		if(is.na(x$bp.S[1])) {
			0
		} else {
			length(x$bp.S)
		}
	}
})
z2 <- which(z1 > 0)
while(length(z2 > 0)) {
	for(i in z2) {
		bfastOut[[i]]$St <- bfastOut[[i+1]]$St
		bfastOut[[i]]$bp.S <- bfastOut[[i+1]]$bp.S
	}
	z1 <- sapply(bfastOut, function(x) {
		if(length(x$bp.S)==0) {
			NA
		} else {
			if(is.na(x$bp.S[1])) {
				0
			} else {
				length(x$bp.S)
			}
		}
	})
	z2 <- which(z1 > 0)
}
# ################################################## PLACEHOLDER TEMPORARIO

# 1) trend - number of bps
so1 <- s[[1]]
values(so1) <- sapply(bfastOut, function(x) {
	ifelse(is.null(x$bp.T),
		   NA,
		   ifelse(sum(is.na(x$bp.T) > 0), 0, length(x$bp.T)))
})

# 2) trend - median magnitude of all bps
so2 <- s[[1]]
values(so2) <- sapply(bfastOut, function(x) {
	ifelse(is.null(x$mag.T),
		   NA,
		   ifelse(is.null(dim(x$mag.T)), 0, median(x$mag.T[,3])))
})

# 3) season - periodogram
so3 <- s[[1]]
values(so3) <- sapply(bfastOut, function(x) {
	if(is.na(x$St[1])) {
		NA
	} else {
		a <- x$St
		p <- pgm(a)
		mp <- sort(p$P, TRUE)[1:5]
		# weighted mean of frequencies (months) using peridogram values
		sum((length(a)/as.numeric(names(mp)) * mp)/sum(mp))
	}
})

soz <- so3
values(soz)[values(soz) <= 9] <- 9

# 4) phenology parameters
p <- 0.2

	# 4.1) greening season - start of season
	so4_1 <- s[[1]]
	values(so4_1) <- round(sapply(1:length(bfastOut), function(i) {
		if(is.na(bfastOut[[i]]$St[1]) | is.null(bfastOut[[i]]$St[1])) {
			NA
		} else {
			z <- round(bfastOut[[i]]$St, 5)

			zz <- (z - min(z))
			zz <- zz/max(zz)

			ind_m <- which(z == min(z))[1]
			j <- 1
			ind_M <- which(z == max(z))[j]
			while(ind_M <= ind_m) {
				j <- j + 1
				ind_M <- which(z == max(z))[j]
			}

			a <- as.numeric(zz[ind_m:ind_M])
			b <- time(z)[ind_m:ind_M]
			b <- (b - floor(b))*12 + 1

			ind_a <- which(a >= p)[1]

			b[ind_a]
		}
	}), 0)

	# 4.2) greening season - length
	so4_2 <- s[[1]]
	values(so4_2) <- round(sapply(1:length(bfastOut), function(i) {
		if(is.na(bfastOut[[i]]$St[1]) | is.null(bfastOut[[i]]$St[1])) {
			NA
		} else {
			z <- round(bfastOut[[i]]$St, 5)

			zz <- (z - min(z))
			zz <- zz/max(zz)

			ind_m1 <- which(z == min(z))[1]

			i <- 1
			ind_M <- which(z == max(z))[i]
			while(ind_M <= ind_m1) {
				i <- i + 1
				ind_M <- which(z == max(z))[i]
			}

			j <- 2
			ind_m2 <- which(z == min(z))[j]
			while(ind_m2 <= ind_M) {
				j <- j + 1
				ind_m2 <- which(z == min(z))[j]
			}

			a <- as.numeric(zz[ind_m1:ind_m2])

			sum(!(a < p))
		}
	}), 0)

	# 4.3) greening season - amplitude
	so4_3 <- s[[1]]
	values(so4_3) <- sapply(1:length(bfastOut), function(i) {
		if(is.na(bfastOut[[i]]$St[1]) | is.null(bfastOut[[i]]$St[1])) {
			NA
		} else {
			z <- round(bfastOut[[i]]$St, 5)
			max(z)-min(z)
		}
	})

# -------------------------------------------------------------- PLOTS ----

saveFile <- T

plotRes(so = so1,
		elev = elev2,
		col = brewer.pal(9, "RdPu"),
		legTitle = "Frequency of Breakpoints",
		saveFile = saveFile,
		fileName = "eco-fray-jorge-results-breaks-freq",
		legType = "legend")

plotRes(so = so2,
		elev = elev2,
		col = c("chocolate2", "black", "chartreuse2"),
		legTitle = "Magnitude of Breakpoints",
		saveFile = saveFile,
		fileName = "eco-fray-jorge-results-breaks-mag",
		legType = "colorbar",
		legBias = 0.87)

plotRes(so = soz,
		elev = elev2,
		col = c("cyan", "black"),
		legTitle = "Seasonal Periodogram Frequency",
		saveFile = saveFile,
		fileName = "eco-fray-jorge-results-periodogram",
		legType = "colorbar")

plotRes(so = so4_1,
		elev = elev2,
		col = hsvPalette[-1],
		legTitle = "Start of Seasonal Period",
		saveFile = saveFile,
		fileName = "eco-fray-jorge-results-sos",
		legType = "legend")

plotRes(so = so4_3,
		elev = elev2,
		col = c("black", "green"),
		legTitle = "Seasonal Amplitude",
		saveFile = saveFile,
		fileName = "eco-fray-jorge-results-seasonal-amp",
		legType = "colorbar")

# ------------------------------------------------------ PREC ANALYSIS ----

cruP <- nc_open(paste0(p_nc, "prec_1961_2014_SA.nc"), readunlim = F)

prec <- ncvar_get(cruP, "pre")
lon <- ncvar_get(cruP, "lon")
lat <- ncvar_get(cruP, "lat")
tim <- ncvar_get(cruP, "time")
tim <- strftime(as.Date(tim, origin = "1900-01-01"), "%Y-%m")

ext <- rowMeans(bbox(shapefile(paste0(p_shape, f_area, ".shp"))))

lon_ext <- which.min(abs(lon-ext[1]))
lat_ext <- which.min(abs(lat-ext[2]))

precTs <- prec[lon_ext, lat_ext, 1:dim(prec)[3]]

mYear <- as.numeric(substr(head(tim,1), 1, 4))
MYear <- as.numeric(substr(tail(tim,1), 1, 4))
mMonth <- as.numeric(substr(head(tim,1), 6, 7))
MMonth <- as.numeric(substr(tail(tim,1), 6, 7))

precTs <- ts(precTs,
			  start = c(mYear, mMonth),
			  end = c(MYear, MMonth),
			  frequency = 12)

whichNotNA <- which(sapply(lapply(bfastOut, function(x) {x$Yt[1]}), function(x) !is.na(x)))[1]

precTs <- window(precTs,
				 tsp(bfastOut[[whichNotNA]]$Yt)[1],
				 tsp(bfastOut[[whichNotNA]]$Yt)[2])

serie <- matrix(sapply(bfastOut, function(x) {x$St}),
				ncol = length(bfastOut))

serie <- ts(serie,
			start = tsp(bfastOut[[whichNotNA]]$Yt)[1],
			end = tsp(bfastOut[[whichNotNA]]$Yt)[2],
			frequency = 12)

# 1) maior escala (pixel medio, regiao inteira)
	px1 <- ts(rowMeans(serie, na.rm = T),
			  start = tsp(serie[,1])[1],
			  end = tsp(serie[,1])[2],
			  frequency = 12)

	serie_ccf <- ccf(px1, precTs, plot = F)
	serie_acf <- serie_ccf$acf[serie_ccf$lag>=0]
	serie_lag <- serie_ccf$lag[serie_ccf$lag>=0]

	# maior correlacao
	serie_1 <- max(na.omit(abs(serie_acf)))
	serie_1

	# lag com maior correlacao
	serie_2 <- serie_lag[which.max(abs(serie_acf))]*12
	serie_2

# 2) menor escala (pixel by pixel)
	px2_1 <- s[[1]]
	px2_2 <- s[[1]]
	serie_1 <- NULL
	serie_2 <- NULL
	for(i in 1:ncol(serie)) {
		if(is.na(serie[1,i])) {
			serie_1[i] <- serie_2[i] <- NA
		} else {
			serie_ccf <- ccf(serie[,i], precTs, plot = F, lag.max = 5)

			serie_acf <- serie_ccf$acf[serie_ccf$lag>=0]
			serie_lag <- serie_ccf$lag[serie_ccf$lag>=0]

			if(max(serie_acf) <= qnorm(0.95+0.05/2)/sqrt(serie_ccf$n.used)) {
				serie_1[i] <- 0
				serie_2[i] <- NA
			} else {
				serie_1[i] <- max(serie_acf)
				serie_2[i] <- serie_lag[which.max(serie_acf)]*12
			}
		}
	}
	values(px2_1) <- serie_1
	values(px2_2) <- serie_2

	plotRes(so = px2_1,
			elev = elev2,
			col = c("black", "cyan"),
			legTitle = "Highest Correlation",
			saveFile = saveFile,
			fileName = "eco-fray-jorge-results-prec-cor",
			legType = "colorbar")

	plotRes(so = px2_2,
			elev = elev2,
			col = rev(c("green", "black")),
			legTitle = "Lag",
			saveFile = saveFile,
			fileName = "eco-fray-jorge-results-prec-lag",
			legType = "legend",
			legNCol = 6)

	# plotRes(so = px2_2 + 6,
	# 		col = hsvPalette[6:11],
	# 		legTitle = "Month",
	# 		saveFile = saveFile,
	# 		legCycle = T,
	# 		fileName = paste0(f_area, "prec_raster"),
	# 		legType = "legend")

# monthly boxplot
	precMon <- matrix(c(NA, precTs, NA), ncol = 12, byrow = T)
	precMon <- data.frame(c(precMon), rep(1:12, each = nrow(precMon)))
	colnames(precMon) <- c("prec", "mon")
	precMon <- precMon[-which(is.na(precMon$prec)),]

	# qual o mes cuja precipitacao acumulada media atinge pelo menos 20%?
	precMonCum <- aggregate(precMon$prec, list(precMon$mon), mean)
	precMonCum <- cumsum(precMonCum$x)/max(cumsum(precMonCum$x))
	which(precMonCum > .2)[1] # mes 6

	precPalette <- rep("lightgray", 12)
	precPalette[which(precMonCum > .2)[1]] <- "#73b1ef"

	gg <- ggplot(data = precMon, mapping = aes(x = mon, y = prec, group = mon)) +
		stat_boxplot(geom = "errorbar") +
		geom_boxplot(fill = precPalette,
					 outlier.shape = NA) +
		labs(x = "Month",
			 y = "Precipitation (mm)"
			 # title = "Fray Jorge\nMonthly precipitation from 1990 to 2011"
			 ) +
		stat_summary(
			fun.y = mean,
			geom = "point",
			shape = 12,
			size = 2
		) +
		scale_x_continuous(breaks = 1:12) +
		scale_y_continuous(limits = c(0,100)) +
		theme_light() +
		theme(
			text = element_text(size = 14),
			plot.title = element_text(size = 14),
			axis.text = element_text(size = 12)
		)
	gg
	ggsave(
		file = paste0(p_down, "prec_boxplot.pdf"),
		plot = gg,
		width = 10,
		height = 7
	)

# ------------------------------------------------------- FOG ANALYSIS ----

# filter only zip files from wd
setwd(paste0(p_raster, f_area, "/MOD09GA"))
files <- dir(pattern = ".zip", full.names = T)

# remove .tmp files if they exist
tmp <- dir(pattern = ".tmp", full.names = T)
tmp <- file.remove(tmp)

# output folder
dir.create("../extracted-fog/", showWarnings = F)
outPath <- paste0("../extracted-fog")

# progress bar to follow up loop progress
print("Extracting files...")
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

# list all .tif files inside the extracted folder
filesMOD09 <- dir(path = paste0(p_raster, f_area, "/extracted-fog/"),
			 pattern = paste0("MOD09.tif"))
filesMOD35 <- dir(path = paste0(p_raster, f_area, "/extracted-fog/"),
				  pattern = paste0("MOD35.tif"))

# create rasters for all files
rastersMOD09 <- lapply(paste0(p_raster, f_area, "/extracted-fog/", filesMOD09),
				  FUN = raster)
rastersMOD35 <- lapply(paste0(p_raster, f_area, "/extracted-fog/", filesMOD35),
					   FUN = raster)

# create a stack of rasters
fogStack1 <- raster::stack(rastersMOD09)
fogStack2 <- raster::stack(rastersMOD35)

# create a time variable (Z) for both fogStacks
fogStack1 <- setZ(fogStack1, as.Date(paste0(substr(names(fogStack1), 2, 8), ".01"), format = "%Y.%m.%d"))
fogStack2 <- setZ(fogStack2, as.Date(paste0(substr(names(fogStack2), 2, 8), ".01"), format = "%Y.%m.%d"))

# extract time series from fogStacks
dFog1 <- zooExtract(x = fogStack1, 1:ncell(fogStack1))
dFog2 <- zooExtract(x = fogStack2, 1:ncell(fogStack2))

# convert raster pixels to polygons
fogPolys <- rasterToPolygons(projectRaster(fogStack1[[1]], crs = proj_utm(fogStack1)))

# create an list of extents of all polygons
fogCoords <- list()
for(i in 1:nrow(fogPolys@data)) {
	fogCoords[[i]] <- fogPolys@polygons[[i]]@Polygons[[1]]@coords[c(1,3),]
}

# find all corresponding pixels from the ROI in fog
z <- projectRaster(s[[1]])
values(z) <- 0
sPoints <- rasterToPoints(z)[,1:3]
for(i in 1:nrow(sPoints)) {
	tmp <- which(sapply(fogCoords, function(x) {
		sPoints[i,1] >= x[1,1] & sPoints[i,1] <= x[2,1] & sPoints[i,2] >= x[2,2] & sPoints[i,2] <= x[1,2]
	}))
	if(length(tmp) == 0) sPoints[i,3] <- NA
	else sPoints[i,3] <- tmp
}
colnames(sPoints)[3] <- "fogPixel"

whichNotNA <- which(sapply(lapply(bfastOut, function(x) {x$Yt[1]}), function(x) !is.na(x)))[1]

# analise na menor excala (pixel by pixel) da correlacao com o fog
serie <- matrix(sapply(bfastOut, function(x) { x$St }),
				ncol = length(bfastOut))

serie <- ts(serie,
			start = tsp(bfastOut[[whichNotNA]]$Yt)[1],
			end = tsp(bfastOut[[whichNotNA]]$Yt)[2],
			frequency = 12)

fogTs1 <- ts(dFog1,
			start = c(as.numeric(substr(head(getZ(fogStack1), 1), 1, 4)),
					  as.numeric(substr(head(getZ(fogStack1), 1), 6, 7))),
			end = c(as.numeric(substr(tail(getZ(fogStack1), 1), 1, 4)),
					  as.numeric(substr(tail(getZ(fogStack1), 1), 6, 7))),
			frequency = 12)

fogTs2 <- ts(dFog2,
			 start = c(as.numeric(substr(head(getZ(fogStack1), 1), 1, 4)),
			 		  as.numeric(substr(head(getZ(fogStack1), 1), 6, 7))),
			 end = c(as.numeric(substr(tail(getZ(fogStack1), 1), 1, 4)),
			 		as.numeric(substr(tail(getZ(fogStack1), 1), 6, 7))),
			 frequency = 12)

# match serie to fog start period
serie <- window(serie,
				start = c(as.numeric(substr(head(getZ(fogStack1), 1), 1, 4)),
						  as.numeric(substr(head(getZ(fogStack1), 1), 6, 7))))

# match fog to serie end period
fogTs1 <- window(fogTs1,
				 end = tsp(serie)[2])

fogTs2 <- window(fogTs2,
				 end = tsp(serie)[2])

fx1 <- s[[1]]
fx2 <- s[[1]]
serie_1 <- NULL
serie_2 <- NULL
for(i in 1:ncol(serie)) {
	if(is.na(serie[1,i])) {
		serie_1[i] <- serie_2[i] <- NA
	} else {
		serie_ccf <- ccf(serie[,i], fogTs2[,sPoints[i, 3]], plot = F, lag.max = 6)

		serie_acf <- serie_ccf$acf[serie_ccf$lag >= 0]
		serie_lag <- serie_ccf$lag[serie_ccf$lag >= 0]

		if(max(serie_acf) <= qnorm(0.95 + 0.05 / 2)/sqrt(serie_ccf$n.used)) {
			serie_1[i] <- 0
			serie_2[i] <- NA
		} else {
			serie_1[i] <- max(serie_acf)
			serie_2[i] <- serie_lag[which.max(serie_acf)]*12
		}
	}
}
values(fx1) <- serie_1
values(fx2) <- serie_2

plotRes(so = fx1,
		elev = elev2,
		col = c("black", "cyan"),
		legTitle = "Highest Correlation",
		saveFile = saveFile,
		fileName = "eco-fray-jorge-results-fog-cor",
		legType = "colorbar")

plotRes(so = fx2,
		elev = elev2,
		col = c("black", "green"),
		legTitle = "Lag",
		saveFile = saveFile,
		fileName = "eco-fray-jorge-results-fog-lag",
		legNCol = 7,
		legType = "legend")

# delete temporary directory with unzipped images
unlink(paste0(p_raster, f_area, "/extracted-fog/"),
	   recursive = T)

# monthly boxplot
fogMon <- matrix(c(NA, fogTs2[,23]), ncol = 12, byrow = T)
fogMon <- data.frame(c(fogMon), rep(1:12, each = nrow(fogMon)))
colnames(fogMon) <- c("fog", "mon")
fogMon <- fogMon[-which(is.na(fogMon$fog)),]

# qual o mes cuja precipitacao acumulada media atinge pelo menos 20%?
fogMonCum <- aggregate(fogMon$fog, list(fogMon$mon), mean)
fogMonCum <- cumsum(fogMonCum$x)/max(cumsum(fogMonCum$x))
z[i] <- which(fogMonCum > .2)[1]

fogPalette <- rep("white", 12)
fogPalette[10:11] <- "royalblue3"

gg <- ggplot(data = fogMon, mapping = aes(x = mon, y = fog, group = mon)) +
	stat_boxplot(geom = "errorbar") +
	geom_boxplot(fill = fogPalette) +
	labs(x = "Month",
		 y = "Fog",
		 title = "Fray Jorge\nMonthly fog from 1990 to 2011\nFog pixel with highest rain correlation at months 10/11") +
	stat_summary(
		fun.y = mean,
		geom = "point",
		shape = 6,
		size = 2
	) +
	scale_x_continuous(breaks = 1:12) +
	scale_y_continuous(limits = c(0,100)) +
	theme_dark() +
	theme(
		text = element_text(size = 16),
		plot.title = element_text(size = 16),
		axis.text = element_text(size = 12)
	)
gg
ggsave(
	file = paste0(p_down, "fog_boxplot_3.pdf"),
	device = "pdf",
	plot = gg,
	width = 12,
	height = 9,
	pointsize = 10
)


# ---------------------------------------------------------- ELEVATION ----

# unzipping tif from elevation zip
outPath <- paste0(p_raster, f_area, "/SRTMGL1_003/")
auxList <- unzip(paste0(outPath, f_area, "_elevation.zip"),
				 list = T,
				 overwrite = F)$Name

# list .tif files indexes inside the zips
toUnzip <- grep(".tif", auxList)

# extracts .tif files from zips
if(!prod(auxList[toUnzip] %in% dir(outPath))) {
	unzip(paste0(outPath, f_area, "_elevation.zip"),
		  files = auxList[toUnzip],
		  exdir = outPath)
}

elev <- raster(paste0(outPath, "SRTMGL1_003.elevation.tif"))
elev2 <- mask(elev, shp)
slope <- terrain(elev, opt = "slope")
aspect <- terrain(elev, opt = "aspect")
hill <- hillShade(slope, aspect, 45, 120)

plotRes(so = mask(elev, shp),
		col = rev(c("white", "wheat4", "black", "olivedrab", "olivedrab1")),
		legTitle = "Altitude (m)",
		saveFile = saveFile,
		fileName = paste0(f_area, "_elev"),
		legStep = 50)

############# TESTS

### FLOODING

z <- mask(elev, shp)
z[values(z) > 300] <- NA
z[!is.na(values(z))] <- 1
raster::plot(
	z,
	col = "blue",
	legend = FALSE,
	add = T
)

### ALTITUDE SLOPE/ASPECT

saveGIF(movie.name = paste0(p_down, "teste.gif"),
		interval = 0.2,
		for (i in seq(0, 360, 45)) {
			hill <- hillShade(slope, aspect, 45, i)

			# fray jorge
			plotRes(
				so = projectRaster(so2, crs = proj_ll),
				col = c("chocolate2", "black", "chartreuse2"),
				legTitle = f_band,
				saveFile = saveFile,
				fileName = paste0(f_area, "_", f_sat, "_trend_bpmag"),
				legStep = 0.02,
				legBias = 0.87,
				legRound = 2
			)

			raster::plot(
				mask(hill, shp),
				col = grey(c(0:35, 65:100) / 100, alpha = 0.4),
				legend = FALSE,
				add = T
			)
		})

##############

so <- so4_3

# first let both rasters have the same projection, the do a resampling
z <- resample(projectRaster(so, crs = proj_ll, method = "ngb"), elev, method = "ngb")

# par(mfrow = c(1,2))
#
# # visualizing both rasters, after reproject and resampling
# plotRes(so = z,
# 		col = c("black", "green"),
# 		legTitle = f_band,
# 		saveFile = saveFile,
# 		fileName = paste0(f_area, "_", f_sat, "_season_amp"),
# 		legStep = 0.02,
# 		legRound = 2)
# raster::plot(
# 	elev,
# 	col = grey(c(0:35, 65:100) / 100, alpha = 0.4),
# 	legend = FALSE,
# 	add = T
# )

# par(mar = c(5, 4, 4, 2) + 0.1)

# finding correspondence between rasters pixels
p1 <- rasterToPoints(z)
p2 <- rasterToPoints(elev)
m <- merge(p1, p2)

plot(m[, 3], m[, 4], pch = 20, col = rgb(0, 0, 0, 0.1), xlab = "Seasonal amplitude", ylab = "Elevation")
dev.off()
boxplot(
	m[, 4] ~ m[, 3],
	pch = 20,
	col = rgb(0, 0, 0, 0.1),
	ylab = "Elevation",
	xlab = "Start of season (month)"
)

# plot(m[, 3], floor(m[, 4]/100), pch = 20, col = rgb(0, 0, 0, 0.1), xlab = "medida", ylab = "elevacao")
# boxplot(
# 	m[, 3] ~ floor(m[, 4] / 100),
# 	pch = 20,
# 	col = rgb(0, 0, 0, 0.1),
# 	ylab = "medida",
# 	xlab = "elevacao"
# )
