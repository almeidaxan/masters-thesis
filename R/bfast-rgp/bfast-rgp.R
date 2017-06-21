# -------------------------------------------------------- GP ANALYSIS ----

# define rdata file to be loaded
f_script <- "bfast-csv"
f_region <- "dantas" # "dantas" "eucalipto"
f_sat <- "MCD13Q1" # "LS5e7e8_SR" "MCD13Q1" "modis"
load(paste0(p_rdata, f_script, "/", f_region, "-", f_sat, ".RData"))

# Xtable ----
out <- list()
# for(j in 1:2) {
	j <- 1
	out[[j]] <- matrix(nrow = length(tsDistances$acronym) + 4, ncol = 1)
	# if(j == 1) {
	# 	f_sat <- "MCD13Q1"
	# } else {
	# 	f_sat <- "LS5e7e8_SR"
	# }
	iCol <- 1
	for(i in c(250, 500, 1000)) {
		# if(i%%3 == 1) {
		# 	viName <- "ndvi"
		# }
		# if(i%%3 == 2) {
		# 	viName <- "evi"
		# }
		# if(i%%3 == 0) {
		# 	viName <- "evi2"
		# }

		load(paste0(p_rdata, scriptName, "/", f_region, "-", f_sat, "-", i, "-", viName, "-2.RData"))

		baseAcc <- sapply(gpRes, function(x) {
			sapply(x, function(y) {
				y$baseAcc
			}) %>% rowMeans()
		})
		decAcc <- sapply(gpRes, function(x) {
			sapply(x, function(y) {
				y$decAcc
			}) %>% rowMeans()
		})
		rawAcc <- sapply(gpRes, function(x) {
			sapply(x, function(y) {
				y$rawAcc
			}) %>% rowMeans()
		})
		resMean <- c(rowMeans(decAcc), rowMeans(rawAcc), rowMeans(baseAcc))
		resSd <- c(apply(decAcc, MARGIN = 1, sd), apply(rawAcc, MARGIN = 1, sd), apply(baseAcc, MARGIN = 1, sd))
		resMean[1:4] <- resMean[1:4] + 2*resSd[1:4]
		resMean[which(resMean > 1)] <- 1
		if(iCol == 1) {
			out[[j]] <- paste(sprintf("%.3f", round(resMean, 3)), "$//pm$", sprintf("%.3f", round(resSd, 3)))
		} else {
			out[[j]] <- cbind(out[[j]], paste(sprintf("%.3f", round(resMean, 3)), "$//pm$", sprintf("%.3f", round(resSd, 3))))
		}

		iCol <- iCol + 1
	}

	rownames(out[[j]]) <- c(
		"GP//textsubscript{dec}",
		"GP//textsubscript{dec.val}",
		"GP//textsubscript{raw}",
		"GP//textsubscript{raw.val}",
		tsDistances$acronym[which(tsDistances$type %in% c("eq", "un"))] %>% as.character
	)

	colnames(out[[j]]) <- c(
		"250",
		"500",
		"1000"
		# "NDVI"
		# "EVI",
		# "EVI2"
	)
# }
print(xtable(out[[1]]), sanitize.text.function = identity)
# print(xtable(out[[2]]), sanitize.text.function = identity)

# Syntatic trees ----
viName <- "ndvi"
load(paste0(p_rdata, scriptName, "/", f_region, "-", f_sat, "-1000-", viName, "-2.RData"))
pdf(file = paste0(p_down, "syntree.pdf"), width = 12, height = 8, pointsize = 14)
	gpPlotSynTree(res = gpRes[[1]][[3]], type = "dec", val = T)
dev.off()

# Importance of distance measures ----
rawDists <- as.character(tsDistances$acronym)
decDists <- c(
	paste0(tsDistances$acronym, "_T"),
	paste0(tsDistances$acronym, "_S"),
	paste0("DEC", 1:3)
)
if(f_region == "eucalipto") {
	viName <- "ndvi"
	for(smpSize in c(250, 500, 1000)) {
		load(paste0(p_rdata, scriptName, "/", f_region, "-", f_sat, "-", smpSize, "-", viName, "-2.RData"))
		a1 <- a2 <- a3 <- a4 <- NULL
		for(seed in 1:5) {
			for(fold in 1:5) {
				result <- gpRes[[seed]][[fold]]

				b1 <- body(result$decBestIndVal)
				b2 <- body(result$decBestIndTrain)
				b3 <- body(result$rawBestIndVal)
				b4 <- body(result$rawBestIndTrain)

				for(i in 1:4) {
					assign(paste0("a", i), {
						c(get(paste0("a", i)), exprToGraph(get(paste0("b", i)))$vertices[exprToGraph(get(paste0("b", i)))$vertices %in% c(rawDists, decDists)])
					})
				}
			}
		}
		for(i in 1:4) {
			tab <- prop.table(table(get(paste0("a", i))))

			if(i %in% c(1,2)) {
				pad <- which(!(decDists %in% (table(get(paste0("a", i))) %>% names())))
				padVec <- rep(0, pad %>% length())
				names(padVec) <- as.character(decDists[pad])
			} else {
				pad <- which(!(rawDists %in% (table(get(paste0("a", i))) %>% names())))
				padVec <- rep(0, pad %>% length())
				names(padVec) <- as.character(rawDists[pad])
			}

			tab <- c(tab, padVec)
			tab <- tab[order(names(tab))]

			assign(paste0("h", i, "_", smpSize), tab)
		}
	}

	outRaw <- cbind(h3_250, h4_250, h3_500, h4_500, h3_1000, h4_1000) %>% multiply_by(100) %>% round(1)
	outDec <- cbind(h1_250, h2_250, h1_500, h2_500, h1_1000, h2_1000) %>% multiply_by(100) %>% round(1)

	colnames(outDec) <- colnames(outRaw) <- rep(c("Train", "Val."), 3)

	print(xtable(outRaw[,seq(1,5,2)], digits = 1), sanitize.text.function = identity)
	print(xtable(outDec[,seq(1,5,2)], digits = 1), sanitize.text.function = identity)

	apply(outRaw[,seq(1,5,2)], MARGIN = 2, FUN = {function(x) {sort(x, decreasing = T) %>% head(5)}})
	apply(outDec[,seq(1,5,2)], MARGIN = 2, FUN = {function(x) {sort(x, decreasing = T) %>% head(5)}})
}
if(f_region == "dantas") {
	outRaw <- outDec <- NULL
	for(f_sat in c("LS5e7e8_SR", "MCD13Q1")) {
		for(viName in names(bfastOut)) {
			load(paste0(p_rdata, scriptName, "/", f_region, "-", f_sat, "-", viName, "-2.RData"))
			a1 <- a2 <- a3 <- a4 <- NULL
			for(seed in 1:10) {
				for(fold in 1:5) {
					result <- gpRes[[seed]][[fold]]

					b1 <- body(result$decBestIndVal)
					b2 <- body(result$decBestIndTrain)
					b3 <- body(result$rawBestIndVal)
					b4 <- body(result$rawBestIndTrain)

					for(i in 1:4) {
						assign(paste0("a", i), {
							c(get(paste0("a", i)), exprToGraph(get(paste0("b", i)))$vertices[exprToGraph(get(paste0("b", i)))$vertices %in% c(rawDists[which(tsDistances$type == "un")], decDists[c(which(tsDistances$type == "un"), which(tsDistances$type == "un") + 20, c(41, 42, 43))])])
						})
					}
				}
			}
			for(i in 1:4) {
				tab <- prop.table(table(get(paste0("a", i))))

				if(i %in% c(1,2)) {
					pad <- which(!(decDists[c(which(tsDistances$type == "un"), which(tsDistances$type == "un")+20)] %in% (table(get(paste0("a", i))) %>% names())))
					padVec <- rep(0, pad %>% length())
					names(padVec) <- as.character(decDists[pad])
				} else {
					pad <- which(!(rawDists[which(tsDistances$type == "un")] %in% (table(get(paste0("a", i))) %>% names())))
					padVec <- rep(0, pad %>% length())
					names(padVec) <- as.character(rawDists[pad])
				}

				tab <- c(tab, padVec)
				tab <- tab[order(names(tab))]

				assign(paste0("h", i, "_", viName), tab)
			}
		}
		outRaw <- cbind(outRaw, h3_ndvi, h4_ndvi, h3_evi, h4_evi, h3_evi2, h4_evi2)
		outDec <- cbind(outDec, h1_ndvi, h2_ndvi, h1_evi, h2_evi, h1_evi2, h2_evi2)
	}

	outRaw <- outRaw %>% multiply_by(100) %>% round(1)
	outDec <- outDec %>% multiply_by(100) %>% round(1)

	colnames(outRaw) <- colnames(outDec) <- rep(c("Train", "Val."), 6)

	print(xtable(outRaw[,seq(1,11,2)], digits = 1), sanitize.text.function = identity)
	print(xtable(outDec[,seq(1,11,2)], digits = 1), sanitize.text.function = identity)

	apply(outRaw[,seq(1,11,2)], MARGIN = 2, FUN = {function(x) {sort(x, decreasing = T) %>% head(5)}})
	apply(outDec[,seq(1,11,2)], MARGIN = 2, FUN = {function(x) {sort(x, decreasing = T) %>% head(5)}})
}

# Fitness vs time ----
tmp <- {function(so) {
	tmp2 <- unlist(lapply(so, function(x) {
		c(max(x$fitnessEvo[,1]),
		  min(x$fitnessEvo[,2]),
		  max(x$fitnessEvo[,2]))
	}))
	tmp2 <- matrix(unlist(tmp2), ncol = 3, byrow = T)
	tmp2
}}
tmp2 <- rbind(tmp(z1), tmp(z2), tmp(z3))
tmp3 <- {function(so) {
	lapply(so, function(x) {
		x$fitnessEvo
	})
}}

# Leaflet map ----
	# savfor
	indA <- 306 # landsat 281, modis 306
	indB <- 78 # landsat 70, modis 78

	z1 <- sapply(bfastOut[[1]], function(x){ c(x$coords) }) %>% t()
	z2 <- sapply(bfastOut[[1]], function(x){ x$class })
	z2 <- sapply(bfastOut[[1]], function(x){ x$class })
	z3 <- sapply(bfastOut[[1]], function(x){ ifelse(is.na(x$bp.T[1]), 0, x$bp.T %>% length()) })
	z <- data.frame(lat = unlist(z1[,1]), lon = unlist(z1[,2]), class = z2, bp = z3)
	z$color <- ifelse(z$class == "f",
					  "limegreen",
					  ifelse(z$class == "s",
					  	   "darkorange",
					  	   NA))
	z <- z[!(z$color %>% is.na()),]
	m <- leaflet(options = list(zoomControl = F, attributionControl = F)) %>%
		addProviderTiles(providers$CartoDB.Positron) %>%
		addScaleBar(position = "bottomright") %>%
		addSimpleGraticule(interval = 10, showOriginLabel = F) %>%
		addCircleMarkers(
			lng = z$lon,
			lat = z$lat,
			radius = 5,
			weight = 1,
			color = "black",
			opacity = .3,
			fillColor = z$color,
			fillOpacity = 1
		) %>%
		# addLabelOnlyMarkers(
		# 	# lng = z$lon,
		# 	# lat = z$lat,
		# 	# label = rownames(z),
		# 	# label = as.character(z$bp),
		# 	lng = z$lon[which(rownames(z) %in% as.character(c(306, 78)))],
		# 	lat = z$lat[which(rownames(z) %in% as.character(c(306, 78)))],
		# 	label = as.character(c("B", "A")),
		# 	labelOptions = labelOptions(noHide = T,
		# 								textsize = "30px",
		# 								direction = "top",
		# 								offset = c(12, -30),
		# 								textOnly = T)
		# ) %>%
		addCircleMarkers(
			lng = z$lon[which(rownames(z) %in% as.character(c(indA, indB)))],
			lat = z$lat[which(rownames(z) %in% as.character(c(indA, indB)))],
			radius = 6,
			weight = 3,
			color = "black",
			opacity = 1,
			fillColor = c("darkorange", "limegreen"),
			fillOpacity = 1
		) %>%
		addLegend(
			position = "topright",
			bins = 2,
			labels = c("Forest", "Savanna"),
			colors = c("limegreen", "darkorange"),
			opacity = 1
		)
	m

	# euc
	pointsCoords <- rasterToPoints(rasEuc)
	m <- leaflet(options = list(zoomControl = F, attributionControl = F)) %>%
		addProviderTiles(providers$CartoDB.Positron) %>%
		addScaleBar(position = "bottomright") %>%
		addSimpleGraticule(interval = 10, showOriginLabel = F) %>%
		addPolygons(
			data = spTransform(shpArea, proj_ll),
			color = "black",
			fillOpacity = 0,
			opacity = 1,
			weight = 2
		) %>%
		addCircleMarkers(
			lng = pointsCoords[smpOther, 1],
			lat = pointsCoords[smpOther, 2],
			radius = 3,
			fillOpacity = 1,
			fillColor = "green",
			weight = 0
		) %>%
		addCircleMarkers(
			lng = pointsCoords[smpEuc, 1],
			lat = pointsCoords[smpEuc, 2],
			radius = 3,
			fillOpacity = 1,
			fillColor = "red",
			weight = 0
		) %>%
		addCircleMarkers(
			lng = pointsCoords[c(smpEuc[16], smpOther[202]), 1],
			lat = pointsCoords[c(smpEuc[16], smpOther[202]), 2],
			radius = 6,
			weight = 3,
			color = "black",
			opacity = 1,
			fillColor = c("red", "green"),
			fillOpacity = 1
		) %>%
		addLegend(
			position = "topright",
			bins = 2,
			labels = c("Eucalytpus", "Non Eucalyptus"),
			colors = c("red", "green"),
			opacity = 1
		)
	# addPolygons(
	# 	data = shpEuc,
	# 	color = "red",
	# 	fillOpacity = .8,
	# 	opacity = 0,
	# 	weight = 2
	# )
	m

# Time series ----
	# savfor
	indA <- 306 # landsat 281, modis 306
	indB <- 78 # landsat 70, modis 78

	pdf(file = paste0(p_down, "AB-ls-timeseries.pdf"), width = 7, height = 7, pointsize = 13)
		layout(matrix(1:3, ncol = 1), widths = 1, heights = c(1,.82,1), respect = FALSE)
		par(mar = c(0, 4.1, 4.1, 2.1))
		plot(
			bfastOut$ndvi[[indA]]$Yt,
			ylim = c(0, 1),
			ylab = "NDVI",
			xlab = "",
			col = "white",
			xaxt = "n"
		)
		grid()
		box()
		lines(bfastOut$ndvi[[indA]]$Yt, col = "limegreen")
		lines(bfastOut$ndvi[[indB]]$Yt, col = "darkorange")
		legend(x = 2005.3, y = .25, horiz = T, col = c("limegreen", "darkorange"),
			   legend = c("Forest", "Savanna"), bty = "n", lwd = 1, cex = 1.3)
		par(mar = c(0, 4.1, 0, 2.1))
		plot(
			bfastOut$evi[[indA]]$Yt,
			ylim = c(0, 1),
			ylab = "EVI",
			xlab = "",
			yaxt = "n",
			col = "white",
			xaxt = "n"
		)
		grid()
		box()
		lines(bfastOut$evi[[indA]]$Yt, col = "limegreen")
		lines(bfastOut$evi[[indB]]$Yt, col = "darkorange")
		axis(side = 4)
		par(mar = c(4.1, 4.1, 0, 2.1))
		plot(
			bfastOut$evi2[[indA]]$Yt,
			ylim = c(0, 1),
			ylab = "EVI2",
			col = "white"
		)
		grid()
		box()
		lines(bfastOut$evi2[[indA]]$Yt, col = "limegreen")
		lines(bfastOut$evi2[[indB]]$Yt, col = "darkorange")
	dev.off()

	# euc
	indA <- 16
	indB <- 452

	pdf(file = paste0(p_down, "AB-ls-timeseries.pdf"), width = 10, height = 7, pointsize = 13)
		par(mar = c(4.1, 4.1, 0.1, 0.1))
		plot(bfastOut[[indB]]$Yt,
			 col = "white",
			 ylim = c(0, 1),
			 ylab = "NDVI")
		grid()
		box()
		lines(bfastOut[[indB]]$Yt,
			 col = "forestgreen")
		lines(bfastOut[[indA]]$Yt,
			  col = "red")
		legend(x = 2003.3, y = .1, horiz = T, col = c("red", "forestgreen"),
			   legend = c("Eucalyptus", "Non Eucalyptus"), bty = "n", lwd = 1, cex = 1.1)
	dev.off()

# Bfast output ----
	# savfor
	indA <- 281 # landsat 281, modis 306
	indB <- 70 # landsat 70, modis 78

	hMon <- attributes(bfastOut)$hMon

	tsA <- bfastOut$evi2[[indA]] # landsat 281, modis 306
	resA <- bfast(Yt = tsA$Yt,
				  max.iter = 1,
				  h = hMon/length(tsA$Yt))
	tsB <- bfastOut$ndvi[[indB]] # landsat 70, modis 78
	resB <- bfast(Yt = tsB$Yt,
				  max.iter = 1,
				  h = hMon/length(tsB$Yt))

	pdf(file = paste0(p_down, "A-ls-evi2-bfast.pdf"), width = 9, height = 7, pointsize = 13)
	plot(resA,
		 main = "")
	dev.off()

	pdf(file = paste0(p_down, "B-landsat-ndvi-bfast.pdf"), width = 9, height = 7, pointsize = 13)
	plot(resB,
		 main = "")
	dev.off()

	# euc
	indA <- 16
	indB <- 452

	hMon <- 48

	tsA <- bfastOut[[indA]]
	resA <- bfast(Yt = tsA$Yt,
				  max.iter = 1,
				  h = hMon/length(tsA$Yt))
	tsB <- bfastOut[[indB]]
	resB <- bfast(Yt = tsB$Yt,
				  max.iter = 1,
				  h = hMon/length(tsB$Yt))

	pdf(file = paste0(p_down, "A-modis-ndvi-bfast.pdf"), width = 9, height = 7, pointsize = 13)
	plot(resA,
		 main = "")
	dev.off()

	pdf(file = paste0(p_down, "B-modis-ndvi-bfast.pdf"), width = 9, height = 7, pointsize = 13)
	plot(resB,
		 main = "")
	dev.off()

# Dissimilarity matrices images ----
par(mar = rep(0, 4))
n <- 17
z <- dDecTest[[n]] %>% as.matrix()
z <- z/max(z)
image(
	x = z %>% apply(2, rev) %>% t(),
	col = rainbow(256, end = .68) %>% rev(),
	axes = F
)















# Error map ----

f_script <- "bfast-csv"
f_region <- "dantas" # "dantas" "eucalipto"
f_sat <- "MCD13Q1" # "LS5e7e8_SR" "MCD13Q1" "modis"
load(paste0(p_rdata, f_script, "/", f_region, "-", f_sat, ".RData"))

# k <- 5 # folds
# set.seed(1)
# folds <- as.character(classif)
# folds[which(folds == levels(classif)[1])] <- sample(cut(1:length(which(folds == levels(classif)[1])), breaks = k, labels = F))
# folds[which(folds == levels(classif)[2])] <- sample(cut(1:length(which(folds == levels(classif)[2])), breaks = k, labels = F))
# folds <- as.numeric(folds)


f_sat <- "MCD13Q1"
viName <- c("ndvi","evi","evi2")
g <- NULL

for(j in 1:3){
	load(paste0(p_rdata, scriptName, "/", f_region, "-", f_sat, "-", viName[j], "-2.RData"))

	g <- c(g,unlist(lapply(gpRes, function(x){
		lapply(x, function(y){
			c(y$rawWrongTrain,y$rawWrongVal,y$decWrongTrain,y$decWrongVal)
		})
	})))
}

a <- 20 # 20 mais errados
pos <- as.numeric(names(sort(table(g))[-c(1:(length(table(g))-a))]))

table_w <- data.frame(point=0, lat=0, lon=0)
for(i in 1:a){
	table_w[i,] <- c(pos[i],bfastOut$ndvi[[pos[i]]]$coords)
}

m <- leaflet(options = list(zoomControl = F, attributionControl = F)) %>%
	addProviderTiles(providers$CartoDB.Positron) %>%
	addScaleBar(position = "bottomright") %>%
	addSimpleGraticule(interval = 10, showOriginLabel = F) %>%
	addCircleMarkers(
		lng = table_w$lon,
		lat = table_w$lat,
		radius = 5,
		weight = 1,
		color = "black",
		opacity = .3,
		fillColor = "red",
		fillOpacity = 1
	)
m






# Dantas LANDSAT
f_sat <- "LS5e7e8_SR"
viName <- c("ndvi","evi","evi2")
g <- NULL

for(j in 1:3){
	load(paste0(p_rdata, scriptName, "/", f_region, "-", f_sat, "-", viName[j], "-2.RData"))

	g <- c(g,unlist(lapply(gpRes, function(x){
		lapply(x, function(y){
			c(y$rawWrongTrain,y$rawWrongVal,y$decWrongTrain,y$decWrongVal)
		})
	})))
}

a <- 20 # 20 mais errados
pos <- as.numeric(names(sort(table(g))[-c(1:(length(table(g))-a))]))

table_w2 <- data.frame(point=0, lat=0, lon=0)
for(i in 1:a){
	table_w2[i,] <- c(pos[i],bfastOut$ndvi[[pos[i]]]$coords)
}
















# Remote sensing images ----
	# savfor
	f_script <- "bfast-csv"
	f_region <- "dantas"
	f_sat <- "MCD13Q1"
	load(paste0(p_rdata, f_script, "/", f_region, "-", f_sat, ".RData"))

	indA <- 306 # landsat 281, modis 306
	indB <- 78 # landsat 70, modis 78

	rsLsA <- raster("D:/Documents/Dropbox/Unicamp/Projeto de Mestrado/Programas/Data/raster/savfor-A/LT5_SR/LT50060622005266/LT50060622005266.evi2.tif")
	rsMdA <- raster("D:/Documents/Dropbox/Unicamp/Projeto de Mestrado/Programas/Data/raster/savfor-A/MOD13Q1/MOD13Q1_005_2005_09_30/e973dd767ce1175ad8593bc7ce99bc5c.evi2.tif")
	rsLsB <- raster("D:/Documents/Dropbox/Unicamp/Projeto de Mestrado/Programas/Data/raster/savfor-B/LT5_SR/LT52190632005142/LT52190632005142.evi2.tif")
	rsMdB <- raster("D:/Documents/Dropbox/Unicamp/Projeto de Mestrado/Programas/Data/raster/savfor-B/MOD13Q1/MOD13Q1_005_2005_04_23/fca8c690a29991637a1a2b3f96f167c7.evi2.tif")

	values(rsLsA)[values(rsLsA) == 0] <- NA
	values(rsLsB)[values(rsLsB) == 0] <- NA

	customColPal <- colorNumeric(
		palette = viridis(128),
		na.color = "#80808000",
		alpha = T,
		domain = c(-1,1)
	)

	customColPalLeg <- colorNumeric(
		palette = viridis(128),
		domain = c(1,-1)
	)

	m <- leaflet(options = list(zoomControl = F, attributionControl = F)) %>%
		addTiles(
			urlTemplate = "http://{s}.google.com/vt/lyrs=s&x={x}&y={y}&z={z}",
			attribution = paste(
				"Map data &copy;",
				substr(Sys.Date(), 1, 4),
				" Google, Imagery &copy;",
				substr(Sys.Date(), 1, 4),
				" TerraMetrics"
			),
			options = list(
				maxZoom = 20,
				subdomains = c("mt0", "mt1", "mt2", "mt3")
			)
		) %>%
		addScaleBar(position = "bottomright") %>%
		addSimpleGraticule(interval = 10, showOriginLabel = F) %>%
		addRasterImage(
			x = projectRaster(from = rsMdB, crs = proj_ll),
			colors = customColPal
		) %>%
		addLegend(
			position = "topright",
			pal = customColPalLeg,
			values = c(-1,1),
			title = "EVI2",
			opacity = 1
		)
	m

	# euc
	ras <- raster(paste0(p_raster, "eucalipto/MOD13Q1/MOD13Q1_005_2012_08_28/d69e24f856c9dbfeb279a459ad134250.ndvi.tif"))
	shpArea <- paste0(p_shape, "emprise_classification2014.shp") %>% shapefile()
	shpArea <- spTransform(shpArea, proj_ll)
	shpEuc <- paste0(p_shape, "eucalipto.shp") %>% shapefile()
	shpOther <- gDifference(shpArea, shpEuc)

	areaEuc <- sapply(shpEuc@polygons[[1]]@Polygons, function(x) {
		x@area
	})

	top <- order(areaEuc, decreasing = T)[1:150]
	shpEucTop <- shpEuc
	shpEucTop@polygons[[1]]@Polygons <- shpEucTop@polygons[[1]]@Polygons[top]

	suppressWarnings({
		rasEuc <- mask(ras, shpEucTop, updatevalue = -1, updateNA = -1)
		rasOther <- mask(ras, shpOther, updatevalue = -1, updateNA = -1)
	})

	pointsCoords <- rasterToPoints(rasEuc)
	pointsEuc <- which(values(rasEuc) != -1, arr.ind = T)
	pointsOther <- which(values(rasOther) != -1, arr.ind = T)

	set.seed(1)
	smpEuc <- sample(pointsEuc, 250)
	smpOther <- sample(pointsOther, 1000)

	rsMd <- raster("D:/Documents/Dropbox/Unicamp/Projeto de Mestrado/Programas/Data/raster/eucalipto/MOD13Q1/MOD13Q1_005_2012_08_28/d69e24f856c9dbfeb279a459ad134250.ndvi.tif")
	rsMd <- mask(projectRaster(from = rsMd, crs = proj_ll), spTransform(shapefile(paste0(p_shape, "emprise_classification2014.shp")), proj_ll))

	customColPal <- colorNumeric(
		palette = viridis(128),
		na.color = "#80808000",
		alpha = T,
		domain = c(-1,1)
	)

	customColPalLeg <- colorNumeric(
		palette = viridis(128),
		domain = c(1,-1)
	)

	m <- leaflet(options = list(zoomControl = F, attributionControl = F)) %>%
		addTiles(
			urlTemplate = "http://{s}.google.com/vt/lyrs=s&x={x}&y={y}&z={z}",
			attribution = paste(
				"Map data &copy;",
				substr(Sys.Date(), 1, 4),
				" Google, Imagery &copy;",
				substr(Sys.Date(), 1, 4),
				" TerraMetrics"
			),
			options = list(
				maxZoom = 20,
				subdomains = c("mt0", "mt1", "mt2", "mt3")
			)
		) %>%
		addScaleBar(position = "bottomright") %>%
		addSimpleGraticule(interval = 10, showOriginLabel = F) %>%
		addRasterImage(
			x = rsMd,
			colors = customColPal
		) %>%
		addCircleMarkers(
			lng = pointsCoords[c(smpEuc[16], smpOther[202]), 1],
			lat = pointsCoords[c(smpEuc[16], smpOther[202]), 2],
			radius = 4,
			weight = 2,
			color = "black",
			opacity = 1,
			fillOpacity = 0
		) %>%
		addLegend(
			position = "topright",
			pal = customColPalLeg,
			values = c(-1,1),
			title = "NDVI",
			opacity = 1
		)
	m

