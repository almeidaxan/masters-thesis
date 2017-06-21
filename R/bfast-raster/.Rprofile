# -------------------------------------------------------------------------
# Author:
# 		Alexandre Esteves Almeida (2016)
# ------------------------------------------------------------- IMPDEF ----

scriptName <- "bfast-raster"

# load packages, install if necessary
packs <- c("bfastSpatial",
		   "lubridate",
		   "strucchange",
		   "raster",
		   "ncdf4",
		   "Cairo",
		   "animation",
		   "directlabels",
		   "ggplot2",
		   "robustbase",
		   "RColorBrewer")

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

# define a HSV pallete spectrum
hsvPalette <- NULL
for(i in seq(0, 1, length.out = 12)) {
	hsvPalette <- c(hsvPalette, grDevices::hsv(i, 1, 1))
}

# periodogram function
pgm <- function(serie) {
	n <- length(serie)
	m <- n/2 - 1
	I <- abs(fft(serie)/sqrt(n))^2  # periodogram
	P <- (4/n)*I                    # scaled-periodogram
	f <- 0:m/n
	list(f=f, P=P[1:(m+1)])
}

# function to plot results in a raster with legend
plotRes <- function(so,
					elev,
					col,
					legTitle = "",
					legType = c("legend", "colourbar"),
					legBias = 1,
					legNCol = NULL,
					saveFile = F,
					fileName = NULL) {
	legType <- match.arg(legType, legType)

	colFunc <- colorRampPalette(col, bias = legBias)

	sotmp <- projectRaster(so, crs = proj_ll, method = "ngb")
	sotmp <- rasterToPoints(sotmp) %>% data.frame()
	names(sotmp)[3] <- "fill"

	soelev <- rasterToPoints(elev) %>% data.frame()
	names(soelev)[3] <- "z"

	gg <- ggplot() +
		geom_raster(data = sotmp, mapping = aes(x = x, y = y, fill = fill)) +
		stat_contour(
			data = soelev,
			mapping = aes(x = x,
						  y = y,
						  z = z,
						  colour = ..level..),
			colour = colorRampPalette(c(rgb(.5,.5,.5), "white"))(6545),
			breaks = c(100, 200, 350, 600),
			show.legend = F,
			size = 1.1
		) +
		geom_dl(
			data = soelev,
			aes(x = x,
				y = y,
				z = z,
				label = ..level..),
			colour = colorRampPalette(c(rgb(.5,.5,.5), "white"))(6545),
			breaks = c(100, 200, 350, 600),
			method = list(
				"bottom.pieces",
				cex = .9,
				show.legend = F
			),
			stat = "contour"
		) +
		coord_equal() +
		scale_x_continuous(expand = c(0, 0)) +
		scale_y_continuous(expand = c(0, 0)) +
		theme_light() +
		theme(
			legend.position = "top"
		)
	if(legType == "legend") {
		b <- min(sotmp$fill):max(sotmp$fill)
		gg <- gg +
			scale_fill_gradientn(breaks = b,
								 colors = colFunc(128),
								 guide = legType) +
			guides(fill = guide_legend(
				label.position = "bottom",
				label.hjust = .5,
				keywidth = 1.3,
				keyheight = 1,
				title = legTitle,
				title.position = "top",
				title.hjust = .5,
				ncol = legNCol,
				nrow = 1
			))
	} else {
		gg <- gg +
			scale_fill_gradientn(colors = colFunc(128),
								 guide = legType) +
			guides(fill = guide_colorbar(
				barwidth = 20,
				barheight = .6,
				title = legTitle,
				title.position = "top",
				title.hjust = .5
			))
	}
	gg <- gg +
		labs(x = expression(paste("Longitude (", degree, ")")),
			 y = expression(paste("Latitude (", degree, ")"))) +
		theme(text = element_text(size = 18),
			  axis.text = element_text(size = 14))

	if (saveFile) {
		ggsave(
			file = paste0(p_down, fileName, ".pdf"),
			plot = gg,
			width = 9,
			height = 9
		)
	} else {
		print(gg)
	}
}
