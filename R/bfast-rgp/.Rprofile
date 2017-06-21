# ------------------------------------------------------------- IMPDEF ----

scriptName <- "bfast-rgp"

# load packages, install if necessary
packs <- c(
	"bfast",
	"car",
	"caret",
	"cluster",
	"colorspace",
	"FastKNN",
	"igraph",
	"leaflet",
	"pryr",
	"raster",
	"rgp",
	"robust",
	"stringr",
	"TSdist",
	"xtable"
)

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

tsDistances <- data.frame(
	name = c(
		"cor", "fourier", "sts", "spec.llr", "int.per", "acf", "euclidean",
		"manhattan", "infnorm", "cid", "cort", "mindist.sax", "pacf",
		"ccor", "per", "ar.mah.statistic", "ar.pic", "erp", "edr", "dtw"
	),
	acronym = c(
		"PEC", "FOU", "STS", "LLR", "INP", "ACF", "EUC",
		"MAN", "INF", "CID", "CRT", "SAX", "PCF", "COR",
		"PER", "ARM", "ARP", "ERP", "EDR", "DTW"
	),
	type = c(
		"eq", "eq", "eq", "eq", "eq", "un", "eq", "eq", "eq", "eq", "eq",
		"un", "un", "un", "eq", "un", "un", "un", "un", "un"
	)
)

gpFunctions <- c("+", "-", "*", "sSqrt", "sDiv", "sLn", "abs")

# ---------------------------------------------------------- FUNCTIONS ----

# GP raw vs decomposition comparison
gpPlotFold <- function(res, titleText = "") {
	nFolds <- length(res)

	defaultPar <- par()
	par(mar = c(5, 4, 4, 10) + .1)

	plot(0, 0,
		 xlim = c(1, nFolds),
		 ylim = c(0, 1),
		 col = "white",
		 xlab = "fold",
		 ylab = "balanced accuracy",
		 yaxp = c(0,1,10),
		 main = titleText)
	abline(h = seq(0,1,.1), v = 1:5, col = "gray", lty = 3)
	box()

	par(xpd = T)

	# baseline
		basAcc <- sapply(res, function(x) { x$baseline })
		for(i in 1:nrow(basAcc)) {
			lines(
				y = basAcc[i,],
				x = 1:nFolds,
				lwd = 1,
				col = rainbow(nrow(basAcc))[i]
			)
		}
	# GPdec
		decAcc <- sapply(res, function(x) { x$decAcc })
		lines(decAcc,
			  col = "black",
			  lwd = 2,
			  lty = 1)
	# GPraw
		rawAcc <- sapply(res, function(x) { x$rawAcc })
		lines(rawAcc,
			  col = "black",
			  lwd = 2,
			  lty = 5)

	legend(
		"right",
		inset = c(-0.41, 0),
		legend = c(
			expression("GP"["dec."]),
			expression("GP"["raw"]),
			as.character(tsDistances$acronym)
		),
		col = c(rep("black", 2), rainbow(9)),
		lty = c(1, 5, rep(1, 9)),
		lwd = c(2, 2, rep(1, 9)),
		bty = "n"
	)

	suppressWarnings(par(defaultPar))
}

gpPlotSynTree <- function(res, type = c("raw", "dec"), val = NULL) {
	type <- match.arg(type)

	if(type == "dec") {
		if(val) {
			z <- body(res$decBestIndVal)
		} else {
			z <- body(res$decBestIndTrain)
		}
	}
	if(type == "raw") {
		if(val) {
			z <- body(res$rawBestIndVal)
		} else {
			z <- body(res$rawBestIndTrain)
		}
	}

	g1 <- exprToGraph(z)

	x1 <- suppressWarnings(which(!is.na(as.numeric(g1$vertices))))
	x2 <- suppressWarnings(round(as.numeric(g1$vertices), 1))
	g1$vertices[x1] <- as.character(x2[!is.na(x2)])
	g2 <- graph(g1$edges, n = length(g1$vertices))
	V(g2)$label <- g1$vertices
	V(g2)$col <-
		ifelse(
			g1$vertices %in% gpFunctions,
			"white",
			ifelse(!is.na(x2), "lightgoldenrod1", "SkyBlue")
		)
	plot(g2,
		 layout = layout.reingold.tilford(g2, root = 1),
		 edge.arrow.size = 0.7,
		 asp = 0,
		 margin = 0,
		 vertex.label.cex = 0.6,
		 vertex.size = 6,
		 vertex.color = V(g2)$col,
		 vertex.label.color = "black")
}
