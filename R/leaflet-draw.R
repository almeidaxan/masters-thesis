# ------------------------------------------------------------- README ----

# Verificar se pasta
# 		leaflet/htmlwidgets/lib/ contem a pasta Leaflet.draw/
# E se a pasta
# 		leaflet/htmlwidgets/plugins/ contem a pasta leaflet-draw-plugin/
#
# Se nao, baixar manualmente as pastas do GitHub abaixo e coloca-las nos
# respectivos locais.
#
# GitHub:
# https://github.com/RCura/leaflet/tree/master/inst/htmlwidgets
#
# Apos isso, rodar o codigo abaixo.

# ------------------------------------------------------------- IMPDEF ----

# load packages, install if necessary
packs <- c("leaflet",
		   "devtools",
		   "viridis",
		   "rgeos",
		   "raster")

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

# --------------------------------------------------------------- MAIN ----

# leafletDrawDependencies <- function() {
# 	list(
# 		htmltools::htmlDependency(
# 			"Leaflet.draw",
# 			"0.2.3",
# 			system.file("htmlwidgets/lib/Leaflet.draw/dist/", package = "leaflet"),
# 			script = "leaflet.draw.js",
# 			stylesheet = "leaflet.draw.css"
# 		),
# 		htmltools::htmlDependency(
# 			"leaflet-draw-plugin",
# 			"0.0.1",
# 			system.file("htmlwidgets/plugins/leaflet-draw-plugin/", package = "leaflet"),
# 			script = "leaflet-draw-plugin.js"
# 		)
# 	)
# }
#
# addDrawToolbar <- function(map,
# 						   layerID = "drawLayer",
# 						   position = c("topleft", "topright", "bottomleft", "bottomright"),
# 						   polyline = F,
# 						   polygon = F,
# 						   rectangle = F,
# 						   circle = F,
# 						   marker = T,
# 						   edit = F,
# 						   remove = F) {
# 	position = match.arg(position)
# 	map$dependencies <- c(map$dependencies, leafletDrawDependencies())
# 	map$drawToolbar <- T
# 	invokeMethod(map, getMapData(map), method = "addDrawToolbar", layerID,
# 				 position, polyline, polygon, rectangle, circle,
# 				 marker, edit, remove)
# }

# sat <- "LT50060621986310"
#
# z <-
# 	raster(
# 		paste0(
# 			"D:/Documents/Dropbox/Unicamp/Projeto de Mestrado/Programas/Data/raster/gpA/LT5_SR/",
# 			sat,
# 			"/",
# 			sat,
# 			".ndvi.tif"
# 		)
# 	)

# z <-
# 	raster(
# 		paste0(
# 			"D:/Documents/Dropbox/Unicamp/Projeto de Mestrado/Programas/Data/raster/gpB/MOD13Q1/MOD13Q1_005_2000_03_21/9865fb97181fa030c30547b06a20172d.evi.tif"
# 		)
# 	)
#
# shp <- shapefile(paste0(p_shape, "gpB.shp"))
#
# cusPal <- colorNumeric(
# 	palette = rev(terrain.colors(128)),
# 	domain = c(0,1)
# )

shp <- shapefile("D:/Documents/Dropbox/Unicamp/Projeto de Mestrado/Programas/Data/shape/fray-jorge-park.shp")
rst <- raster("D:/Documents/Dropbox/Unicamp/Projeto de Mestrado/Programas/Data/raster/fray-jorge-park/LT5_L1T_TOA_FMASK/LT50010812011104COA00/LT50010812011104COA00.ndvi.tif")
rst <- projectRaster(from = rst, crs = proj_ll)
ele <- raster("D:/Documents/Dropbox/Unicamp/Projeto de Mestrado/Programas/Data/raster/fray-jorge-park/SRTMGL1_003/SRTMGL1_003.elevation.tif")

rst[values(rst) == 0] <- NA

slope <- terrain(ele, opt = "slope")
aspect <- terrain(ele, opt = "aspect")
hill <- hillShade(slope, aspect, 45, 120)

# y <- raster("D:/Documents/Dropbox/Unicamp/Projeto de Mestrado/Projetos/Catarina_Alexandre_2016/Dados_Loic/raster/LC80010622014272_ndmi_sub.tif")
# y <- projectRaster(from = y, crs = proj_ll)
# x <- raster("D:/Documents/Dropbox/Unicamp/Projeto de Mestrado/Projetos/Catarina_Alexandre_2016/Dados_Loic/classif_1984/landCover_tefe_1984.tif")
# x <- projectRaster(from = x, crs = proj_ll, method = "ngb")
# x[values(x) %in% c(7,8)] <- NA

# shpExt <- extent(a3)
# shpExt <- as(shpExt, 'SpatialPolygons')

customColPal <- colorNumeric(
	palette = viridis(128),
	na.color = "#80808000",
	alpha = T,
	domain = c(-1,1)
)

customColEle <- colorNumeric(
	palette = colorRampPalette(c("#000000", "#FFFFFF"))(256),
	na.color = "#80808000",
	alpha = T,
	domain = c(minValue(elev2), maxValue(elev2))
)

customColPalLeg <- colorNumeric(
	palette = viridis(128),
	domain = c(-1,1)
)

customColEleLeg <- colorNumeric(
	palette = colorRampPalette(c("#000000", "#FFFFFF"))(256),
	domain = c(minValue(elev2), maxValue(elev2))
)

m <- leaflet(options = list(zoomControl = F, attributionControl = F)) %>%
	# addProviderTiles(providers$CartoDB.Positron) %>%
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
	# addRasterImage(
	# 	x = elev2,
	# 	colors = customColEle
	# ) #%>%
	addRasterImage(
		x = mask(rst, shp),
		colors = customColPal,
		maxBytes = 10 * 1024 * 1024
	)#%>%
	# addPolylines(
	# 	data = shp,
	# 	weight = 3,
	# 	color = "red",
	# 	fillOpacity = 1,
	# 	opacity = 1
	# )
m

m %>%
addScaleBar(position = "bottomright") %>%
addSimpleGraticule(interval = 10, showOriginLabel = F) %>%
addLegend(
	position = "topright",
	pal = customColPalLeg,
	values = c(-1, 1),
	title = "NDVI",
	opacity = 1
)
# addLegend(
# 	position = "topright",
# 	colors = rainbow(6) %>% str_sub(1,7),
# 	labels = c("Agriculture", "Bare Soil", "Cloud", "Forest", "Cloud Shadow", "Urban"),
# 	title = "Class",
# 	opacity = 1
# )

a2 <- SpatialPolygons(shp@polygons, proj4string = crs(shp))
a3 <- gUnaryUnion(a2)
a3 <- spTransform(a3, proj_ll)



a4 <- a2
a4@polygons[setdiff(1:length(a4@polygons), cla1984$id[cla1984$claMode == 4])] <- NULL
a4@plotOrder <- 1:(a4@polygons %>% length())
a5 <- gUnaryUnion(a4)
a5 <- spTransform(a5, proj_ll)
