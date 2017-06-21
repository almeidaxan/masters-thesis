# ------------------------------------------------------------- IMPORT ----

# define default CRAN mirror
local({r <- getOption("repos")
	r["CRAN"] <- "https://cran.rstudio.com/"
	options(repos = r)
})

if(!require("pacman")) {
	install.packages("pacman")
	require("pacman")
}

# always load the following packages from CRAN
packs <- c(packs, "foreach", "parallel", "beepr", "magrittr", "insertPipe")

# always load the following packages from GitHub
packsGitHub <- c(
	"dutri001/bfastSpatial",
	"vqv/ggbiplot",
	"ASBecker/insertPipe"
)

# load parallel backend package, based on OS
if(Sys.info()[["sysname"]] == "Windows") {
	# windows OS
	packs <- c(packs, "doSNOW")
} else {
	# unix-based OSs
	packs <- c(packs, "doMC")
}

# load packages using 'pacman' functions
p_load(char = packs)
p_load_gh(char = packsGitHub)

# ------------------------------------------------------------- DEFINE ----

# useful paths
if(Sys.info()[["sysname"]] == "Windows") {
	# windows OS
	path <- "D:/Documents/Dropbox/Unicamp/Projeto de Mestrado/Programas/Data"
	p_down <- "C:/Users/Alexandre/Downloads/"
} else {
	# unix-based OSs
	if(Sys.getenv("USER") == "almeida") {
		path <- "/Users/almeida/Dropbox/Unicamp/Projeto de Mestrado/Programas/Data"
		p_down <- "/Users/almeida/Downloads/"
	}
	if(Sys.getenv("USER") == "Menini") {
		path <- "/Users/almeida/Dropbox/Unicamp/Pos/Programas/Data"
		p_down <- "/Users/almeida/Downloads/"
	}
}
p_csv <- paste0(path, "/csv/")
p_raster <- paste0(path, "/raster/")
p_shape <- paste0(path, "/shape/")
p_stack <- paste0(path, "/stack/")
p_rdata <- paste0(path, "/rdata/")
p_nc <- paste0(path, "/nc/")

proj_ll <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
proj_utm <- function(shape) {
	p <- "+proj=utm +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 +units=m"
	paste0(p, " +zone=", zoneCalc(extent(shape)[1]))
}

# calculates UTM zone based on long coordinate
zoneCalc <- function(long) {
	(floor((long + 180)/6) %% 60) + 1
}

# lists the size (in MB) of all variables in the environment
sizeCalc <- function() {
	cbind(
		sort(
			sapply(
				ls(envir = .GlobalEnv),
				function(x) {
					round(object.size(get(x))/(10^6),1)
				}
			)
		)
	)
}

# check if objects have 'Date' class
is.date <- function(x) inherits(x, "Date")

# garbage cleaning
rm(path)
