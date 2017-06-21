# ---------------------------------------------------------- CSV MERGE ----

# load Rdata with Juan's GPVI functions
load(paste0(p_rdata, "gpvi.RData"))

# choose which csvs to merge
f_region <- "dantas"
f_vi <- c("MOD13Q1", "MYD13Q1") #c("MOD09A1", "MYD09A1")
#f_fire <- "MCD45A1"
out_vi <- paste0("MCD", substr(f_vi[1], 4, 7))
setwd(paste0(p_csv, f_region))

# create merged data folder
dir.create(paste0("final-", out_vi, "/"), showWarnings = F)
#dir.create(paste0("final-", f_fire, "/"), showWarnings = F)

# number of csvs to process
nFiles <- length(dir(f_vi[1]))

# loop through all csvs
print("Merging CSVs and preprocessing...")
pb <- txtProgressBar(min = 0,
					 max = nFiles,
					 style = 3)
dfList <- list()
for (i in 1:nFiles) {
	# read all csvs with the same name
	file_vi1 <- read.csv(dir(f_vi[1], full.names = T)[i])
	file_vi2 <- read.csv(dir(f_vi[2], full.names = T)[i])
	#file_fire <- read.csv(dir(f_fire[1], full.names = T)[i])

	# format date column
	file_vi1$date <- as.Date(file_vi1$date, format = "%Y_%m_%d")
	file_vi2$date <- as.Date(file_vi2$date, format = "%Y_%m_%d")

	# if the product ends with "09A1", calculate all VIs
	if(grep(glob2rx("*09A1"), f_vi[1]) %>% length() > 0) {
		file_vi1$ndvi <- (file_vi1$b2 - file_vi1$b1)/(file_vi1$b2 + file_vi1$b1)
		file_vi2$ndvi <- (file_vi2$b2 - file_vi2$b1)/(file_vi2$b2 + file_vi2$b1)
		file_vi1$evi <- 2.5*(file_vi1$b2 - file_vi1$b1)/(file_vi1$b2 + 6*file_vi1$b1 - 7.5*file_vi1$b3 + 1)
		file_vi2$evi <- 2.5*(file_vi2$b2 - file_vi2$b1)/(file_vi2$b2 + 6*file_vi2$b1 - 7.5*file_vi2$b3 + 1)
		file_vi1$evi2 <- 2.5*(file_vi1$b2 - file_vi1$b1)/(file_vi1$b2 + 2.4*file_vi1$b1 + 1)
		file_vi2$evi2 <- 2.5*(file_vi2$b2 - file_vi2$b1)/(file_vi2$b2 + 2.4*file_vi2$b1 + 1)
		file_vi1$gpvi <- apply(file_vi1[,2:8], MARGIN = 1, FUN = gpviM)
		file_vi2$gpvi <- apply(file_vi2[,2:8], MARGIN = 1, FUN = gpviM)
	}

	# if(nrow(file_fire) > 0) {
	# 	file_fire$date <- year(as.Date(file_fire$date, format = "%Y_%m_%d"))
	# 	file_fire <- unique(file_fire)
	# 	for(j in as.numeric(names(which(table(file_fire$date) > 1)))) {
	# 		if(sum(file_fire$qa[file_fire$date == j] == min(file_fire$qa[file_fire$date == j])) == 1) {
	# 			ind <- which.min(file_fire$qa[file_fire$date == j])
	# 			bd <- file_fire$burndate[file_fire$date == j][ind]
	# 		} else {
	# 			bd <- floor(mean(file_fire$burndate[file_fire$date == j]))
	# 		}
	# 		file_fire <- file_fire[!(file_fire$date == j), ]
	# 		file_fire[nrow(file_fire)+1, ] <- c(j, bd, 0)
	# 	}
	# 	rownames(file_fire) <- NULL
	# 	file_fire$qa <- NULL
	# 	file_fire <- file_fire[order(file_fire$date), ]
	# 	colnames(file_fire)[2] <- "burnday"
	# 	file_fire$burndate <- sort(as.Date(paste(file_fire$date, file_fire$burnday), format="%Y %j"))
	#
	# 	# save fire data
	# 	write.csv(file_fire,
	# 			  file = paste0("final-", f_fire, "/", dir(f_vi[1])[i]),
	# 			  row.names = F,
	# 			  quote = F)
	# }

	# merge two vi csvs
	final_vi <- rbind(file_vi1, file_vi2)
	final_vi <- final_vi[order(final_vi$date),]

	# reduce date to months
	timeYM <- strftime(final_vi$date, "%Y-%m")
	final_vi <- data.frame(time = timeYM, final_vi)
	final_vi$date <- NULL
	final_vi <- aggregate(x = final_vi[,2:ncol(final_vi)],
						  by = list(final_vi[,1]),
						  FUN = median)
	names(final_vi)[1] <- "time"

	# create an empty continuous year-month df
	mYear <- as.numeric(substr(head(final_vi$time, 1), 1, 4))
	MYear <- as.numeric(substr(tail(final_vi$time, 1), 1, 4))
	mMonth <- as.numeric(substr(head(final_vi$time, 1), 6, 7))
	MMonth <- as.numeric(substr(tail(final_vi$time, 1), 6, 7))
	fullYear <- rep(mYear:MYear, each = 12)
	fullMonth <- rep(1:12, times = length(mYear:MYear))
	# correcting the first year
	if(mMonth>1) {
		fullYear <- fullYear[-(1:(mMonth-1))]
		fullMonth <- fullMonth[-(1:(mMonth-1))]
	}
	# correcting the last year
	fullYear <- fullYear[-((length(fullYear)-12+MMonth+1):length(fullYear))]
	fullMonth <- fullMonth[-((length(fullMonth)-12+MMonth+1):length(fullMonth))]
	dfull <- data.frame(time = paste0(fullYear,"-",sprintf("%02d", fullMonth)))

	# merging final_vi df with the full year-month df
	d <- merge(dfull, final_vi, all = T)
	d$time <- NULL

	dts <- na.approx(
		ts(d,
		   start = c(mYear, mMonth),
		   end = c(MYear, MMonth),
		   frequency = 12),
		rule = 2
	)

	# preparing df to save as csv
	d <- data.frame(dfull, round(dts,4))

	# save merged vi data
	# write.csv(d,
	# 		  file = paste0(paste0("final-", out_vi, "/"), dir(f_vi[1])[i]),
	# 		  row.names = F)

	dfList[[i]] <- dts

	# update progress bar
	setTxtProgressBar(pb, i)
}
close(pb)

# subsetting time series period
#dfList <- lzpply(dfList, function(x) { window(x, start = c(2005, 1)) } )

# -------------------------------------------------------------- BFAST ----

# loads csv with:
#	- coordiantes, in lat, long order
#	- classes: f=forest, s=savanna, t=transition, m=mixed
d <- read.csv(paste0(p_csv, f_region, ".csv"))

# minimun segment size (months) used to calculate h parameter in bfast
hMon <- 18

# precipitation data
cruP <- nc_open(paste0(p_nc, "prec_1961_2014_SA.nc"), readunlim = F)
prec <- ncvar_get(cruP, "pre")
lon <- ncvar_get(cruP, "lon")
lat <- ncvar_get(cruP, "lat")
tim <- ncvar_get(cruP, "time")
tim <- strftime(as.Date(tim, origin = "1900-01-01"), "%Y-%m")

# VIs columns indexes
viCol <- grep(glob2rx("*vi*"), colnames(dfList[[1]]))

# run bfast for each vi time series
bfastOut <- vector("list", length(viCol))
pb <- txtProgressBar(min = 0,
					 max = length(dfList)*length(viCol),
					 style = 3)
progressBar <- function(n) setTxtProgressBar(pb, n)
for(j in viCol) {
	l <- j - viCol[1] + 1
	bfastOut[[l]] <- list()
	for(i in 1:length(dfList)) {
		bfastOut[[l]][[i]] <- list()

		res <- bfast(Yt = dfList[[i]][,j],
					 max.iter = 1,
					 h = hMon/nrow(dfList[[i]]))

		# 1) original time series (Yt)
			bfastOut[[l]][[i]]$Yt <- res$Yt
		# 2) trend component (Tt)
			bfastOut[[l]][[i]]$Tt <- res$output[[1]]$Tt
		# 3) seasonal component (St)
			bfastOut[[l]][[i]]$St <- res$output[[1]]$St
		# 4) trend breakpoints (bp.T)
			if(!res$nobp$Vt) {
				bfastOut[[l]][[i]]$bp.T <- res$output[[1]]$bp.Vt$breakpoints
			} else {
				bfastOut[[l]][[i]]$bp.T <- NA
			}
		# 5) seasonal breakpoints (bp.S)
			if(!res$nobp$Wt) {
				bfastOut[[l]][[i]]$bp.S <- res$output[[1]]$bp.Wt$breakpoints
			} else {
				bfastOut[[l]][[i]]$bp.S <- NA
			}
		# 6) trend bp magnitudes matrix, columns described below:
			# [,1] -> value before bp
			# [,2] -> value after bp
			# [,3] -> magnitude of the bp
			bfastOut[[l]][[i]]$mag.T <- res$Mags
		# 7) centroid of latlong coordinates
			bfastOut[[l]][[i]]$coords <- data.frame(lat = d$lat[i], lon = d$long[i], row.names = "")
		# 8) precipitation time series
			lon_ext <- which.min(abs(lon - bfastOut[[l]][[i]]$coords$lon))
			lat_ext <- which.min(abs(lat - bfastOut[[l]][[i]]$coords$lat))

			precTs <- prec[lon_ext, lat_ext, 1:dim(prec)[3]]

			mYear <- as.numeric(substr(head(tim,1), 1, 4))
			MYear <- as.numeric(substr(tail(tim,1), 1, 4))
			mMonth <- as.numeric(substr(head(tim,1), 6, 7))
			MMonth <- as.numeric(substr(tail(tim,1), 6, 7))

			precTs <- ts(precTs,
						 start = c(mYear, mMonth),
						 end = c(MYear, MMonth),
						 frequency = 12)

			# supress 'window' function warnings
			options(warn = -1)
			precTs <- window(precTs,
							 tsp(bfastOut[[l]][[i]]$Yt)[1],
							 tsp(bfastOut[[l]][[i]]$Yt)[2])
			options(warn = 0)

			bfastOut[[l]][[i]]$prec <- precTs
		# 9) class (f=forest, s=savanna, t=transition)
			bfastOut[[l]][[i]]$class <- d$class[i]
		# 10) observation/comments
			bfastOut[[l]][[i]]$obs <- d$obs[i]

		# update progress bar
		setTxtProgressBar(pb, i+(l-1)*length(dfList))
	}
}
beep()
close(pb)

# defining bfastOut attributes
attributes(bfastOut) <- list(hMon = hMon)
names(bfastOut) <- colnames(dfList[[1]])[viCol]

# save rdata only with bfastOut list
dir.create(paste0(p_rdata, scriptName), showWarnings = F)
save(file = paste0(p_rdata, scriptName, "/", f_region, "-MCD", substr(f_vi[1], 4, 7), ".RData"),
	  list = "bfastOut",
	  envir = .GlobalEnv)
