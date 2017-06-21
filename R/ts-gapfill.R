# -------------------------------------------------------------------------(#HEADER)
# Script:
# 		ts-gapfill
# Author:
# 		Alexandre Esteves Almeida (2016)
# Description:
# 		- Script to gap-fill univariate time series
#		- Based on Pascual-Granado et al. (2014) MIARMA paper
#		- Application in e-Phenology project data

# -------------------------------------------------------------------------(#IMPDEF)

# load packages, install if necessary
packs <- c("forecast",
			  "plyr",
			  "beepr")

# source import-define script to load required libraries and define paths
if(Sys.info()[["sysname"]] == "Windows") {
	# windows OS
	source("D:/Documents/Dropbox/Unicamp/Projeto de Mestrado/Programas/R/import-define.R")
} else {
	# unix-based OSs
	if(Sys.getenv("RSTUDIO_USER_IDENTITY") == "almeida") {
		source("/Users/almeida/Dropbox/Unicamp/Projeto de Mestrado/Programas/R/import-define.R")
	}
	if(Sys.getenv("RSTUDIO_USER_IDENTITY") == "Menini") {
		source("/Users/almeida/Dropbox/Unicamp/Pos/Programas/R/import-define.R")
	}
}

# -------------------------------------------------------------------------(#FUNCS)

# gap-filling function
gap_filling <- function(serie,
								peso='v',
								plot=T,
								plot_ind=F,
								lag_val=5,
								beep_ind=T) {
	# find NA data
	pos_na <- diff(append(which(is.na(serie)),0))
	pos_na[pos_na==1] <- 0
	pos_na[length(which(is.na(serie)))] <- 1
	pos_na[pos_na>0] <- 1
	pos_na <- pos_na*which(is.na(serie))

	len_na <- diff(append(which(pos_na>0), 0, after=0))
	pos_na <- which(is.na(serie))

	par_na <- data.frame(l=rep(NA,length(len_na)),r=rep(NA,length(len_na)))

	aux_1 <- 1
	for(i in 1:length(len_na)) {
		aux_2 <- aux_1 + len_na[i] - 1
		par_na[i,1] <- pos_na[aux_1]
		par_na[i,2] <- pos_na[aux_2]
		aux_1 <- aux_2 + 1
	}

	serie_gf <<- serie
	set.seed(1)
	for (i in 1:nrow(par_na)) {
		if (i > 1) {
			aux_f1 <- par_na[i - 1, 2] + 1
		} else {
			aux_f1 <- 1
		}
		aux_f2 <- par_na[i, 1] - 1
		aux_b1 <- par_na[i, 2] + 1
		if (i < nrow(par_na)) {
			aux_b2 <- par_na[i + 1, 1] - 1
		} else {
			aux_b2 <- length(serie)
		}

		# f = forward, b = backward
		fit_f <<- auto.arima(serie_gf[aux_f1:aux_f2], ic="bic", allowdrift = F)
		# b indexes are reversed, because it's a "backwards" projection
		fit_b <<- auto.arima(serie_gf[aux_b2:aux_b1], ic="bic", allowdrift = F)

		tam_na <- par_na[i,2] - par_na[i,1] + 1

		pred_f <- predict(fit_f, tam_na)$pred
		pred_b <- predict(fit_b, tam_na)$pred
		se_f <-  predict(fit_f, tam_na)$se
		se_b <-  predict(fit_b, tam_na)$se

		# adjusting b time to equal f
		tsp(pred_f) <- c(par_na[i,1], par_na[i,2], 1)
		tsp(pred_b) <- tsp(pred_f)

		est_f <- pred_f + sapply(c(se_f), function(x) rnorm(1,0,x))
		est_b <- pred_b + sapply(c(se_b), function(x) rnorm(1,0,x))

		# defining a weight function [START]
		# ----------------------------------
		## information criterion (IC) weight
		p_f1 <- abs(fit_f$bic)/(abs(fit_f$bic) + abs(fit_b$bic))
		p_b1 <- abs(fit_b$bic)/(abs(fit_f$bic) + abs(fit_b$bic))
		if (is.na(p_f1)) {
			p_b1 <- 1
			p_f1 <- 0
		}
		if (is.na(p_b1)) {
			p_b1 <- 0
			p_f1 <- 1
		}

		if (tam_na > 1) {
			## distance weight (only makes sense for more than an obs)
			if(peso == 'd') {
				p_f2 <- rep(0, tam_na)
				half <- ceiling(tam_na/2)
				p_f2[1:half] <- seq(1, p_f1, length.out = half)
				p_f2[(half + 1):tam_na] <- seq(p_f1, 0, length.out = tam_na - half)
				p_b2 <- 1 - p_f2

				serie_gf[par_na[i,1]:par_na[i,2]] <<- as.numeric(est_f*p_f2 + est_b*p_b2)
			}
			## variability weight (only makes sense for more than an obs)
			if(peso == 'v') {
				quant_f <- serie_gf[aux_f1:aux_f2]
				quant_f2 <- quant_f[quant_f >= quantile(quant_f,.05) & quant_f <= quantile(quant_f,.95)]
				quant_b <- serie_gf[aux_b1:aux_b2]
				quant_b2 <- quant_b[quant_b >= quantile(quant_b,.05) & quant_b <= quantile(quant_b,.95)]

				if(length(quant_f2)>1) quant_f <- quant_f2
				if(length(quant_b2)>1) quant_b <- quant_b2

				# verifies if f&b, f, and b are just one obs
				if(length(quant_f)==1 & length(quant_b)==1) {
					p_f3 <- p_b3 <- 0.5
				} else if (length(quant_f)==1) {
					p_f3 <- 0
					p_b3 <- 1
				} else if (length(quant_b)==1) {
					p_f3 <- 1
					p_b3 <- 0
				} else {
					sd_mean_fb <- mean(c(sd(quant_f),sd(quant_b)))
					sd_est_f <- sd(est_f)
					sd_est_b <- sd(est_b)
					p_f3 <- 1-abs(sd_est_f-sd_mean_fb)/(abs(sd_est_f-sd_mean_fb) + abs(sd_est_b-sd_mean_fb))
					p_b3 <- 1-abs(sd_est_b-sd_mean_fb)/(abs(sd_est_f-sd_mean_fb) + abs(sd_est_b-sd_mean_fb))
				}

				if((length(quant_f)==1 | length(quant_b)==1) & len_na[i]>1 & i>1 & i<nrow(par_na)) {
					# go to the next for iteration and leaves the gap open
					next
				} else {
					trend <- seq(median(quant_f), median(quant_b), length.out=(par_na[i,2]-par_na[i,1]+1))

					if(length(names(fit_f$coef[1]))>0) {
						if(names(fit_f$coef[1]) == 'intercept') {
							untrend_f <- est_f-fit_f$coef[1]
						} else untrend_f <- est_f-mean(est_f)
					} else untrend_f <- est_f-mean(est_f)
					if(length(names(fit_b$coef[1]))>0) {
						if(names(fit_b$coef[1]) == 'intercept') {
							untrend_b <- est_b-fit_b$coef[1]
						} else untrend_b <- est_b-mean(est_b)
					} else untrend_b <- est_b-mean(est_b)

					serie_gf[par_na[i,1]:par_na[i,2]] <<- as.numeric(untrend_f*p_f3 + untrend_b*p_b3 + trend)
				}
			}
		} else {
			# if the gap has only one obs, only applies the IC weight
			serie_gf[par_na[i,1]:par_na[i,2]] <<- as.numeric(est_f*p_f1 + est_b*p_b1)
		}
		# --------------------------------
		# defining a weight function [END]

		par(mfrow = c(1,1))

		if(plot_ind) {
			plot_y <- serie_gf[(par_na[i,1] - lag_val):(par_na[i,2] + lag_val)]
			plot_x <- (par_na[i,1] - lag_val):(par_na[i,2] + lag_val)

			plot(y = plot_y, x = plot_x, type = "l")
			points(y = plot_y, x = plot_x, pch = 20, cex = 0.6)
			lines(est_f, lty = 2,
					col = rgb(1, 0, 0, .4))
			lines(est_b, lty=2,
					col = rgb(0, 0, 1, .4))
			if (tam_na > 1) {
				if(peso == 'd') {
					lines(est_f*p_f2 + est_b*p_b2, col = rgb(.8,0,.8))
					points(est_f*p_f2 + est_b*p_b2, col = rgb(.8,0,.8), pch = 20, cex = 0.6)
				}
				if(peso == 'v') {
					lines(untrend_f*p_f3 + untrend_b*p_b3 + trend, lwd = 2, col = rgb(.8,0,.8))
					points(untrend_f*p_f3 + untrend_b*p_b3 + trend, col = rgb(.8,0,.8), pch = 20, cex = 0.6)
				}
			} else {
				points(est_f*p_f1 + est_b*p_b1,
						 col = rgb(.8,0,.8),
						 pch = 20, cex = 0.6)
			}
		}
	}

	if(plot) {
		par(xpd=T)
		plot(serie_gf, type = "l", xlab = "Time (year)", ylab = "measurement", col=rgb(.4,.4,1))
		lines(serie)
		legend("top", inset=c(0,-.2), lwd=c(2,2), legend=c("Original","Gap-filling"), horiz=T, col=c("black",col=rgb(.4,.4,1)), cex=.8)
		par(xpd=F)
	}

	if(beep_ind) beep()
}

# function to verify if a year is a leap-year
is.leapyear <- function(year) {
	return(((year %% 4 == 0) & (year %% 100 != 0)) | (year %% 400 == 0))
}

# -------------------------------------------------------------------------(#MAIN)

file <- "d.csv"
p_file <- paste0(p_csv, file)
p_file_out <- paste0(p_csv, "gapfilled_", file)

# input phenology data
d <- read.csv(p_file, header=T, sep=",")

# exclude days without photos
d <- d[(d[,3] != 0) & !is.na(d[,3]),]

# create days with no obs that were excluded in the previous step
days_count <- NULL
full_day <- NULL
full_year <- NULL
for(i in min(d$year):max(d$year)) {
	j <- i - min(d$year) + 1
	# quantity of days in the full year
	days_count[j] <- 365 + is.leapyear(i)
	# fix the quantity of days of years from the extremes
	if(i == min(d$year)) {
		full_day <- c(full_day, min(d$day[d$year==i]):days_count[j])
		days_count[j] <- days_count[j] - min(d$day[d$year==i]) + 1
	} else if(i == max(d$year)) {
		full_day <- c(full_day, 1:max(d$day[d$year==i]))
		days_count[j] <- max(d$day[d$year==i])
	} else {
		full_day <- c(full_day, 1:days_count[j])
	}
}

full_year <- rep(min(d$year):max(d$year), days_count)
full <- data.frame(year = full_year, day = full_day)
d_merge <- merge(full, d, all=T)
# ignore additional day (366) from leap-years
# ts functions doesn't handle well these days
d_merge <- d_merge[-which(d_merge$day==366),]

# creating the ts
ts_meanR <- ts(d_merge$meanR, start=c(d_merge$year[1],d_merge$day[1]), end=c(tail(d_merge$year,1),tail(d_merge$day,1)), frequency=365)
ts_meanG <- ts(d_merge$meanG, start=c(d_merge$year[1],d_merge$day[1]), end=c(tail(d_merge$year,1),tail(d_merge$day,1)), frequency=365)
ts_meanB <- ts(d_merge$meanG, start=c(d_merge$year[1],d_merge$day[1]), end=c(tail(d_merge$year,1),tail(d_merge$day,1)), frequency=365)
ts_relR <- ts(d_merge$relR, start=c(d_merge$year[1],d_merge$day[1]), end=c(tail(d_merge$year,1),tail(d_merge$day,1)), frequency=365)
ts_relG <- ts(d_merge$relG, start=c(d_merge$year[1],d_merge$day[1]), end=c(tail(d_merge$year,1),tail(d_merge$day,1)), frequency=365)
ts_relB <- ts(d_merge$relG, start=c(d_merge$year[1],d_merge$day[1]), end=c(tail(d_merge$year,1),tail(d_merge$day,1)), frequency=365)
ts_excG <- ts(d_merge$excG, start=c(d_merge$year[1],d_merge$day[1]), end=c(tail(d_merge$year,1),tail(d_merge$day,1)), frequency=365)

# applying gap-filling function to ts
gap_filling(ts_meanR, peso='v')
while(is.na(prod(serie_gf))) gap_filling(serie_gf, peso='v')
final_meanR <- serie_gf
title("meanR")

gap_filling(ts_meanG, peso='v')
while(is.na(prod(serie_gf))) gap_filling(serie_gf, peso='v')
final_meanG <- serie_gf
title("meanG")

gap_filling(ts_meanB, peso='v')
while(is.na(prod(serie_gf))) gap_filling(serie_gf, peso='v')
final_meanB <- serie_gf
title("meanB")

gap_filling(ts_relR, peso='v')
while(is.na(prod(serie_gf))) gap_filling(serie_gf, peso='v')
final_relR <- serie_gf
title("relR")

gap_filling(ts_relG, peso='v')
while(is.na(prod(serie_gf))) gap_filling(serie_gf, peso='v')
final_relG <- serie_gf
title("relG")

gap_filling(ts_relB, peso='v')
while(is.na(prod(serie_gf))) gap_filling(serie_gf, peso='v')
final_relB <- serie_gf
title("relB")

gap_filling(ts_excG, peso='v')
while(is.na(prod(serie_gf))) gap_filling(serie_gf, peso='v')
final_excG <- serie_gf
title("excG")

# -------------------------------------------------------------------------(#OUTPUT)

# exporting each ts plot individually
serie = ts_relR
serie_gf = final_relR

png("relR_gf.png", type="cairo", width=780, height=540, pointsize = 20)
	par(xpd=T)
	plot(serie_gf, type = "l", xlab = "Time (year)", ylab="", col=rgb(.4,.4,1))
	lines(serie)
	legend("top", inset=c(0,-.12), lwd=c(2,2), legend=c("Original","Gap-filling"), horiz=T, col=c("black",col=rgb(.4,.4,1)), cex=.8)
	par(xpd=F)
dev.off()

# export a gap-filled csv with all the ts
output <- data.frame(ano = d_merge[,1],
                     dia = d_merge[,2],
                     relR = final_relR,
                     relG = final_relG,
                     relB = final_relB,
                     excG = final_excG )
write.csv2(output, p_file_out, row.names=F, quote=F)