library(INLA, verbose=FALSE)
library(reshape2)
library(RColorBrewer)

xWide = read.table(paste0("https://www.stat.gouv.qc.ca/statistiques/", 
                          "population-demographie/deces-mortalite/", 
                          "WeeklyDeaths_QC_2010-2020_AgeGr.csv"), 
                   sep = ";", skip = 7, col.names = c("year", "junk",  
                                                      "age", paste0("w", 1:53)))
xWide = xWide[grep("^[[:digit:]]+$", xWide$year), ] 


x = reshape2::melt(xWide, id.vars = c("year", "age"),
                   measure.vars = grep("^w[[:digit:]]+$", colnames(xWide))) 
x$dead = as.numeric(gsub("[[:space:]]", "", x$value)) 
x$week = as.numeric(gsub("w", "", x$variable))
x$year = as.numeric(x$year)
x = x[order(x$year, x$week, x$age), ] 
newYearsDay = as.Date(ISOdate(x$year, 1, 1)) 
x$time = newYearsDay + 7 * (x$week - 1)
x = x[!is.na(x$dead), ]
x = x[x$week < 53, ]
plot(x[x$age == "Total", c("time", "dead")], type = "o", log = "y")
xWide2 = reshape2::dcast(x, week + age ~ year, value.var = "dead") 
Syear = grep("[[:digit:]]", colnames(xWide2), value = TRUE)


Scol = RColorBrewer::brewer.pal(length(Syear), "Spectral") 
matplot(xWide2[xWide2$age == "Total", Syear], type = "l",lty = 1, col = Scol)
legend("topright", col = Scol, legend = Syear, bty = "n", lty = 1, lwd = 3)
dateCutoff = as.Date("2020/3/1")
xPreCovid = x[x$time < dateCutoff, ]
xPostCovid = x[x$time >= dateCutoff, ]
toForecast = expand.grid(age = unique(x$age), time = unique(xPostCovid$time),
                         dead = NA)
xForInla = rbind(xPreCovid[, colnames(toForecast)],
                 toForecast)
xForInla = xForInla[order(xForInla$time, xForInla$age),]
xForInla$timeNumeric = as.numeric(xForInla$time)
xForInla$timeForInla = (xForInla$timeNumeric - as.numeric(as.Date("2015/1/1")))/365.25 
xForInla$timeIid = xForInla$timeNumeric
xForInla$sin12 = sin(2 * pi * xForInla$timeNumeric/365.25)
xForInla$sin6 = sin(2 * pi * xForInla$timeNumeric *
                      2/365.25)
xForInla$cos12 = cos(2 * pi * xForInla$timeNumeric/365.25)
xForInla$cos6 = cos(2 * pi * xForInla$timeNumeric *2/365.25)

xForInlaTotal= xForInla[xForInla$age == 'Total', ] 


res = inla(dead ~ sin12 + sin6 + cos12 + cos6 +
             f(timeIid, prior='pc.prec', param= c(log(1.2), 0.5)) + f(timeForInla, 
                                                                      model = 'rw2', 
                                                                      prior='pc.prec', 
                                                                      param= c(0.01, 0.5)),
           data=xForInlaTotal,
           control.predictor = list(compute=TRUE, link=1), 
           control.compute = list(config=TRUE),
           # control.inla = list(fast=FALSE, strategy='laplace'),
           family='poisson')

qCols = paste0(c(0.5, 0.025, 0.975), "quant") 
rbind(res$summary.fixed[, qCols], Pmisc::priorPostSd(res)$summary[,qCols])

matplot(xForInlaTotal$time, res$summary.fitted.values[,qCols], type= "l", 
        ylim = c(1000, 1800), lty = c(1,2,2),col="black",log="y")


points(x[x$age == "Total", c("time", "dead")], col = "red", cex = 0.4)


matplot(xForInlaTotal$time, res$summary.random$timeForInla[, 
                                                           c("0.5quant", "0.975quant", 
                                                             "0.025quant")], type = "l", 
        lty = c(1, 2, 2), col = "black", ylim = c(-1, 1) *0.1)

sampleList = INLA::inla.posterior.sample(30, res, selection = list(Predictor = 0)) 
sampleIntensity = exp(do.call(cbind, Biobase::subListExtract(sampleList,"latent")))
sampleDeaths = matrix(rpois(length(sampleIntensity),
                            sampleIntensity), nrow(sampleIntensity), 
                      ncol(sampleIntensity))

matplot(xForInlaTotal$time, sampleDeaths, col = "#00000010", lwd = 2, lty = 1, 
        type = "l", log = "y")
points(x[x$age == "Total", c("time", "dead")], col = "red", cex = 0.5)
matplot(xForInlaTotal$time, sampleDeaths, col = "#00000010",
        lwd = 2, lty = 1, type = "l", log = "y", xlim = as.Date(c("2019/6/1",
                                                                  "2020/11/1")), 
        ylim = c(1, 2.3) * 1000)
points(x[x$age == "Total", c("time", "dead")], col = "red", cex = 0.5)
title(main = "Sample Deaths in Total Population")

xPostCovidTotal = xPostCovid[xPostCovid$age == "Total", ]
xPostCovidForecast = sampleDeaths[match(xPostCovidTotal$time, xForInlaTotal$time), ]
excessDeaths = xPostCovidTotal$dead - xPostCovidForecast

matplot(xPostCovidTotal$time, xPostCovidForecast, type = "l", ylim = c(1000, 2200), 
        col = "black")
points(xPostCovidTotal[, c("time", "dead")], col = "red") 
matplot(xPostCovidTotal$time, excessDeaths, type = "l",
        lty = 1, col = "#00000030")
excessDeathsSub = excessDeaths[xPostCovidTotal$time > as.Date("2020/03/01") & 
                                 xPostCovidTotal$time < as.Date("2020/06/01"), ]
excessDeathsInPeriod = apply(excessDeathsSub, 2, sum)
round(quantile(excessDeathsInPeriod))
round(quantile(excessDeaths[nrow(excessDeaths), ]))

head(x)

#FIT MODEL FOR DEATHS OF AGE >70
xForInla70= xForInla[xForInla$age == '70 years old and over', ] 

res70 = inla(dead ~ sin12 + sin6 + cos12 + cos6 +
             f(timeIid, prior='pc.prec', param= c(log(1.2), 0.5)) 
           + f(timeForInla, model = 'rw2', prior='pc.prec', param= c(0.01, 0.5)),
           data=xForInla70,
           control.predictor = list(compute=TRUE, link=1), 
           control.compute = list(config=TRUE),
           # control.inla = list(fast=FALSE, strategy='laplace'),
           family='gamma')

qCols70 = paste0(c(0.5, 0.025, 0.975), "quant") 
rbind(res$summary.fixed[, qCols], Pmisc::priorPostSd(res)$summary[,qCols])

matplot(xForInla70$time, res70$summary.fitted.values[,qCols70], type= "l", 
        ylim = c(500, 1800), lty = c(1,2,2),col="black",log="y")


points(x[x$age == "70 years old and over", c("time", "dead")], col = "red", cex = 0.4)


matplot(xForInla70$time, res70$summary.random$timeForInla[, 
                                                           c("0.5quant", "0.975quant", 
                                                             "0.025quant")], type = "l", 
        lty = c(1, 2, 2), col = "black", ylim = c(-1, 1) *0.1)

sampleList70 = INLA::inla.posterior.sample(30, res70, selection = list(Predictor = 0)) 
sampleIntensity70 = exp(do.call(cbind, Biobase::subListExtract(sampleList70,
                                                             "latent")))
sampleDeaths70 = matrix(rpois(length(sampleIntensity70),
                            sampleIntensity70), nrow(sampleIntensity70), 
                      ncol(sampleIntensity70))

matplot(xForInla70$time, sampleDeaths70, col = "#00000010", lwd = 2, lty = 1, 
        type = "l", log = "y")
points(x[x$age == "70 years old and over", c("time", "dead")], col = "red", cex = 0.5)
title(main = "Deaths in Population of Age > 70")
matplot(xForInla70$time, sampleDeaths70, col = "#00000010",
        lwd = 2, lty = 1, type = "l", log = "y", xlim = as.Date(c("2019/6/1",
                                                                  "2020/11/1")), 
        ylim = c(1, 2.3) * 700)
points(x[x$age == "70 years old and over", c("time", "dead")], col = "red", cex = 0.5)
title(main = "Sample Deaths in Population of Age > 70")

xPostCovid70 = xPostCovid[xPostCovid$age == "70 years old and over", ]
xPostCovidForecast70 = sampleDeaths70[match(xPostCovid70$time, xForInla70$time), ]
excessDeaths70 = xPostCovid70$dead - xPostCovidForecast70

matplot(xPostCovid70$time, xPostCovidForecast70, type = "l", ylim = c(700, 2200), 
        col = "black")
points(xPostCovid70[, c("time", "dead")], col = "red") 

title(main = "Projected Deaths in Population of Age > 70")

matplot(xPostCovid70$time, excessDeaths70, type = "l",
        lty = 1, col = "#00000030")

excessDeathsSub70 = excessDeaths70[xPostCovid70$time > as.Date("2020/03/01") & 
                                 xPostCovid70$time < as.Date("2020/06/01"), ]
excessDeathsInPeriod70 = apply(excessDeathsSub70, 2, sum)
round(quantile(excessDeathsInPeriod70))
round(quantile(excessDeaths70[nrow(excessDeaths70), ]))


#FIT MODEL FOR DEATHS OF AGE <50
xForInla50= xForInla[xForInla$age == '0-49 years old', ] 

res50 = inla(dead ~ sin12 + sin6 + cos12 + cos6 +
               f(timeIid, prior='pc.prec', param= c(log(1.2), 0.5)) 
             + f(timeForInla, model = 'rw2', prior='pc.prec', param= c(0.01, 0.5)),
             data=xForInla50,
             control.predictor = list(compute=TRUE, link=1), 
             control.compute = list(config=TRUE),
             # control.inla = list(fast=FALSE, strategy='laplace'),
             family='gamma')

qCols50 = paste0(c(0.5, 0.025, 0.975), "quant") 
rbind(res$summary.fixed[, qCols], Pmisc::priorPostSd(res)$summary[,qCols])

matplot(xForInla50$time, res50$summary.fitted.values[,qCols50], type= "l", 
        ylim = c(50, 110), lty = c(1,2,2),col="black",log="y")


points(x[x$age == "0-49 years old", c("time", "dead")], col = "red", cex = 0.4)


matplot(xForInla50$time, res50$summary.random$timeForInla[, 
                                                          c("0.5quant", "0.975quant", 
                                                            "0.025quant")], type = "l", 
        lty = c(1, 2, 2), col = "black", ylim = c(-1, 1) *0.1)

sampleList50 = INLA::inla.posterior.sample(30, res50, selection = list(Predictor = 0)) 
sampleIntensity50 = exp(do.call(cbind, Biobase::subListExtract(sampleList50,
                                                               "latent")))
sampleDeaths50 = matrix(rpois(length(sampleIntensity50),
                              sampleIntensity50), nrow(sampleIntensity50), 
                        ncol(sampleIntensity50))

matplot(xForInla50$time, sampleDeaths50, col = "#00000010", lwd = 2, lty = 1, 
        type = "l", log = "y")
points(x[x$age == "0-49 years old", c("time", "dead")], col = "red", cex = 0.5)
title(main = "Deaths in Population of Age < 50")
matplot(xForInla50$time, sampleDeaths50, col = "#00000010",
        lwd = 2, lty = 1, type = "l", log = "y", xlim = as.Date(c("2019/6/1",
                                                                  "2020/11/1")), 
        ylim = c(1, 2.5) * 45)
points(x[x$age == "0-49 years old", c("time", "dead")], col = "red", cex = 0.5)
title(main = "Sample Deaths in Population of Age < 50")

xPostCovid50 = xPostCovid[xPostCovid$age == "0-49 years old", ]
xPostCovidForecast50 = sampleDeaths50[match(xPostCovid50$time, xForInla50$time), ]
excessDeaths50 = xPostCovid50$dead - xPostCovidForecast50

matplot(xPostCovid50$time, xPostCovidForecast50, type = "l", ylim = c(30, 100), 
        col = "black")
points(xPostCovid50[, c("time", "dead")], col = "red") 

title(main = "Projected Deaths in Population of Age < 50")

matplot(xPostCovid50$time, excessDeaths50, type = "l",
        lty = 1, col = "#00000030")

excessDeathsSub50 = excessDeaths50[xPostCovid50$time > as.Date("2020/03/01") & 
                                     xPostCovid50$time < as.Date("2020/06/01"), ]
excessDeathsInPeriod50 = apply(excessDeathsSub50, 2, sum)
round(quantile(excessDeathsInPeriod50))
round(quantile(excessDeaths50[nrow(excessDeaths50), ]))


