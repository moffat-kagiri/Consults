library('INLA', verbose=FALSE)
library('Biobase')
library('Pmisc')

cUrl = paste0('https://scrippsco2.ucsd.edu/assets/data/atmospheric/',
              'stations/flask_co2/daily/daily_flask_co2_mlo.csv')
cFile = basename(cUrl)
if(!file.exists(cFile)) download.file(cUrl, cFile)

co2s = read.table(cFile,header = FALSE, sep = ",",skip = 69, stringsAsFactors = FALSE, 
                  col.names = c("day", "time", "junk1", "junk2", "Nflasks", "quality", 
                                "co2"))
co2s$date = strptime(paste(co2s$day, co2s$time), format = "%Y-%m-%d %H:%M", tz = "UTC")

# remove low-quality measurements
co2s = co2s[co2s$quality == 0, ]
plot(co2s$date, co2s$co2, log = "y", cex = 0.3, 
     col = "#00000040", xlab = "time", ylab = "ppm")
plot(co2s[co2s$date > ISOdate(2015, 3, 1, tz = "UTC"), 
          c("date", "co2")], log = "y", type = "o", xlab = "time", 
     ylab = "ppm", cex = 0.5)
title(main = "Progress since 2015")


co2s$day = as.Date(co2s$date)
toAdd = data.frame(day = seq(max(co2s$day) + 3, 
                             as.Date("2025/1/1"),by = "10 days"), co2 = NA)

co2ext = rbind(co2s[, colnames(toAdd)], toAdd)

# Create time variables to use

derivLincomb = inla.make.lincombs(timeRw2 = D[-1, ])
names(derivLincomb) = gsub("^lc", "time", names(derivLincomb))

timeOrigin = ISOdate(1980, 1, 1, 0, 0, 0, tz = "UTC")
co2s$days = as.numeric(difftime(co2s$date, timeOrigin, units = "days"))

co2ext$timeInla = round(as.numeric(co2ext$day - timeOrigin)/365.25,2)
co2ext$cos12 = cos(2 * pi * co2ext$timeInla) 
co2ext$sin12 = sin(2 * pi * co2ext$timeInla)
co2ext$cos6 = cos(2 * 2 * pi * co2ext$timeInla) 
co2ext$sin6 = sin(2 * 2 * pi * co2ext$timeInla)

#Predictions
timeBreaks = seq(min(co2s$date), ISOdate(2025, 1, 1, tz = "UTC"), by = "14 days")
timePoints = timeBreaks[-1]
co2s$timeRw2 = as.numeric(cut(co2s$date, timeBreaks))

timeBreaks1 = seq(min(co2s$date), ISOdate(2030, 1, 1, tz = "UTC"), by = "14 days")
timePoints1 = timeBreaks1[-1]
StimePred = as.numeric(difftime(timePoints1, timeOrigin, units = "days"))/365.35
predLincomb = inla.make.lincombs(timeRw2 = Diagonal(length(timePoints1)), 
                                 `(Intercept)` = rep(1, length(timePoints1)), 
                                 sin12 = sin(2* pi * StimePred), 
                                 cos12 = cos(2 * pi * StimePred), 
                                 sin6 = sin(2 * 2 * pi * StimePred), 
                                 cos6 = cos(2 * 2 * pi * StimePred))
names(predLincomb) = gsub("^lc", "pred", names(predLincomb))
StimeIndex = seq(1, length(timePoints1))

# disable some error checking in INLA
mm = get("inla.models", INLA:::inla.get.inlaEnv()) 
if(class(mm) == 'function') mm = mm() 
mm$latent$rw2$min.diff = NULL
assign("inla.models", mm, INLA:::inla.get.inlaEnv())
co2res = inla(co2 ~ sin12 + cos12 + sin6 + cos6 + f(timeInla, model = 'rw2',
                                                    prior='pc.prec', 
                                                    param = c(0.1, 0.5)), data = co2ext, 
              family='gamma', control.family = list(hyper=list(prec=list(
                                                      prior='pc.prec', 
                                                      param=c(0.1, 0.5)))),
              control.inla = list(strategy='gaussian'), 
              lincomb = c(derivLincomb,seasonLincomb, predLincomb), 
              control.predictor = list(compute=TRUE, link=1), 
              control.compute=list(config = TRUE), 
              verbose=FALSE)
qCols = c('0.5quant','0.025quant','0.975quant') 
co2resTab1 <- exp(co2res$summary.fixed[, c("mean", "0.025quant", "0.975quant")])

Pmisc::priorPost(co2res)$summary[,qCols]


sampleList = INLA::inla.posterior.sample(30, co2res, selection = list(timeInla = 0))
sampleMean = do.call(cbind, Biobase::subListExtract(sampleList, "latent"))
sampleDeriv = apply(sampleMean, 2, diff)/diff(co2res$summary.random$timeInla$ID)
matplot(co2ext$day, co2res$summary.fitted.values[, qCols], type = "l", col = "black", 
        lty = c(1, 2, 2), log = "y", xlab = "time", ylab = "ppm")
Stime = timeOrigin + round(365.25 * co2res$summary.random$timeInla$ID) 
matplot(Stime, co2res$summary.random$timeInla[, qCols],type = "l", col = "black", 
        lty = c(1, 2, 2), xlab = "time", ylab = "y")

matplot(Stime[-1], sampleDeriv, type = "l", lty= 1, xaxs = "i", 
        col = "#00000020", xlab = "time", ylim = quantile(sampleDeriv, c(0.01, 0.995)))



forX = as.Date(c("2018/1/1", "2021/1/1"))
forX = seq(forX[1], forX[2], by = "6 months")
toPlot = which(Stime > min(forX) & Stime < max(forX))

matplot(Stime[toPlot], sampleDeriv[toPlot, ], type = "l",
        lty = 1, lwd = 2, xaxs = "i", col = "#00000050",
        xlab = "time", ylab = "deriv", xaxt = "n", 
        ylim = quantile(sampleDeriv[toPlot,], c(0.01, 0.995)))
axis(1, as.numeric(forX), format(forX, "%b%Y"))


derivPred = co2res$summary.lincomb.derived[grep("time", 
                                                rownames(co2res$summary.lincomb.derived)), 
                                           c("0.5quant", "0.025quant","0.975quant")]
scaleTo10Years = (10 * 365.25/as.numeric(diff(timePoints, units = "days")))
matplot(timePoints[-1], scaleTo10Years * derivPred, type = "l", col = "black", 
        lty = c(1, 2, 2), ylim = c(0, 0.1), xlim = range(as.numeric(co2s$date)), 
        xaxs = "i", xaxt = "n", xlab = "time", ylab = "log ppm, change per 10yr")
axis(1, xax, format(xax, "%Y"))
abline(v = ISOdate(1989, 11, 1, tz = "UTC"), col = "red")
abline(v = ISOdate(1992, 1, 1, tz = "UTC"), col = "red")
abline(v = ISOdate(2021, 1, 1, tz = "UTC"), col = "purple")
abline(v = ISOdate(2020, 1, 1, tz = "UTC"), col = "purple")
legend("topright", legend = c("1989-91 USSR Collapse", "2020+ COVID-19 pandemic"), 
       col = c("red","purple"))







