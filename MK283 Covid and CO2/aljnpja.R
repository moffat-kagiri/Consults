
co2s$day = as.Date(co2s$date)
toAdd = data.frame(day = seq(max(co2s$day) + 3, as.Date("2025/1/1"), 
                             by = "10 days"), co2 = NA)

co2ext = rbind(co2s[,colnames(toAdd)], toAdd) 
timeOrigin = as.Date("2000/1/1")
co2ext$timeInla = round(as.numeric(co2ext$day - timeOrigin)/365.25, 2)
co2ext$cos12 = cos(2 * pi * co2ext$timeInla) 
co2ext$sin12 = sin(2 * pi * co2ext$timeInla) 
co2ext$cos6 = cos(2 * 2 * pi * co2ext$timeInla) 
co2ext$sin6 = sin(2 * 2 * pi * co2ext$timeInla)

mm = get("inla.models", INLA:::inla.get.inlaEnv()) 
if(class(mm) == 'function') mm = mm() 
mm$latent$rw2$min.diff = NULL
assign("inla.models", mm, INLA:::inla.get.inlaEnv())

co2res = inla(co2 ~ sin12 + cos12 + sin6 + cos6 + f(timeInla, model = 'rw2',
                                                    prior='pc.prec', 
                                                    param = c(0.1, 0.5)), 
              data = co2ext, family='gamma', 
              control.family = list(hyper=list(
                prec=list(prior='pc.prec', param=c(0.1, 0.5)))),
              # add this line if your computer has trouble 
              #	control.inla = list(strategy='gaussian'),
              control.predictor = list(compute=TRUE, link=1), 
              control.compute = list(config=TRUE), verbose=FALSE)
qCols = c('0.5quant','0.025quant','0.975quant') 
Pmisc::priorPost(co2res)$summary[,qCols]

sampleList = INLA::inla.posterior.sample(30, co2res, selection = list(timeInla = 0))
sampleMean = do.call(cbind, Biobase::subListExtract(sampleList, "latent"))
sampleDeriv = apply(sampleMean, 2, diff)/diff(co2res$summary.random$timeInla$ID)



matplot(co2ext$day, co2res$summary.fitted.values[, qCols], type = "l", 
        col = "black", lty = c(1, 2, 2), log = "y", xlab = "time", ylab = "ppm")
Stime = timeOrigin + round(365.25 * co2res$summary.random$timeInla$ID) 

matplot(Stime, co2res$summary.random$timeInla[, qCols],type = "l", 
        col = "black", lty = c(1, 2, 2), xlab = "time", ylab = "y")
matplot(Stime[-1], sampleDeriv, type = "l", lty = 1,
        xaxs = "i", col = "#00000020", xlab = "time", ylab = "deriv", 
        ylim = quantile(sampleDeriv, c(0.01, 0.995)))
forX = as.Date(c("2018/1/1", "2021/1/1")) 
forX = seq(forX[1], forX[2], by = "6 months")
toPlot = which(Stime > min(forX) & Stime < max(forX))

matplot(Stime[toPlot], sampleDeriv[toPlot, ], type = "l", lty = 1, lwd = 2, xaxs = "i", 
        col = "#00000050", xlab = "time", ylab = "deriv", xaxt = "n", 
        ylim = quantile(sampleDeriv[toPlot,], c(0.01, 0.995)))

axis(1, as.numeric(forX), format(forX, " b Y"))