library("tidyverse")
library("janitor")
library("INLA")
library("ggplot2")

cUrl = paste0('https://scrippsco2.ucsd.edu/assets/data/atmospheric/',
              'stations/flask_co2/daily/daily_flask_co2_mlo.csv')
cFile = basename(cUrl)
if(!file.exists(cFile)) download.file(cUrl, cFile)
co2s = read.table(cFile, header=FALSE, sep=',', skip=69, stringsAsFactors=FALSE,
                  col.names=c('dateOrig','time','junk1','junk2', 'Nflasks',
                              'quality','co2'))
co2s$date = as.Date(co2s$dateOrig)
# remove low-quality measurements
co2s = co2s[co2s$quality==0,  ]

co2s[1:4, ]

plot(co2s$date, co2s$co2, log='y', cex=0.3, col='#00000040', xlab='time', ylab='ppm')
plot(co2s[co2s$date > as.Date("2019/1/1"), c('date','co2')], log='y', type='o', 
     xlab='time', ylab='ppm', cex=0.5)


timeOrigin = as.Date("1984/1/1")
co2s$timeNumeric = as.numeric(co2s$date - timeOrigin)
co2s$timeYears = co2s$timeNumeric/365.25

co2s$sin12 = sin(2*pi*co2s$timeYears)
co2s$cos12 = cos(2*pi*co2s$timeYears)

co2s$sin6 = sin(2*pi*co2s$timeYears*2)
co2s$cos6 = cos(2*pi*co2s$timeYears*2)

co2s$logCo2 = log(co2s$co2)

res = lm(logCo2 ~ timeNumeric + sin12 + cos12 + 
           sin6 + cos6, data=co2s)

knitr::kable(summary(res)$coef, digits=5)

StimeYears = seq(-0.1, 1.1, len=1001)
newData = cbind(
  sin12 = sin(2*pi*StimeYears),
  cos12 = cos(2*pi*StimeYears),
  sin6 = sin(2*2*pi*StimeYears),
  cos6 = cos(2*2*pi*StimeYears))

res2 = cbind(
  est = newData %*% res$coef[colnames(newData)],
  se.fit = sqrt(diag(newData %*% 
                       vcov(res)[
                         colnames(newData), colnames(newData)
                       ] %*% t(newData)))
) 

res2ci = res2 %*% Pmisc::ciMat(0.95)

forX = as.Date('2020/1/1') + 365.25*StimeYears
matplot(forX,
        exp(res2ci), 
        xlab='time', ylab='co2 trend',
        type='l', xaxt='n',
        xaxs='i', log='y',
        col='black', 
        lty=c(1,2,2))
forXlab = seq(as.Date('2020/1/1'), len=12, by='3 months')
axis(1,as.numeric(forXlab), format(forXlab, '%b'))



co2s$timeInla = round(co2s$timeYears, 2)

berlinwall<-
  co2s %>%
  filter(date >= as.Date("1987/01/01") & date <= as.Date("1991/01/01"))


mm = get("inla.models", INLA:::inla.get.inlaEnv())
if(class(mm) == 'function') mm = mm()
mm$latent$rw2$min.diff = NULL
assign("inla.models", mm, INLA:::inla.get.inlaEnv())

resInla = inla(logCo2 ~ f(timeInla, model='rw2',
                          prior = 'pc.prec', param = c(0.1, 0.5)) + 
                 sin12 + cos12 + 
                 sin6 + cos6, data=co2s,
               control.family = list(hyper=list(prec=list(
                 prior='pc.prec', param = c(0.1,0.5)
               )))
)
Pmisc::priorPost(resInla)$summary[,c(4,3,5)]

newDataSeason = data.frame(
  cos12 = cos(2*pi*StimeYears),
  sin12 = sin(2*pi*StimeYears),
  cos6 = cos(2*2*pi*StimeYears),
  sin6 = sin(2*2*pi*StimeYears))
forLincombs = do.call(
  inla.make.lincombs,
  newDataSeason)
resInla2 = inla(resInla$.args$formula, 
                data=co2s,
                control.family = resInla$.args$control.family,
                lincomb = forLincombs)
summary(resInla2)

resInlab <- inla(resInla$.args$formula, 
                data=berlinwall,
                control.family = resInla$.args$control.family,
                lincomb = forLincombs)
summary(resInlab)

plot(berlinwall$date, berlinwall$co2, log='y', cex=0.3, col='#00000040', 
     xlab='time', ylab='ppm')
title(main = "Berlin Wall Effect")
lines(fitted(resInlab))

qcols = paste0(c('0.025','0.5','0.975'), 'quant')
matplot(forX, 
        exp(resInla2$summary.lincomb.derived[,qcols]), 
        type='l', xlab='', ylab='relative co2',
        xaxt='n', lty=1, col=c('black','grey','grey'),
        xaxs='i')
axis(1,as.numeric(forXlab), format(forXlab, '%b'))

coronavirus<-
  co2s %>%
  filter(date >= as.Date("2000/01/01") & date <= as.Date("2025/01/01"))


co2s$estimate <-
  predict(resInla2, newdata = coronavirus, type = "response")