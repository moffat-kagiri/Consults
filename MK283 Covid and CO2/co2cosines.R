
#+ header, result='asis', echo=FALSE
knitr::knit_hooks$set(
    marginsp = function(before, options, envir){  
      if(!before) return()
# use small margins       
      par(mar=c(1.5+0.9*options$marginsp,
              1.5+0.9*options$marginsp,0.2,0.2),
          mgp=c(1.45, 0.45, 0),cex=1)
    }   
)
knitr::opts_chunk$set(
  fig.height=4, fig.width=4, echo=TRUE, marginsp=TRUE, out.width=Pmisc::out.width(0.9))
  knitr::knit_hooks$set(plot=Pmisc::hook_plot_mdsubfig)
  knitr::opts_chunk$set(dev='png')
Pmisc::markdownHeader('CO2 and seasonal sines and cosines')

#' 
#' 
#' 
#' 

#' Mauna Loa Observatory, Hawaii
#' [scrippsco2.ucsd.edu/data/atmospheric_co2/mlo.html](https://scrippsco2.ucsd.edu/data/atmospheric_co2/mlo.html)

#+ co2Data
cUrl = paste0('https://scrippsco2.ucsd.edu/assets/data/atmospheric/',
	'stations/flask_co2/daily/daily_flask_co2_mlo.csv')
cFile = basename(cUrl)
if(!file.exists(cFile)) download.file(cUrl, cFile)
co2s = read.table(cFile, header=FALSE, sep=',', skip=69, stringsAsFactors=FALSE,
  col.names=c('dateOrig','time','junk1','junk2', 'Nflasks','quality','co2'))
co2s$date = as.Date(co2s$dateOrig)
# remove low-quality measurements
co2s = co2s[co2s$quality==0,  ]
#'
#+ showData
co2s[1:4, ]
#' 

#+ co2plot, fig.height=4, fig.width=4, fig.cap='CO2 at Mauna Loa Observatory, Hawaii', fig.pos='htbp', fig.subcap=c('all','recent'), fig.ncol=2, out.width=Pmisc::out.width(0.4)
plot(co2s$date, co2s$co2, log='y', cex=0.3, col='#00000040', xlab='time', ylab='ppm')
plot(co2s[co2s$date > as.Date("2019/1/1"), c('date','co2')], log='y', type='o', xlab='time', ylab='ppm', cex=0.5)
#'
#' 
#' create a numeric time variable
#' number of days since 1 Jan 2000
#' then divide by 365.25, the number of days it takes for the earth to go around the sun.
#' to get a time variable in years
#+ setupCo2, echo=TRUE
timeOrigin = as.Date("2000/1/1")
co2s$timeNumeric = as.numeric(co2s$date - timeOrigin)
co2s$timeYears = co2s$timeNumeric/365.25
#'



#' # Modelling seasonal effects
#'
#'
#' - We'll want a model like
#' 
#' $$
#' \begin{aligned}
#' Y_i \sim & \operatorname{Poisson}[\lambda(t_i)]\\
#' \log[\lambda(t)] = & \mu + s(t) + f(t) \\
#' \end{algined}
#' $$
#' 
#' - $s(t)$ is a seasonal cycle
#' - $f(t)$ is a smoothly varying function
#' $$
#' - use months as indicator variables?
#' - $s(t) = \theta_j$ if $t$ is in month $j$?
#' - abrupt changes on the first of every month.
#' - The annual cycles looks sinusoidal
#' - or $s(t)= \rho \cos(2 \pi t/365.25 + \phi)$?
#' - the later has 2 parameters instead of 12.
#' - ... and is smooth
#' - $\rho$ is the amplitude
#' - $\phi$ is the phase
#' - $1/365.25$ is the frequency, one cycle every 365.25 days
#' 
#' A trick:
#' $$
#' \beta_1 \cos(2 \pi t/365.25) + \beta_2 \sin(2 \pi t/365.25) = 
#' \sqrt{\beta_1 + \beta_2}
#' \cos\left[2 \pi t/365.25 + 
#' \operatorname{arctan} \left(-{\frac {\beta_2}{\beta_1}}\right)
#' \right]
#' $$
#' 
#' - use $X_{i1} = \cos(2 \pi t_i/365.25)$ and $X_{i2}(t) = \sin(2 \pi t_i/365.25)$ as linear covariates
#' - where $t_i$ is the time (in months) of observation $i$
#'


#' # First harmonics

#' - The monthly effect isn't perfectly sinusoidal
#' - use a 12 month and a 6 month frequency
#' - these four sinusoidal basis functions can capture a wide range of seasonal effects

#' \begin{align*}
#' Y_i \sim & \text{Poisson}(O_i \lambda_i)\\
#' \log(\lambda_i) = & X_i \beta + f(t_i)\\
#' X_{i0} = & 1\\
#' X_{i1} = & \cos(2 \pi t_i/365.25)\\
#' X_{i2} = & \sin(2 \pi t_i/365.25)\\
#' X_{i3} = & \cos(2 \pi t_i/182.625)\\
#' X_{i4} = & \sin(2 \pi t_i/182.625)
#' \end{align*}

#+ monthHarmonicPlot, echo=FALSE, fig.cap='first harmonic plot', fig.ncol=1
set.seed(0)
Stime = seq(-18,18, len=1001)
Xmat = cbind(cos(2*pi*Stime/12),sin(2*pi*Stime/12),cos(2*pi*Stime/6),sin(2*pi*Stime/6))
coefMat = matrix(rnorm(ncol(Xmat)*3), nrow=ncol(Xmat))
matplot(Stime, Xmat %*% coefMat, type='l', lty=1, xaxs='i', xlab='time', ylab='X beta')
#'

#+ coefMatHarmonic,echo=FALSE, warning=FALSE
rownames(coefMat) =c('cos12','sin12','cos6','sin6')
knitr::kable(coefMat, digits=2)
#'
#' # Add seasonal term to co2 model
#' 
#' Create 12 month cycles
#+ createSines, echo=TRUE
co2s$sin12 = sin(2*pi*co2s$timeYears)
co2s$cos12 = cos(2*pi*co2s$timeYears)
#'
#' add 6 month cycles
#+ createSinesHarmoic, echo=TRUE
co2s$sin6 = sin(2*pi*co2s$timeYears*2)
co2s$cos6 = cos(2*pi*co2s$timeYears*2)
#'
#' log transform 
#+ logCo2, echo=TRUE
co2s$logCo2 = log(co2s$co2)
#'
#' fit model with $s(t) =$ sum of 4 sinusoids
#+ co2gam, echo=TRUE, cache=TRUE
res = lm(logCo2 ~ timeNumeric + sin12 + cos12 + 
                  sin6 + cos6, data=co2s)
#'
#+ modelparameters
knitr::kable(summary(res)$coef, digits=5)
#'
#' create new data to predict seasonal effect
#+ newData, echo=TRUE
StimeYears = seq(-0.1, 1.1, len=1001)
newData = cbind(
 sin12 = sin(2*pi*StimeYears),
 cos12 = cos(2*pi*StimeYears),
 sin6 = sin(2*2*pi*StimeYears),
 cos6 = cos(2*2*pi*StimeYears))
#' 
#' feed it to the predict function
#' ask for only the smooth effect, not the sinusoids
#' and create 95% prediction interval
#+ plotTrend, cache=TRUE
res2 = cbind(
	est = newData %*% res$coef[colnames(newData)],
	se.fit = sqrt(diag(newData %*% 
		vcov(res)[
			colnames(newData), colnames(newData)
			] %*% t(newData)))
	) 

res2ci = res2 %*% Pmisc::ciMat(0.95)
#'
#' 
#' plot predictions
#+ predPlot, fig.cap='predicted trend', fig.subcap=' ',fig.ncol=1, echo=TRUE
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
#' 
#' now with INLA
#+ cyclesInla
library(INLA)
co2s$timeInla = round(co2s$timeYears, 2)
# disable some error checking in INLA
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
#' Predict seasonal effect
#' with 99.9% interval
#+ predSeason

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
#'



#' plot seasonal predictions
#+ plotSeason, fig.cap='seasonal predictions', fig.ncol=1 
qcols = paste0(c('0.025','0.5','0.975'), 'quant')
matplot(forX, 
	exp(resInla2$summary.lincomb.derived[,qcols]), 
	type='l', xlab='', ylab='relative co2',
	xaxt='n', lty=1, col=c('black','grey','grey'),
	xaxs='i')
axis(1,as.numeric(forXlab), format(forXlab, '%b'))
#' 

