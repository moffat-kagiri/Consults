

dataDir = '../data' 
smokeFile = file.path(dataDir, 'smoke2014.RData') 
if(!file.exists(smokeFile)){  
  download.file(
      'http://pbrown.ca/teaching/astwo/data/smoke2014.RData',
      smokeFile) 
} 
load(smokeFile) 
smoke[1:3,c('Age','ever_cigarettes','Sex','Race', 
        'state','school', 'RuralUrban')]




forInla = smoke[,c('Age','ever_cigarettes','Sex','Race', 
        'state','school', 'RuralUrban')]
forInla = na.omit(forInla)
forInla$y = as.numeric(forInla$ever_cigarettes)
forInla$ageFac = relevel(factor(forInla$Age), '14')

  toPredict = expand.grid(
    ageFac = levels(forInla$ageFac),
    RuralUrban = levels(forInla$RuralUrban),
    Sex = levels(forInla$Sex)
    )
forLincombs = do.call(inla.make.lincombs, 
  as.data.frame(model.matrix( ~ ageFac:RuralUrban + Sex, 
    data=toPredict)))

library("INLA")



fitS2 = inla(y ~ Sex + ageFac:RuralUrban + 
    f(state, model='iid', hyper=list(
        prec=list(prior='pc.prec', param=c(99, 0.05)))
    ),
  data=forInla, family='binomial',
  lincomb = forLincombs)

dim(toPredict)
dim(fitS2$summary.lincomb.derived)
toPredict[1:2,]
fitS2$summary.lincomb.derived[1:2,]

rbind(
    fitS2$summary.fixed[, c('mean','0.025quant','0.975quant')],
    Pmisc::priorPostSd(fitS2)$summary[, c('mean','0.025quant','0.975quant')]
  )


theCoef = exp(fitS2$summary.lincomb.derived[,
    c('0.5quant','0.025quant','0.975quant')])
theCoef = theCoef/(1+theCoef)
theSd = fitS2$summary.lincomb.derived[,'sd']

toPredict$Age = as.numeric(as.character(toPredict$ageFac))

isMale = toPredict$Sex == 'M'
shiftRural = 0.1*(toPredict$RuralUrban == 'Rural')

theCex = min(theSd)/theSd
plot(toPredict[isMale,'Age'] + shiftRural[isMale], 
  theCoef[isMale,'0.5quant'], 
  xlab='age', ylab='probability', ylim = c(0.015, 0.7),
  pch = 15, log='y', 
  cex = 2*theCex,
  col = mapmisc::col2html(
    c(Urban = 'red', Rural = 'green')[as.character(toPredict[isMale,'RuralUrban'])],
    0.4)
  )


segments(toPredict[isMale,'Age']+ shiftRural[isMale], 
  theCoef[isMale,'0.025quant'], 
  y1=theCoef[isMale,'0.975quant'],
  col = c(Urban = 'red', Rural = 'green')[as.character(toPredict[isMale,'RuralUrban'])])

legend('bottomright', pch=16, col=c('red','green'), legend = c('Urban','Rural'),
  bty='n')
