#Load data
smokeFile=load("C:/Users/sa/Documents/MI205/smoke2014.RData")
smoke[1:3, c("Age", "ever_cigarettes", "Sex", "Race", "state", 
             "school", "RuralUrban")]

#Prepare data
forInla = smoke[smoke$Age > 10, c("Age", "ever_cigarettes", 
                                  "Sex", "Race", "state", "school", 
                                  "RuralUrban", "Harm_belief_of_chewing_to")]
forInla = na.omit(forInla)
forInla$y = as.numeric(forInla$ever_cigarettes)
forInla$ageFac = factor(as.numeric(as.character(forInla$Age))) 
forInla$chewingHarm = factor(forInla$Harm_belief_of_chewing_to,levels = 1:4, 
                             labels = c("less", "equal", "more", "dunno"))

#INLA
library("INLA")

toPredict = expand.grid(ageFac = levels(forInla$ageFac), 
                        RuralUrban = levels(forInla$RuralUrban), 
                        chewingHarm = levels(forInla$chewingHarm), 
                        Sex = levels(forInla$Sex))
forLincombs = do.call(inla.make.lincombs, 
                      as.data.frame(model.matrix(~Sex +ageFac * RuralUrban * chewingHarm, 
                                                 data = toPredict))) 
fitS2 = inla(y ~ Sex + ageFac * RuralUrban * chewingHarm + 
               f(state, model = "iid", hyper = list(prec = list(prior = "pc.prec", 
                                                                param = c(99, 0.05)))), 
             data = forInla, family = "binomial",control.inla = list(strategy = "gaussian"), 
             lincomb = forLincombs)

rbind(fitS2$summary.fixed[, c("mean", "0.025quant", "0.975quant")], 
      Pmisc::priorPostSd(fitS2)$summary[, c("mean", "0.025quant", "0.975quant")])
theCoef = exp(fitS2$summary.lincomb.derived[, c("0.5quant", "0.025quant", "0.975quant")])
theCoef = theCoef/(1 + theCoef)

toPredict$Age = as.numeric(as.character(toPredict$ageFac)) 
toPredict$shiftX = as.numeric(toPredict$chewingHarm)/10 
toPredict$x = toPredict$Age + toPredict$shiftX
toPlot = toPredict$Sex == "M" & toPredict$RuralUrban == "Rural"

plot(toPredict[toPlot, "x"], theCoef[toPlot, "0.5quant"], xlab = "age", 
     ylab = "probability", ylim = c(0,1), pch = 15, col = toPredict[toPlot, "chewingHarm"]) 

segments(toPredict[toPlot, "x"], theCoef[toPlot, "0.025quant"],
         y1 = theCoef[toPlot, "0.975quant"], col = toPredict[toPlot, "chewingHarm"])

legend("topleft", fill = 1:nlevels(toPredict$chewingHarm), 
       legend = levels(toPredict$chewingHarm), bty = "n", title = "chewing harmful")

