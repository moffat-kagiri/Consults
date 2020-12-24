library(INLAutils)
library(INLA)
library(ggplot2)
library(tidyverse)
#load data

sUrl = "http://www.bristol.ac.uk/cmm/media/migrated/jsp.zip" 
dir.create(file.path("..", "data"), showWarnings = FALSE) 
(Pmisc::downloadIfOld(sUrl, file.path("..", "data")))

#create dataset from the info

school = read.fwf("../data/JSP.DAT", widths = c(2, 1, 1, 1, 2, 4, 2, 2, 1), 
                  col.names = c("school", "class", "gender", "socialClass", 
                                "ravensTest", "student", "english", "math", "year"))

#variables

school$socialClass = factor(school$socialClass, 
                            labels = c("I", "II", "IIIn", "IIIm", "IV", "V", 
                                       "longUnemp", "currUnemp", "absent"))
school$gender = factor(school$gender, labels = c("f", "m"))
school$classUnique = paste(school$school, school$class) 
school$studentUnique = paste(school$school, school$class,school$student)
school$grade = factor(school$year)

#generalized linear model

schoolLme = glmmTMB::glmmTMB(math ~ gender + socialClass + grade + (1 | school) + 
                               (1 | classUnique) + (1 | studentUnique), data = school)
summary(schoolLme)
knitr::kable(summary(schoolLme)$coef,digits = 3,caption = 'Regression Result')

#histogram

hist(1 - school$math,  breaks = 100)

#INLA

prec.prior <- list(prec = list(param = c(30, 0.05)))

mathscore = INLA::inla(math ~ gender + socialClass + grade+ f(
  school, model='iid', 
  hyper = prec.prior) +f(
    school, model='iid', 
    hyper = prec.prior), 
  data = school, control.predictor = list(compute = TRUE))

summary(mathscore)
knitr::kable(mathscore$summary.fixed, digits = 2, caption = "Posterior Quantiles")

#Plots of the original data
genderplot <- ggplot(school, aes(x= math, fill=gender, color=gender)) +
  geom_histogram(position="identity", binwidth=1, alpha=0.5) + labs(title = "Gender")

soclassplot <- ggplot(school, aes(x= math, fill=socialClass, color=socialClass)) +
  geom_histogram(position="identity", binwidth=1, alpha=0.5) + labs(title = "Socialclass")

gradeplot <- ggplot(school, aes(x= math, fill=grade, color=grade)) +
  geom_histogram(position="identity", binwidth=1, alpha=0.5) + labs(title = "Grade")

genderplot

soclassplot

gradeplot

#PLot of prior and posterior

sdRes = Pmisc::priorPostSd(mathscore)
do.call(matplot, sdRes$matplot)
do.call(legend, sdRes$legend)


