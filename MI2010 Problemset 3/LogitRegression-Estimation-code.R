mydata <- read.csv("https://stats.idre.ucla.edu/stat/data/binary.csv")

## Standard Logistic Regression
mylogit<-glm(admit ~ gre + gpa +
               as.factor(rank), data=mydata, family="binomial")
summary(mylogit)

## log(p/(1-p)) = -3.99 + 0.002*gre + 0.804*gpa-0.675*x3
##                           -1.34*x4 - 1.55*x5


## Survey Estimation for Logistic Regression
n=length(mydata$admit)
N=6000

#install.packages("survey")
library(survey)
## Using the Survey Library
fpc.srs = rep(N, n)

ucla.design <- svydesign(id=~1, data=mydata, fpc=fpc.srs)

mysvyglm <- svyglm(admit ~ gre + gpa + as.factor(rank), 
                   ucla.design, family="binomial")
summary(mysvyglm)
