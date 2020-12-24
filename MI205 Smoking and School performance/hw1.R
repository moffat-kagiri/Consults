library(statsr)
df = data('Affairs',package = 'AER')

Affairs$ever = Affairs$affair > 0
Affairs$religious = factor(Affairs$religiousness,
                           levels = c(2,1,3,4,5), labels = c('no','anti','low','med','high'))

Affairs$age_normalized = (Affairs$age - mean(Affairs$age))/sd(Affairs$age)
Affairs$yearsmarried_normalized = (Affairs$yearsmarried - mean(Affairs$yearsmarried))/sd(Affairs$yearsmarried)
Affairs$education_normalized = (Affairs$education - mean(Affairs$education))/sd(Affairs$education)

res = glm(ever ~ gender:children + age_normalized + yearsmarried_normalized + religious + education_normalized,
          data = Affairs, family = binomial(link = 'logit'))

knitr::kable(summary(res)$coef,digits = 3,caption = 'Regression Result')
