library(HDInterval)
library(tidyverse)

## compare whether two posterior distributions to see if they are sigf diff from each other
## adapted from BEST method developed by Kruschke 2013
########################################################
## load data ##
## flyVac data - mean of posterior dist of parameters
setwd("~/de Bivort 2/Enrichment Paper/Stan Output/CS Means and Vars")
load('lightTime_CSposterior.RData')

#switch factor levels to Low Med High for CS data
post_dist_mu$trt = post_dist_mu$trt %>% recode(a='Low',b1='Med',b2='High')
post_dist_var$trt = post_dist_var$trt %>% recode(v0='Low',v1='Med',v2='High')

#create coeff of var (sd/mean) tibble
post_dist_cv = tibble(trt = post_dist_mu$trt, cv = sqrt(post_dist_var$var)/post_dist_mu$mu)


# calc HDI of differences #
calcTwoSampleDiff = function(post_dist, trt1, trt2){
y1 = post_dist %>% 
    filter(trt == trt1)
  
y2 = post_dist %>% 
    filter(trt == trt2)

HDinterval = as.numeric(hdi(y1[,2]-y2[,2], credMass = 0.99))
percBelowZero = sum((y1[,2]-y2[,2])<=0)/dim(y1)[1]
percAboveZero = sum((y1[,2]-y2[,2])>=0)/dim(y1)[1]
return(c(HDinterval,percBelowZero,percAboveZero))
} #99%HDI

calcTwoSampleDiff(post_dist_mu, 'Low','High')
calcTwoSampleDiff(post_dist_var, 'Med','High')
calcTwoSampleDiff(post_dist_cv, 'Low','Med')

########################################################

## test whether 0 falls in the 99% HDI of the posterior distribution of the effect size
########################################################
#data - use Stan_Enrichment_PaperFigures to modify the contrasts and do CV coeffs calc prior to running this analysis

mu.HDInterval = hdi(post_coef_mu, credMass = 0.99)
var.HDInterval = hdi(post_coef_var, credMass = 0.99)
cv.HDInterval = hdi(post_coef_cv, credMass = 0.99)

########################################################
