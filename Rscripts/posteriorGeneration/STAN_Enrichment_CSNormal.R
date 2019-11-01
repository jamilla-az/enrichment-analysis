##############################
## Enrichment STAN analysis ##
##############################

library(tidyverse)
library(rstan)
library(bayesplot)
library(gridExtra)
library(shinystan)

# CS only #
# turn bias, clump, switch, light pref #

###############
## load data ##
###############
## set the working directory
setwd("~/de Bivort 2/Enrichment Paper/Stan Input")

##CS only ##
#outcome vector
y = read.csv('turnBias_CS_Ymaze.csv', header = F)
y = as.numeric(y$V1)

#predictor matrix
x = read.csv('predMatrix_CS_Ymaze_all.csv', header = F)
X = data.matrix(x)
X = X[,c(3:4)]

# stan data input
pheno_dat = list(N = length(y),
                    L = dim(X)[2],
                    y, X)
##############
## run STAN ##
##############

setwd("~/de Bivort 2/Enrichment Paper/stan scripts")

# CS only
fit <- stan(file = 'CS_NormalModel.stan', data = pheno_dat, 
            iter = 100000, chains = 4, control = list(adapt_delta = 0.8, max_treedepth = 10))

#launch shiny stan
my_sso <- launch_shinystan(fit)

########################
## extract raw values ##
########################

post_dist_raw = as.matrix(fit, pars = c('a', 'be', 'v0','ve')) 

#CS only
post_dist_mu = post_dist_raw %>% 
  as_tibble %>%
  mutate(b1 = a+`be[1]`, b2 = a+`be[2]`) %>%
  select(a,b1,b2) %>%
  gather(key = trt, value = mu)

post_dist_var = post_dist_raw %>% 
  as_tibble %>%
  mutate(v1 = v0+`ve[1]`, v2 = v0+`ve[2]`) %>%
  select(v0,v1,v2) %>%
  gather(key = trt, value = var)

#extract simulated response
sim_turnBias = as.matrix(fit, pars = c("y_rep"))

#save data
setwd("~/de Bivort 2/Enrichment Paper/Stan Output")
save(post_dist_mu, post_dist_var, file = 'turnBias_CSposterior.RData')


#################
## DIAGNOSTICS ##
#################

## comparing original dist to simulated dist
library(bayesplot)
library(gridExtra)

#density plots
## set the working directory
setwd("~/de Bivort 2/Enrichment Paper/Stan Input")
#predictor matrix
x = read.csv('predMatrix_CS_YMaze_all.csv', header = F)
X = data.matrix(x)

## CS ONLY
color_scheme_set("brightblue")
p1=ppc_dens_overlay(y[as.logical(X[,2])], sim_turnBias[1:50,as.logical(X[,2])])
p2=ppc_dens_overlay(y[as.logical(X[,3])], sim_turnBias[1:50,as.logical(X[,3])])
p3=ppc_dens_overlay(y[as.logical(X[,4])], sim_turnBias[1:50,as.logical(X[,4])])

grid.arrange(p1, p2, p3, nrow=1)

trt = as.factor(c(rep('low', sum(as.logical(X[,2]))),
                  rep('med', sum(as.logical(X[,3]))),
                  rep('high', sum(as.logical(X[,4])))))
trt = relevel(trt,c('low'))
ppc_stat_grouped(y, sim_turnBias, group = trt, stat='sd')


