##############################
## Enrichment STAN analysis ##
##############################
library(tidyverse)
library(rstan)
library(shinystan)
library(bayesplot)
library(gridExtra)

##############
# functions #
##############
calcDeviations3Geno = function(post_dist, labels){
  
  post_coef = matrix(NA, dim(post_dist)[1], dim(post_dist)[2]+2)
  post_coef[,1] = apply(post_dist, 1, mean) #grand param val
  
  #dev of each geno from grand param value
  post_coef[,2] = apply(post_dist[,c(1,4)], 1, mean) - post_coef[,1]#dev of G0
  post_coef[,3] = apply(post_dist[,c(2,5)], 1, mean) - post_coef[,1]#dev of G1
  post_coef[,4] = apply(post_dist[,c(3,6)], 1, mean) - post_coef[,1]#dev of G3
  
  #dev of mild enrich from grand mean
  post_coef[,5] = apply(post_dist[,c(4:6)], 1, mean) - post_coef[,1]#dev of MildEnr
  
  #dev of GenoxMildEnr
  post_coef[,6] = post_dist[,4] - post_coef[,1] - post_coef[,2]#dev of G0xMildEnr
  post_coef[,7] = post_dist[,5] - post_coef[,1] - post_coef[,3]#dev of G1xMildEnr
  post_coef[,8] = post_dist[,6] - post_coef[,1] - post_coef[,4]#dev of G2xMildEnr
  
  colnames(post_coef) = labels
  return(post_coef)
}

###############
## load data ##
###############
## set the working directory
setwd("~/de Bivort 2/Enrichment Paper/Stan Input")

##CS only ##

#outcome vector
y = read.csv('intertrialTime_CS_FlyVac.csv', header = F)
y = as.numeric(y$V1)

#predictor matrix
x = read.csv('predMatrix_CS_FlyVac_all.csv', header = F)
X = data.matrix(x[,2:4])

## num turns for all genotypes and trt combos
data = read.csv('allGenoTrt_FlyVacTime.csv', header = F)
colnames(data) = c('geno','trt','pref')
data$geno = as.factor(data$geno)
data$trt = as.factor(data$trt)
y = as.numeric(data$pref) #create outcome vector

X = model.matrix(y ~ geno*trt, data = data, 
                 contrasts.arg = list(geno = 'contr.sum', trt = 'contr.sum'))

# stan data input
trialTime_dat = list(N = length(y),
                    L = dim(X)[2],
                    y, X)

##############
## run STAN ##
##############

setwd("~/de Bivort 2/Enrichment Paper/stan scripts")

# CS light time only
fit <- stan(file = 'CSFlyVac_Gamma.stan', data = trialTime_dat, 
            iter = 100000, chains = 4, control = list(adapt_delta = 0.8, max_treedepth = 10))

# all geno light time only
fit <- stan(file = 'allGeno_GammaModel.stan', data = trialTime_dat, 
            iter = 50000, chains = 4, cores = 4,
            control = list(adapt_delta = 0.8, max_treedepth = 10))


########################
## extract raw values ##
########################

post_dist_raw = as.matrix(fit, pars = c('m', 'v'))
idx = t(post_dist_raw)[!duplicated(t(post_dist_raw))]

## CS ONLY ##
post_dist_raw = t(matrix(idx, 6, 200000))
post_dist_mu = post_dist_raw %>% 
  as_tibble %>%
  select(V1,V2,V3)%>%
  rename(Low = V1, Med = V2, High = V3) %>%
  gather(key = trt, value = mu)

post_dist_var = post_dist_raw %>% 
  as_tibble %>%
  select(V4,V5,V6)%>%
  rename(Low = V4, Med = V5, High = V6) %>%
  gather(key = trt, value = var)


## all geno x trt combo
#create model matrix for mean and var calc later
newdata = expand.grid(geno=levels(data$geno), trt=levels(data$trt))
Xmat = model.matrix(~geno*trt,newdata,
                    contrasts.arg = list(geno = 'contr.sum', trt = 'contr.sum'))

#mean coef
post_dist_raw = as.matrix(fit, pars = c('beta'))
post_dist_raw = post_dist_raw %*% t(Xmat)
post_dist_mu = exp(post_dist_raw) #convert to fitted means for dev calc

post_coef_mu = calcDeviations3Geno(post_dist_mu, c('int','bge[1]','bge[2]','bge[6]','bge[3]',
'bge[4]','bge[5]','bge[7]'))

#var coef
post_dist_raw = as.matrix(fit, pars = c('gamma'))
post_dist_raw = post_dist_raw %*% t(Xmat)
post_dist_var = exp(post_dist_raw) #convert to fitted var for dev calc

post_coef_var = calcDeviations3Geno(post_dist_var, c('int','vge[1]','vge[2]','vge[6]','vge[3]','vge[4]','vge[5]','vge[7]'))

#create fitted CV (sd/mean)
post_dist_cv = sqrt(post_dist_var)/post_dist_mu
post_coef_cv = calcDeviations3Geno(post_dist_cv, c('int','bge[1]','bge[2]','bge[6]','bge[3]','bge[4]','bge[5]','bge[7]'))


#extract simulated response - don't do this if you have more than 2000 iter
sim_trialTime = as.matrix(fit, pars = c("y_rep"))

#save data
setwd("~/de Bivort 2/Enrichment Paper/Stan Output/all Geno Coefs")
save(post_coef_mu, post_coef_var, post_coef_cv, file = 'lightTime_allGenoposteriorCoef.RData')

#################
## DIAGNOSTICS ##
#################
## comparing original dist to simulated dist


#predictor matrix
## CS ONLY
###################
## set the working directory
setwd("~/de Bivort 2/Enrichment Paper/Stan Input")
x = read.csv('predMatrix_CS_FlyVac_all.csv', header = F)
X = data.matrix(x)

color_scheme_set("brightblue")
p1=ppc_dens_overlay(y[as.logical(X[,4])], sim_trialTime[1:50,as.logical(X[,4])])
p2=ppc_dens_overlay(y[as.logical(X[,5])], sim_trialTime[1:50,as.logical(X[,5])])
p3=ppc_dens_overlay(y[as.logical(X[,6])], sim_trialTime[1:50,as.logical(X[,6])])

grid.arrange(p1, p2, p3, nrow=1)

trt = as.factor(c(rep('low', sum(as.logical(X[,2]))),
                  rep('med', sum(as.logical(X[,3]))),
                  rep('high', sum(as.logical(X[,4])))))
trt = relevel(trt,c('low'))
ppc_stat_grouped(y, sim_trialTime, group = trt, stat='sd')
###################

## ALL GENO
library(recipes)
setwd("~/de Bivort 2/Enrichment Paper/Stan Input")
data = read.csv('allGenoTrt_FlyVacTime.csv', header = F)
colnames(data) = c('geno','trt','pref')
data$geno = as.factor(data$geno)
data$trt = as.factor(data$trt)
y = as.numeric(data$pref) #create outcome vector

rec <- recipe(pref ~  geno + trt, data = data) %>%
  step_dummy(geno, one_hot = T) %>% #create dummy vars 
  step_dummy(trt, one_hot = T) %>% #create dummy vars 
  step_interact(terms = ~ matches('geno'):matches('trt')) %>% #introduce interaction
  prep(data = data, retain = TRUE) #perform the above operations to estimate the model
X <- juice(rec, all_predictors(), composition = "matrix")

#density overlay using violin plot
###################
library(RColorBrewer)
library(scales)
#display.brewer.all()
pal = brewer.pal(11,"Spectral")

plot_data = data %>% 
  unite(geno_trt, geno, trt)

breaks = c('0_0','0_1','1_0','1_1','3_0','3_1')
labels = c('CS-Low','CS-Mild','RAL45-Low','RAL45-Mild',
           'RAL535-Low','RAL535-Mild')


q = ggplot(data = plot_data, aes(y = pref, x = geno_trt))

q + 
  geom_violin(aes(x=plot_data$geno_trt, y=sim_trialTime[1,]), color = pal[10], trim=T) +
  geom_jitter(width=0.2, color = pal[10], alpha = 0.3) +
  theme_classic() +
  theme(axis.text.y = element_text(size=14),
        axis.text.x = element_text(angle = 45,size=10, margin = margin(t = 10)),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.length = unit(.25, "cm"),
        axis.line.x = element_blank(),
        axis.line.y = element_line(lineend = 'butt'),
        plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt")) +
  scale_x_discrete(breaks = breaks,
                   labels = labels) +
  scale_y_continuous(breaks = seq(0,120,120/2),
                     labels = c(0,60,120), 
                     limits = c(0,120), expand = c(0,0))
###################


# density overlay with bayesplot
##########################
color_scheme_set("brightblue")
p1=ppc_dens_overlay(y[as.logical(X[,6])], sim_trialTime[1:50,as.logical(X[,6])])
p2=ppc_dens_overlay(y[as.logical(X[,7])], sim_trialTime[1:50,as.logical(X[,7])])
p3=ppc_dens_overlay(y[as.logical(X[,8])], sim_trialTime[1:50,as.logical(X[,8])])
p4=ppc_dens_overlay(y[as.logical(X[,9])], sim_trialTime[1:50,as.logical(X[,9])])
p5=ppc_dens_overlay(y[as.logical(X[,10])], sim_trialTime[1:50,as.logical(X[,10])])
p6=ppc_dens_overlay(y[as.logical(X[,11])], sim_trialTime[1:50,as.logical(X[,11])])
grid.arrange(p1, p2, p3, p4, p5, p6, nrow=2)
##########################

#summary stat comparison
##########################
trt = as.factor(c(rep('g0xlow', sum(as.logical(X[,6]))),
                  rep('g0xmed', sum(as.logical(X[,7]))),
                  rep('g1xlow', sum(as.logical(X[,8]))),
                  rep('g1xmed', sum(as.logical(X[,9]))),
                  rep('g3xlow', sum(as.logical(X[,10]))),
                  rep('g3xmed', sum(as.logical(X[,11])))))

ppc_stat_grouped(y, sim_trialTime, group = trt, stat='sd')
##########################
