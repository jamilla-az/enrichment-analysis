##############################
## Enrichment STAN analysis ##
##############################
library(tidyverse)
library(rstan)
library(shinystan)
library(bayesplot)
library(gridExtra)
        
## set the working directory
setwd("~/de Bivort 2/Enrichment Paper/Stan Input")

##############
# functions #
##############
calcDeviations5Geno = function(post_dist, labels){
  
  post_coef = matrix(NA, dim(post_dist)[1], dim(post_dist)[2]+2)
  post_coef[,1] = apply(post_dist, 1, mean) #grand mean
  
  #dev of each geno from grand param value
  post_coef[,2] = apply(post_dist[,c(1,6)], 1, mean) - post_coef[,1]#dev of G0
  post_coef[,3] = apply(post_dist[,c(2,7)], 1, mean) - post_coef[,1]#dev of G1
  post_coef[,4] = apply(post_dist[,c(3,8)], 1, mean) - post_coef[,1]#dev of G2
  post_coef[,5] = apply(post_dist[,c(4,9)], 1, mean) - post_coef[,1]#dev of G3
  post_coef[,6] = apply(post_dist[,c(5,10)], 1, mean) - post_coef[,1]#dev of G4
  
  #dev of mild enrich from grand mean
  post_coef[,7] = apply(post_dist[,c(6:10)], 1, mean) - post_coef[,1]#dev of MildEnr
  
  #dev of GenoxMildEnr
  post_coef[,8] = post_dist[,6] - post_coef[,1] - post_coef[,2]#dev of G0xMildEnr
  post_coef[,9] = post_dist[,7] - post_coef[,1] - post_coef[,3]#dev of G1xMildEnr
  post_coef[,10] = post_dist[,8] - post_coef[,1] - post_coef[,4]#dev of G2xMildEnr
  post_coef[,11] = post_dist[,9] - post_coef[,1] - post_coef[,5]#dev of G3xMildEnr
  post_coef[,12] = post_dist[,10] - post_coef[,1] - post_coef[,6]#dev of G4xMildEn
  
  colnames(post_coef) = labels
  return(post_coef)
}
###############
## load data ##
###############

## num turns
#outcome vector -CS
y = read.csv('numTurns_CS_YMaze.csv', header = F)
y = as.numeric(y$V1)

#predictor matrix -CS
x = read.csv('predMatrix_CS_FlyVac_all.csv', header = F)
X = data.matrix(x[,2:4])

## num turns for all genotypes and trt combos
setwd("~/de Bivort 2/Enrichment Paper/Stan Input")
data = read.csv('allGenoTrt_YMazeNumTurn.csv', header = F)
colnames(data) = c('geno','trt','pref')
data$geno = as.factor(data$geno)
data$trt = as.factor(data$trt)
y = as.numeric(data$pref) #create outcome vector

X = model.matrix(y ~ geno*trt, data = data, 
                 contrasts.arg = list(geno = 'contr.sum', trt = 'contr.sum'))

# stan data input
numTurn_dat = list(N = length(y),
                   L = dim(X)[2],
                   y, X)

#run STAN
setwd("~/de Bivort 2/Enrichment Paper/stan scripts")

# CS num Turns only
fit <- stan(file = 'CSymazeNumTurn_NegBin.stan', data = numTurn_dat, 
            iter = 100000, chains = 4, control = list(max_treedepth = 15))

# all Geno num Turns only
fit <- stan(file = 'allGeno_NegBinModel.stan', data = numTurn_dat, 
            iter = 50000, chains = 4, cores = 4, control = list(max_treedepth = 15))

#launch shiny stan
my_sso <- launch_shinystan(fit)

########################
## extract raw values ##
########################

## CS ONLY ##
post_dist_raw = as.matrix(fit, pars = c('m', 'v'))
idx = t(post_dist_raw)[!duplicated(t(post_dist_raw))]

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

#save data
setwd("~/de Bivort 2/Enrichment Paper/Stan Output")
save(post_dist_mu, post_dist_var, file = 'numTurns_CSposterior.RData')


## all geno x trt combo
#create model matrix for mean and var calc later
newdata = expand.grid(geno=levels(data$geno), trt=levels(data$trt))
Xmat = model.matrix(~geno*trt,newdata,
                    contrasts.arg = list(geno = 'contr.sum', trt = 'contr.sum'))

#mean coef
post_dist_raw = as.matrix(fit, pars = c('beta'))
post_dist_raw = post_dist_raw %*% t(Xmat)
post_dist_mu = exp(post_dist_raw) #convert to fitted means for dev calc

post_coef_mu = calcDeviations5Geno(post_dist_mu, c('int','bge[1]','bge[2]','bge[3]','bge[4]','bge[10]','bge[5]','bge[6]','bge[7]','bge[8]','bge[9]','bge[11]'))

#var coef
post_dist_raw = as.matrix(fit, pars = c('gamma'))
post_dist_raw = post_dist_raw %*% t(Xmat)
post_dist_phi = exp(post_dist_raw) #convert to fitted dispersion for dev calc

post_dist_var = post_dist_mu + (((post_dist_mu)^2)/post_dist_phi) #convert to var

post_coef_var = calcDeviations5Geno(post_dist_var, c('int','vge[1]','vge[2]','vge[3]','vge[4]','vge[10]','vge[5]','vge[6]','vge[7]','vge[8]','vge[9]','vge[11]') )

#create fitted CV (sd/mean)
post_dist_cv = sqrt(post_dist_var)/post_dist_mu

post_coef_cv = calcDeviations5Geno(post_dist_cv, c('int','bge[1]','bge[2]','bge[3]','bge[4]','bge[10]','bge[5]','bge[6]','bge[7]','bge[8]','bge[9]','bge[11]'))

#extract simulated response - don't do this if more than 2000 iterations
sim_numTurns = as.matrix(fit, pars = c("y_rep"))

#save data
setwd("~/de Bivort 2/Enrichment Paper/Stan Output/all Geno Coefs")
save(post_coef_mu, post_coef_var, post_coef_cv, file = 'numTurns_allGenoposteriorCoef.RData')


#################
## DIAGNOSTICS ##
#################
## comparing original dist to simulated dist

##CS
##################
## set the working directory
setwd("~/de Bivort 2/Enrichment Paper/Stan Input")
x = read.csv('predMatrix_CS_YMaze_all.csv', header = F)
X = data.matrix(x)

color_scheme_set("brightblue")
p1=ppc_dens_overlay(y[as.logical(X[,2])], sim_numTurns[1:50,as.logical(X[,2])])
p2=ppc_dens_overlay(y[as.logical(X[,3])], sim_numTurns[1:50,as.logical(X[,3])])
p3=ppc_dens_overlay(y[as.logical(X[,4])], sim_numTurns[1:50,as.logical(X[,4])])

grid.arrange(p1, p2, p3, nrow=1)

trt = as.factor(c(rep('low', sum(as.logical(X[,2]))),
                  rep('med', sum(as.logical(X[,3]))),
                  rep('high', sum(as.logical(X[,4])))))
trt = relevel(trt,c('low'))
ppc_stat_grouped(y, sim_numTurns, group = trt, stat='sd')
##################

## all geno
setwd("~/de Bivort 2/Enrichment Paper/Stan Input")
library(recipes)
data = read.csv('allGenoTrt_YMazeNumTurn.csv', header = F)
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

breaks = c('0_0','0_1','1_0','1_1','2_0','2_1','3_0','3_1','4_0','4_1')
labels = c('CS-Low','CS-Mild','RAL45-Low','RAL45-Mild','RAL105-Low','RAL105-Mild',
           'RAL535-Low','RAL535-Mild','RAL796-Low','RAL796-Mild')

q = ggplot(data = plot_data, aes(y = pref, x = geno_trt))

q + 
  geom_violin(aes(x=plot_data$geno_trt, y=sim_numTurns[1,]), color = pal[3], trim=T) +
  geom_jitter(width=0.2, color = pal[3], alpha = 0.3) +
  theme_classic() +
  theme(axis.text.y = element_text(size=14),
        axis.text.x = element_text(angle = 45,size=10,
                                   margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.length = unit(.25, "cm"),
        axis.line.x = element_blank(),
        axis.line.y = element_line(lineend = 'butt'),
        plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt")) +
  scale_x_discrete(breaks = breaks,
                   labels = labels) +
  scale_y_continuous(breaks = seq(0,2500,2500/2),
                     labels = c(0,1250,2500), 
                     limits = c(0,2500), expand = c(0,0))
###################



# density overlay with bayes plot
###########################
color_scheme_set("brightblue")
p1=ppc_dens_overlay(y[as.logical(X[,8])], sim_numTurns[1:50,as.logical(X[,8])])
p2=ppc_dens_overlay(y[as.logical(X[,9])], sim_numTurns[1:50,as.logical(X[,9])])
p3=ppc_dens_overlay(y[as.logical(X[,10])], sim_numTurns[1:50,as.logical(X[,10])])
p4=ppc_dens_overlay(y[as.logical(X[,11])], sim_numTurns[1:50,as.logical(X[,11])])
p5=ppc_dens_overlay(y[as.logical(X[,12])], sim_numTurns[1:50,as.logical(X[,12])])
p6=ppc_dens_overlay(y[as.logical(X[,13])], sim_numTurns[1:50,as.logical(X[,13])])
p7=ppc_dens_overlay(y[as.logical(X[,14])], sim_numTurns[1:50,as.logical(X[,14])])
p8=ppc_dens_overlay(y[as.logical(X[,15])], sim_numTurns[1:50,as.logical(X[,15])])
p9=ppc_dens_overlay(y[as.logical(X[,16])], sim_numTurns[1:50,as.logical(X[,16])])
p10=ppc_dens_overlay(y[as.logical(X[,17])], sim_numTurns[1:50,as.logical(X[,17])])
grid.arrange(p1, p2, p3, p4, p5, p6, p7, 
             p8, p9,p10, nrow=3)
###########################

#summary stat comparison
###########################
trt = as.factor(c(rep('g0xlow', sum(as.logical(X[,8]))),
                  rep('g0xmed', sum(as.logical(X[,9]))),
                  rep('g1xlow', sum(as.logical(X[,10]))),
                  rep('g1xmed', sum(as.logical(X[,11]))),
                  rep('g2xlow', sum(as.logical(X[,12]))),
                  rep('g2xmed', sum(as.logical(X[,13]))),
                  rep('g3xlow', sum(as.logical(X[,14]))),
                  rep('g3xmed', sum(as.logical(X[,15]))),
                  rep('g4xlow', sum(as.logical(X[,16]))),
                  rep('g4xmed', sum(as.logical(X[,17])))))

ppc_stat_grouped(y, sim_numTurns, group = trt, stat='mean')
###########################

