##############################
## Enrichment STAN analysis ##
##############################
library(tidyverse)
library(rstan)
library(shinystan)
library(bayesplot)
library(gridExtra)
# all genotypes #
# turn bias, clump, switch, lightPref #
###############
## load data ##
###############
##all genotypes and trt combos ##
setwd("~/de Bivort 2/Enrichment Paper/Stan Input")
data = read.csv('allGenoTrt_YMazePref.csv', header = F)
colnames(data) = c('geno','trt','pref')
data$geno = as.factor(data$geno)
data$trt = as.factor(data$trt)
y = as.numeric(data$pref) #create outcome vector

mat = model.matrix(y ~ geno*trt, data = data, 
                contrasts.arg = list(geno = 'contr.sum', trt = 'contr.sum'))
X = mat[,-1]


# stan data input
pheno_dat = list(N = length(y),
                    L = dim(X)[2],
                    y, X)
##############
## run STAN ##
##############

setwd("~/de Bivort 2/Enrichment Paper/stan scripts")
library(rstan)

#all Geno only
fit <- stan(file = 'allGeno_NormalModel.stan', data = pheno_dat, 
            iter = 100000, chains = 4, cores = 4,
            control = list(adapt_delta = 0.9, max_treedepth = 10))

#launch shiny stan
library(shinystan)
my_sso <- launch_shinystan(fit)

########################
## extract raw values ##
########################

#create model matrix for mean and var calc later
newdata = expand.grid(geno=levels(data$geno), trt=levels(data$trt))
Xmat = model.matrix(~geno*trt,newdata,
                    contrasts.arg = list(geno = 'contr.sum', trt = 'contr.sum'))
#mean coef
post_dist_raw = as.matrix(fit, pars = c('a','bge'))
post_coef_mu = post_dist_raw

#var coef
post_dist_raw = as.matrix(fit, pars = c('v0','vge'))
post_coef_var = post_dist_raw

#extract simulated response
sim_turnBias = as.matrix(fit, pars = c("y_rep"))

#save data
setwd("~/de Bivort 2/Enrichment Paper/Stan Output/all Geno Coefs")
save(post_coef_mu, post_coef_var, Xmat, file = 'turnBias_allGenoposteriorCoef.RData')

#################
## DIAGNOSTICS ##
#################

## comparing original dist to simulated dist
library(bayesplot)
library(gridExtra)

## ALL GENO X TRT COMBO
#predictor matrix
library(recipes)

rec <- recipe(pref ~ geno + trt, data = data) %>%
  step_dummy(geno, one_hot = T) %>%
  step_dummy(trt, one_hot = T) %>%
  step_interact(terms = ~ matches('geno'):matches('trt')) %>% #introduce interaction
  prep(data = data, retain = TRUE)
allX <- juice(rec, all_predictors(), composition = "matrix") #create predictor matrix

X = allX

#density overlay using bayesplot
###################
color_scheme_set("brightblue")
p1=ppc_dens_overlay(y[as.logical(X[,8])], sim_turnBias[1:50,as.logical(X[,8])])
p2=ppc_dens_overlay(y[as.logical(X[,9])], sim_turnBias[1:50,as.logical(X[,9])])
p3=ppc_dens_overlay(y[as.logical(X[,10])], sim_turnBias[1:50,as.logical(X[,10])])
p4=ppc_dens_overlay(y[as.logical(X[,11])], sim_turnBias[1:50,as.logical(X[,11])])
p5=ppc_dens_overlay(y[as.logical(X[,12])], sim_turnBias[1:50,as.logical(X[,12])])
p6=ppc_dens_overlay(y[as.logical(X[,13])], sim_turnBias[1:50,as.logical(X[,13])])
p7=ppc_dens_overlay(y[as.logical(X[,14])], sim_turnBias[1:50,as.logical(X[,14])])
p8=ppc_dens_overlay(y[as.logical(X[,15])], sim_turnBias[1:50,as.logical(X[,15])])
p9=ppc_dens_overlay(y[as.logical(X[,16])], sim_turnBias[1:50,as.logical(X[,16])])
p10=ppc_dens_overlay(y[as.logical(X[,17])], sim_turnBias[1:50,as.logical(X[,17])])
grid.arrange(p1, p2, p3, p4, p5, p6, p7, 
             p8, p9,p10, nrow=3)
###################

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

breaks = c('0_0','0_1','1_0','1_1','3_0','3_1')
labels = c('CS-Low','CS-Mild','RAL45-Low','RAL45-Mild',
           'RAL535-Low','RAL535-Mild')


q = ggplot(data = plot_data, aes(y = pref, x = geno_trt))

q + 
  geom_violin(aes(x=plot_data$geno_trt, y=sim_turnBias[1,]), color = pal[11], trim=T) +
  geom_jitter(width=0.2, color = pal[11], alpha = 0.3) +
  theme_classic() +
  theme(axis.text.y = element_text(size=14),
        axis.text.x = element_text(angle = 45,size=10),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.length = unit(.25, "cm"),
        axis.line.x = element_blank(),
        axis.line.y = element_line(lineend = 'butt'),
        plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt")) +
  scale_x_discrete(breaks = breaks,
                   labels = labels) +
  scale_y_continuous(breaks = seq(0,1.5,1.5/2),
                     labels = c(0,0.75,1.5), 
                     limits = c(0,1.5), expand = c(0,0))
###################

#summary stat comparison
###################
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

ppc_stat_grouped(y, sim_turnBias, group = trt, stat='sd')
###################
