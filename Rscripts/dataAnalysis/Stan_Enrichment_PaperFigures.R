################################
##  Enrichment STAN analysis  ##
##  making figures for paper  ##
################################

library(tidyverse)
library(RColorBrewer)
library(scales)
pal = brewer.pal(11,"Spectral")
show_col(brewer_pal(palette = "Spectral")(11))

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
setwd("~/de Bivort 2/Enrichment Paper/Stan Output/CS Means and Vars")
load('lightTime_CSposterior.RData')

#switch factor levels to Low Med High for CS data
##############
post_dist_mu$trt = post_dist_mu$trt %>% recode(a='Low',b1='Med',b2='High')
post_dist_var$trt = post_dist_var$trt %>% recode(v0='Low',v1='Med',v2='High')
##############

#create fitted means and vars from allGeno coefficient data
#turnBias, light Choice, clump, switch
##############
post_dist_mu = post_coef_mu %*% t(Xmat)
post_dist_var = post_coef_var %*% t(Xmat)

#create fitted CV (sd/mean)
post_dist_cv = sqrt(post_dist_var)/post_dist_mu

#create coefficients - turnBias, clump, switch
post_coef_cv = calcDeviations5Geno(post_dist_cv, c('int','bge[1]','bge[2]','bge[3]','bge[4]','bge[10]','bge[5]','bge[6]','bge[7]','bge[8]','bge[9]','bge[11]'))

#create coefficients - light Choice
post_coef_cv = calcDeviations3Geno(post_dist_cv, c('int','bge[1]','bge[2]','bge[6]','bge[3]','bge[4]','bge[5]','bge[7]'))

##############


## Figure 3 ##
#####################
## VIOLIN PLOTS ##
q = ggplot(post_dist_cv, aes(x=factor(trt, levels = c('Low','Med','High')), y=cv))

q + geom_violin(color = pal[2], fill = pal[2]) +
  #stat_summary(fun.data = 'mean_sdl', fun.args = list(mult = 1.96),
  #color='white', size = 0.5) +
  theme_classic() +
  theme(axis.text.y = element_text(size=18),
        axis.text.x = element_text(size=18),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.length = unit(.25, "cm"),
        axis.line.x = element_blank(),
        axis.line.y = element_line(lineend = 'butt'),
        plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt")) +
  scale_x_discrete(breaks = c('Low','Med','High'),
                   labels = c('none', 'mild','intense'),
                   expand = c(0.25, 0.25)) +
  scale_y_continuous(breaks = seq(0,0.6,0.3),
                     labels = c(0, 0.3, 0.6),limits = c(0,0.6), expand = c(0,0))

#####################

## estimate g, e, gxe coefficients ##
## Figure 4 and 5 ##
####################################

#add in G4, and G4xMild coefficients for mean and var, switch sign to get coef for MildEnr
#turnBias, clump, switch
##############
#mean
`bge[10]` = -1*post_coef_mu[,2] + -1*post_coef_mu[,3] + 
                    -1*post_coef_mu[,4] + -1*post_coef_mu[,5]

`bge[11]` = post_coef_mu[,7] + post_coef_mu[,8] + 
                    post_coef_mu[,9] + post_coef_mu[,10]

post_coef_mu = cbind(post_coef_mu, `bge[10]`, `bge[11]`)
post_coef_mu[,c(2:12)] = post_coef_mu[,c(2:12)]/post_coef_mu[,1] #scale by grand mean
post_coef_mu[,c(6:10)] = -1*post_coef_mu[,c(6:10)] #switch sign to get MildEnr coef

coef_mu = post_coef_mu[,-1] %>% as.tibble %>% gather(key = trt, value = mu)
coef_mu$trt = factor(coef_mu$trt, 
                     levels = c('bge[1]','bge[2]','bge[3]','bge[4]','bge[10]',
                                'bge[5]','bge[6]','bge[7]','bge[8]','bge[9]','bge[11]'))


#var
`vge[10]` = -1*post_coef_var[,2] + -1*post_coef_var[,3] + 
  -1*post_coef_var[,4] + -1*post_coef_var[,5]

`vge[11]` = post_coef_var[,7] + post_coef_var[,8] + 
  post_coef_var[,9] + post_coef_var[,10]

post_coef_var = cbind(post_coef_var, `vge[10]`, `vge[11]`)
post_coef_var[,c(2:12)] = post_coef_var[,c(2:12)]/post_coef_var[,1] #scale by grand var
post_coef_var[,c(6:10)] = -1*post_coef_var[,c(6:10)] #switch sign to get MildEnr coef

coef_var = post_coef_var[,-1] %>% as.tibble %>% gather(key = trt, value = var)
coef_var$trt = factor(coef_var$trt, 
                     levels = c('vge[1]','vge[2]','vge[3]','vge[4]','vge[10]',
                                'vge[5]','vge[6]','vge[7]','vge[8]','vge[9]','vge[11]'))
#cv
post_coef_cv[,c(2:12)] = post_coef_cv[,c(2:12)]/post_coef_cv[,1] #scale by grand cv
coef_cv = post_coef_cv[,-1] %>% as.tibble %>% gather(key = trt, value = cv)
coef_cv$trt = factor(coef_cv$trt, 
                     levels = c('bge[1]','bge[2]','bge[3]','bge[4]','bge[10]',
                                'bge[5]','bge[6]','bge[7]','bge[8]','bge[9]','bge[11]'))
##############

#add in G3, and G3xMild coefficients for mean and var, switch sign to get coef for MildEnr
#light Choice
##############
#mean
`bge[6]` = -1*post_coef_mu[,2] + -1*post_coef_mu[,3]

`bge[7]` = post_coef_mu[,5] + post_coef_mu[,6]

post_coef_mu = cbind(post_coef_mu, `bge[6]`, `bge[7]`)
post_coef_mu[,c(2:8)] = post_coef_mu[,c(2:8)]/post_coef_mu[,1] #scale by grand mean
post_coef_mu[,c(4:6)] = -1*post_coef_mu[,c(4:6)] #switch sign to get MildEnr coef

coef_mu = post_coef_mu[,-1] %>% as.tibble %>% gather(key = trt, value = mu)
coef_mu$trt = factor(coef_mu$trt,
                     levels = c('bge[1]','bge[2]','bge[6]','bge[3]',
                                'bge[4]','bge[5]','bge[7]'))

#var
`vge[6]` = -1*post_coef_var[,2] + -1*post_coef_var[,3]

`vge[7]` = post_coef_var[,5] + post_coef_var[,6]

post_coef_var = cbind(post_coef_var, `vge[6]`, `vge[7]`)
post_coef_var[,c(2:8)] = post_coef_var[,c(2:8)]/post_coef_var[,1] #scale by grand var
post_coef_var[,c(4:6)] = -1*post_coef_var[,c(4:6)] #switch sign to get MildEnr coef

coef_var = post_coef_var[,-1] %>% as.tibble %>% gather(key = trt, value = var)
coef_var$trt = factor(coef_var$trt, 
                      levels = c('vge[1]','vge[2]','vge[6]','vge[3]',
                                 'vge[4]','vge[5]','vge[7]'))
#cv
post_coef_cv[,c(2:8)] = post_coef_cv[,c(2:8)]/post_coef_cv[,1] #scale by grand cv
coef_cv = post_coef_cv[,-1] %>% as.tibble %>% gather(key = trt, value = cv)
coef_cv$trt = factor(coef_cv$trt, 
                     levels = c('bge[1]','bge[2]','bge[6]','bge[3]',
                                'bge[4]','bge[5]','bge[7]'))
##############

#numTurns
##############
post_coef_mu[,c(2:12)] = post_coef_mu[,c(2:12)]/post_coef_mu[,1] #scale by grand mean
coef_mu = post_coef_mu[,-1] %>% as.tibble %>% gather(key = trt, value = mu)
coef_mu$trt = factor(coef_mu$trt,
                     levels = c('bge[1]','bge[2]','bge[3]','bge[4]','bge[10]',
                                'bge[5]','bge[6]','bge[7]','bge[8]','bge[9]','bge[11]'))

post_coef_var[,c(2:12)] = post_coef_var[,c(2:12)]/post_coef_var[,1] #scale by grand var
coef_var = post_coef_var[,-1] %>% as.tibble %>% gather(key = trt, value = var)
coef_var$trt = factor(coef_var$trt, 
                      levels = c('vge[1]','vge[2]','vge[3]','vge[4]','vge[10]',
                                 'vge[5]','vge[6]','vge[7]','vge[8]','vge[9]','vge[11]'))

post_coef_cv[,c(2:12)] = post_coef_cv[,c(2:12)]/post_coef_cv[,1] #scale by grand cv
coef_cv = post_coef_cv[,-1] %>% as.tibble %>% gather(key = trt, value = cv)
coef_cv$trt = factor(coef_cv$trt, 
                     levels = c('bge[1]','bge[2]','bge[3]','bge[4]','bge[10]',
                                'bge[5]','bge[6]','bge[7]','bge[8]','bge[9]','bge[11]'))
##############

#lightTime
##############
post_coef_mu[,c(2:8)] = post_coef_mu[,c(2:8)]/post_coef_mu[,1] #scale by grand mean
coef_mu = post_coef_mu[,-1] %>% as.tibble %>% gather(key = trt, value = mu)
coef_mu$trt = factor(coef_mu$trt,
                     levels = c('bge[1]','bge[2]','bge[6]','bge[3]',
                                'bge[4]','bge[5]','bge[7]'))

post_coef_var[,c(2:8)] = post_coef_var[,c(2:8)]/post_coef_var[,1] #scale by grand var
coef_var = post_coef_var[,-1] %>% as.tibble %>% gather(key = trt, value = var)
coef_var$trt = factor(coef_var$trt, 
                      levels = c('vge[1]','vge[2]','vge[6]','vge[3]',
                                 'vge[4]','vge[5]','vge[7]'))
#cv
post_coef_cv[,c(2:8)] = post_coef_cv[,c(2:8)]/post_coef_cv[,1] #scale by grand cv
coef_cv = post_coef_cv[,-1] %>% as.tibble %>% gather(key = trt, value = cv)
coef_cv$trt = factor(coef_cv$trt, 
                     levels = c('bge[1]','bge[2]','bge[6]','bge[3]',
                                'bge[4]','bge[5]','bge[7]'))
##############

# Plotting #
############

#ymaze
breaks = c('bge[1]','bge[2]','bge[3]','bge[4]','bge[10]',
           'bge[5]','bge[6]','bge[7]','bge[8]','bge[9]','bge[11]')

breaks = c('vge[1]','vge[2]','vge[3]','vge[4]','vge[10]',
           'vge[5]','vge[6]','vge[7]','vge[8]','vge[9]','vge[11]')

labels = c('CS', 'RAL45','RAL105','RAL535','RAL796','MildEnr','CS-Mild',
           'RAL45-Mild','RAL105-Mild','RAL535-Mild','RAL796-Mild')

#fly vac
breaks = c('bge[1]','bge[2]','bge[6]','bge[3]',
             'bge[4]','bge[5]','bge[7]')

breaks = c('vge[1]','vge[2]','vge[6]','vge[3]',
           'vge[4]','vge[5]','vge[7]')

labels = c('CS', 'RAL45','RAL535','MildEnr','CS-Mild',
           'RAL45-Mild','RAL535-Mild')

#plotting fxn
q = ggplot(coef_mu, aes(x=trt, y=mu))
q + geom_violin(color = pal[10], fill = pal[10]) +
  #stat_summary(fun.data = 'mean_sdl', fun.args = list(mult = 1.96),
  #color='white', size = 0.5) +
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
  geom_hline(yintercept=0, linetype = 'dotted')+
  scale_x_discrete(breaks = breaks,
                   labels = labels) +
  scale_y_continuous(breaks = seq(-0.3,0.3,0.6/2),
                     labels = c(-0.3,0, 0.3), 
                     limits = c(-0.3,0.3), expand = c(0,0))

############


# Figure 6: summary plot (A,B) #
#Data
############
#############
gEffects_mu = mean(abs(post_coef_mu[,c(2:5,11)]))
eEffects_mu = mean(abs(post_coef_mu[,6]))
geEffects_mu = mean(abs(post_coef_mu[,c(7:10,12)]))

gEffects_cv = mean(abs(post_coef_cv[,c(2:6)]))
eEffects_cv = mean(abs(post_coef_cv[,7]))
geEffects_cv = mean(abs(post_coef_cv[,c(8:12)]))

dd_mu = tibble(x=c('Switch','Switch','Switch'), 
            lab = factor(c('g','e','ge'),levels = c('g','e','ge')),
            y = c(gEffects_mu,eEffects_mu,geEffects_mu))

dd_cv = tibble(x=c('Switch','Switch','Switch'), 
               lab = factor(c('g','e','ge'),levels = c('g','e','ge')),
               y = c(gEffects_cv,eEffects_cv,geEffects_cv))

muEffects = rbind(muEffects, dd_mu)
cvEffects = rbind(cvEffects, dd_cv)
#Plotting
##########
load('allBehavior_GxE_Effects.RData')

muEffects$x = factor(muEffects$x, 
              levels = c('TurnBias','NumTurns','Switch','Clump','LightPref','LightTime'))
cvEffects$x = factor(cvEffects$x, 
              levels = c('TurnBias','NumTurns','Switch','Clump','LightPref','LightTime'))

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#plotting fxn
q = ggplot(muEffects, aes(x=x, y=y))
q + geom_bar(stat='identity',aes(fill=lab), 
  position = position_stack(reverse = TRUE)) +
  coord_flip() +
  theme_classic() +
  theme(axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14, vjust = 0.5),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.length = unit(.25, "cm"),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"),
        legend.text = element_text(size = 14),
        legend.title = element_blank()) +
  scale_fill_manual(breaks = c('g','e','ge'),
                   labels = c('G','E','GxE'),
                   values = cbPalette) +
  scale_y_continuous(breaks = seq(0,0.4,0.1),
                     labels = c(0,0.10,0.20,0.3,0.4),
                     limits = c(0,0.5))
  
#############

#Figure 6: scatterplot (C)#
#Data
########
gEffects_mu = tibble(lab = rep('g',5), y = apply(post_coef_mu[,c(2:6)],2,mean))
eEffects_mu = tibble(lab = 'e', y = mean(post_coef_mu[,7]))
geEffects_mu = tibble(lab = rep('ge',5), y = apply(post_coef_mu[,c(8:12)],2,mean))

gEffects_cv = tibble(lab = rep('g',5), y = apply(post_coef_cv[,c(2:6)],2,mean))
eEffects_cv = tibble(lab = 'e', y = mean(post_coef_cv[,7]))
geEffects_cv = tibble(lab = rep('ge',5), y = apply(post_coef_cv[,c(8:12)],2,mean))

dd_mu = cbind(x=rep('NumTurns',11), 
               rbind(gEffects_mu,eEffects_mu,geEffects_mu))

dd_cv = cbind(x=rep('NumTurns',11), 
              rbind(gEffects_cv,eEffects_cv,geEffects_cv))

muEffects_geno = rbind(muEffects_geno, dd_mu)
cvEffects_geno = rbind(cvEffects_geno, dd_cv)

save(muEffects_geno, cvEffects_geno, file = 'allBehaviorEffects_scatter.RData')
#######
#Plotting
######
load('allBehaviorEffects_scatter.RData')
scatter_data = cbind(muEffects_geno, cvEffects_geno$y)
colnames(scatter_data) = c('x','lab','mu','cv')
scatter_data$lab = factor(scatter_data$lab, levels = c('g','e','ge'))

q = ggplot(scatter_data, aes(x=mu, y=cv, shape = lab))
q + geom_point(aes(color = lab), size = 3.5) +
  theme_classic() +
  theme(axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14, vjust = 0.5),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.length = unit(.25, "cm"),
        axis.line.y = element_line(lineend = 'butt'),
        axis.line.x = element_line(lineend = 'butt'),
        plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"),
        legend.text = element_text(size = 14),
        legend.title = element_blank()) +
  scale_color_manual(breaks = c('g','e','ge'),
                    labels = c('G','E','GxE'),
                   values = cbPalette) +
  scale_shape_manual(breaks = c('g','e','ge'),
                     labels = c('G','E','GxE'),
                     values = c(15, 16, 17))
  
##########
