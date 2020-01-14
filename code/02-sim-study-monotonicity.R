#######################################################################
# Heaping bias:
# Simulation study to show that imputing the potentially-rounded y's
# decreases bias in the posterior mean estimated responses.
#######################################################################


###########################################################################
# 1. Workspace prep and data setup.
###########################################################################

rm(list=ls())

library(xtable)
library(devtools)
library(data.table)
library(Rcpp)
library(RcppArmadillo)
library(tidyverse)
library(tsbart)
library(gridExtra)
library(sampling)
library(gtools)
library(rpart)       # For CART.
library(rpart.plot)  # For plotting CART.
library(xtable)      # For Table 1 output to txt.

install_github('jestarling/tsbart')
source('./code/helper-functions/ggtheme-publication.R')

nsim=1000; nburn=100; ntree=200

tgrid = seq(1,10,by=1)


#======================================================================
# Create three functions.
#======================================================================

# Shifted sigmoid; looks kind of like example data.
fx1 = function(t,x1,x2){
   2 + .5*x2 + 1 / (1 + exp(-x1 * (t-5)))
}

# Atan.  Could just drop this and do two. (1/3), all three are ok too.
fx2 = function(t,x1,x2){
   2 + atan(t) + .25*x1*(t/5) - .5*x2
}

# Third function, linear in t.
fx3 = function(t,x1,x2){
   2 + t/(8*x1) + x2/2
}

#======================================================================
# Set up functions and visualize data.
#======================================================================

# Plot 25 trajectories.
n=25

# Toy dataset to visualize functions. (All x combos at every t.)
dat = cbind.data.frame(
   't' = rep(tgrid, each=n),
   'x1' = rep(runif(n, .5,1.5), times=length(tgrid)),
   'x2' = rep(runif(n, 0,1), times=length(tgrid)))

dat$Sigmoid = fx1(dat$t, dat$x1, dat$x2)
dat$Arctan = fx2(dat$t, dat$x1, dat$x2)
dat$Linear = fx3(dat$t, dat$x1, dat$x2)


# Plot panel of scenarios.
dat = dat %>% gather(key='scenario',value='fx',-c('t','x1','x2'))

plt_scenarios = ggplot(dat, aes(x=t,colour=factor(x1*x2))) + 
   geom_line(aes(y=fx), alpha=.5) + 
   scale_color_manual(name='', values=rep('black',nrow(dat))) +
   guides(colour=F) + 
   facet_wrap(~scenario,ncol=3) +
   theme_Publication()
   
   ggsave('./output/figure-04.pdf', plot(plt_scenarios), width=9, height=4, dpi=800)
   ggsave('./output/figure-04.tiff', plot(plt_scenarios), width=9, height=4, dpi=800)


#======================================================================
# Repeat the following for 100 replicates of each scenario.
#======================================================================

oos_mse_list = list()
nreps = 50

for(k in 1:nreps){
   
   print(paste0('Replicate ', k, ' of 50.'))
   
   #======================================================================
   # Generate data for analysis.
   #======================================================================
   
   # Number of observations.
   n = 1000 
   tgrid = seq(1,10,by=1)
   
   # Generate dataset.  
   df = cbind.data.frame(
      'id' = 1:n,
      't' = sample(tgrid, size=n, replace=T),
      'x1' = runif(n, .5, 1.5), 
      'x2' = runif(n, 0, 1))
   
   # Set up test/train samples.
   train_idx = sample(1:n, size=.8*n, replace=F)
   df$testtrain = 'test'; df$testtrain[train_idx]='train'
   
   # Calculate fx values for each scenario.
   df$Sigmoid = fx1(df$t, df$x1, df$x2)
   df$Arctan = fx2(df$t, df$x1, df$x2)
   df$Linear = fx3(df$t, df$x1, df$x2)
   
   # Stack dataset for each scenario.
   df = df %>% gather(key='scenario',value='fx',-c('id','t','x1','x2','testtrain'))
   
   # Draw random erros and set up y.
   df$epsilon = rnorm(nrow(df), 0, 1)
   df$y = df$fx + df$epsilon
   
   
   #######################################################################
   # Fit tsBART with and without monotonicity for each scenario.
   #######################################################################
   
   # Separate datasets into test/train for each scenario.
   train_atan = df %>% filter(scenario=='Arctan', testtrain=='train')
   test_atan = df %>% filter(scenario=='Arctan', testtrain=='test')
   
   train_lin = df %>% filter(scenario=='Linear', testtrain=='train')
   test_lin = df %>% filter(scenario=='Linear', testtrain=='test')
   
   train_sig = df %>% filter(scenario=='Sigmoid', testtrain=='train')
   test_sig = df %>% filter(scenario=='Sigmoid', testtrain=='test')
   
   # Fit tsBART and BART models for each scenario.
   fit_atan = tsbart(y=train_atan$y, tgt=train_atan$t, x=train_atan[,c('x1','x2')], 
                     tpred=test_atan$t,xpred=test_atan[,c('x1','x2')], 
                     nsim=nsim, nburn=nburn, ntree=ntree, monotone="no")
   fit_atan_mono = tsbart(y=train_atan$y, tgt=train_atan$t, x=train_atan[,c('x1','x2')], 
                          tpred=test_atan$t,xpred=test_atan[,c('x1','x2')], 
                          nsim=nsim, nburn=nburn, ntree=ntree, monotone="incr")
   
   fit_lin = tsbart(y=train_lin$y, tgt=train_lin$t, x=train_lin[,c('x1','x2')], 
                    tpred=test_lin$t,xpred=test_lin[,c('x1','x2')], 
                    nsim=nsim, nburn=nburn, ntree=ntree, monotone="no")
   fit_lin_mono = tsbart(y=train_lin$y, tgt=train_lin$t, x=train_lin[,c('x1','x2')], 
                         tpred=test_lin$t,xpred=test_lin[,c('x1','x2')], 
                         nsim=nsim, nburn=nburn, ntree=ntree, monotone="incr")
   
   fit_sig = tsbart(y=train_sig$y, tgt=train_sig$t, x=train_sig[,c('x1','x2')], 
                    tpred=test_sig$t,xpred=test_sig[,c('x1','x2')], 
                    nsim=nsim, nburn=nburn, ntree=ntree, monotone="no")
   fit_sig_mono = tsbart(y=train_sig$y, tgt=train_sig$t, x=train_sig[,c('x1','x2')], 
                         tpred=test_sig$t,xpred=test_sig[,c('x1','x2')], 
                         nsim=nsim, nburn=nburn, ntree=ntree, monotone="incr")
   
   #######################################################################
   # Calculate posterior means for each observation.
   #######################################################################
   
   # Acrtan.
   train_atan$pmean = apply(fit_atan$mcmcdraws, 2, mean)
   train_atan$pmean_mono = apply(fit_atan_mono$mcmcdraws, 2, mean)
   
   test_atan$pmean = apply(fit_atan$mcmcdraws_oos, 2, mean)
   test_atan$pmean_mono = apply(fit_atan_mono$mcmcdraws_oos, 2, mean)
   
   # Linear.
   train_lin$pmean = apply(fit_lin$mcmcdraws, 2, mean)
   train_lin$pmean_mono = apply(fit_lin_mono$mcmcdraws, 2, mean)
   
   test_lin$pmean = apply(fit_lin$mcmcdraws_oos, 2, mean)
   test_lin$pmean_mono = apply(fit_lin_mono$mcmcdraws_oos, 2, mean)
   
   # Sigmoid.
   train_sig$pmean = apply(fit_sig$mcmcdraws, 2, mean)
   train_sig$pmean_mono = apply(fit_sig_mono$mcmcdraws, 2, mean)
   
   test_sig$pmean = apply(fit_sig$mcmcdraws_oos, 2, mean)
   test_sig$pmean_mono = apply(fit_sig_mono$mcmcdraws_oos, 2, mean)
   
   #######################################################################
   # Calculate overall MSE, and MSE at each t.
   #######################################################################
   
   # MSE utility function.
   mse = function(y,yhat){
      mean((y-yhat)^2)
   }
   
   overall_oos_mse = cbind.data.frame(
      'scenario' = c('Arctan', 'Linear', 'Sigmoid'),
      'mse' = c(mse(test_atan$y, test_atan$pmean),
                       mse(test_lin$y, test_lin$pmean),
                       mse(test_sig$y, test_sig$pmean)),
      'mse_mono' = c(mse(test_atan$y, test_atan$pmean_mono),
                            mse(test_lin$y, test_lin$pmean_mono),
                            mse(test_sig$y, test_sig$pmean_mono))
   )
   
   oos_mse_list[[k]] = overall_oos_mse
   
}

oos_mse_replicates = do.call(rbind.data.frame, oos_mse_list)
oos_mse_replicates$rep = rep(1:nreps, each=3)
write.csv(oos_mse_replicates, './output/02-sim-study-monotonicity/oos-mse-results-allreplicates.csv')

# Average over replicates.
oos_mse = oos_mse_replicates %>% 
   dplyr::group_by(scenario) %>%
   summarize('mse' = mean(mse), 'mse_mono' = mean(mse_mono))

oos_mse = as.data.frame(oos_mse)
oos_mse[,c(2:3)] = round(oos_mse[,c(2:3)],3)
oos_mse$pct_reduction = (oos_mse$mse - oos_mse$mse_mono) / oos_mse$mse * 100
oos_mse[,4] = round(oos_mse[,4],2)

# Write results to csv file.   
write.csv(oos_mse, './output/table-02.csv')

# Write results to table.
tabl = as.data.frame(oos_mse)
colnames(tabl) = c('Scenario', 'MSE', 'MSE Monotone', 'Percent MSE Reduction')
tabl[,2:4] = round(tabl[,2:4],4)
tabl = xtable(tabl)
print(tabl, file="./output/table-02.txt")

# Save workspace.
#save.image('./output/sim-study-monotonicity-wkspc.RData')


