#######################################################################
# Heaping bias:
# Simulation study to show that imputing the potentially-rounded y's
# decreases bias in the posterior mean estimated responses.
#######################################################################

#######################################################################
# Workspace prep.
#######################################################################

rm(list=ls())

library(FNN)
library(scales)
library(tidyverse)

library(gridExtra)
library(ggthemes)
source('./code/helper-functions/ggtheme-publication.R')

# Define the true underlying function.
#f_true = function(x) log(1 + exp(x))
f_true = function(x) 10/{1 + exp(-x)} - x
curve(f_true, from=-5,to=5)

# Number of Monte Carlo simulations, and sample size.
NMC = 1000
N = 500

# Rounding settings.
cc = 1
n_digits = -log10(cc)

# Read in x's for non-uniform x.
mulago_ga = rescale(unlist(read.csv('./data/mulago_scrubs2.csv') %>% select(gest_age)), to=c(-5,5))
hist(mulago_ga)

x_mulag = rescale(rbeta(N, 10, 2), to=c(-5,5))
hist(x_mulag)

ga_scaled = rescale(mulago_ga, to=c(-5,5)) #Recale between -5, 5.
hist(ga_scaled)

#######################################################################
### Function to create table for mse analysis over range of 
### scales for sig2 vs rounding window (3x to 1/3x).
#######################################################################
#   -  cc = rounding setting (to nearest digit).
#   -  alpha = vector of scales for sig2 vs rounding size; sigma = cc*alpha.
#   -  unif_x = T or F.  If F, mirrors mulago gest age distribution.
#   -  k = k-value for KNN regression.  k=3 is less smooth; increase values to increase smoothing.
#   -  N = sample size.
#   -  NMC = number of monte carlo draws.

mseRatio = function(cc = 1, alpha = seq(.3,3,by=.3), unif_x = T, k=3, N=500, NMC=1000){
   
   ### Set up dataframe to hold results.
   # cc = rounding setting.
   # alpha = rounding-to-sigma ratio.  (cc*alpha = sigma). 
   # sigma = calculated sigma value.
   # rho = actual variance inflation:  var(ytilde) / sigma2
   # mse, mse_tilde:  the mse values fitting model without, with rounded y's.
   # mse_ratio: mse_tilde / mse
   
   df = cbind.data.frame(cc, k, N, NMC, unif_x, alpha)
   df$sigma = df$cc * df$alpha
   df$rho = 0 
   df$mse = 0 
   df$mse_tilde = 0
   df$mse_ratio = 0
   
   # Loop through alpha scenarios, fill in table.
   for(i in 1:nrow(df)){
      
      print(paste0('Scenario ', i, ' of ',nrow(df)))
      
      # For each scenario:
      # Draw x's, calculate fx.
      
      if(unif_x){
         x_grid = runif(N,-5,5)
      } else{
         #x_grid = rnorm(N, 1.1, 2)
         x_grid = rescale(rbeta(N, 10, 2), to=c(-5,5))
      }

      f_grid = f_true(x_grid)
      
      # Temp sigma value for scenario.
      sigma = df$sigma[i]
      
      # Track sums of rho, mse, mse_tilde; averaging these over NMC.
      rho = 0
      mse = 0
      mse_tilde = 0
      mse_ratio = 0
      
      # Loop through Monte Carlo draws.
      for(j in 1:NMC){
         
         # Calculate y, y_tilde.
         y = rnorm(N, f_grid, sigma)
         y_tilde = round(y, n_digits)
         
         # Fit models with y and y_tilde.
         fit = knn.reg(matrix(x_grid, ncol=1), train=matrix(x_grid, ncol=1), y=y, k=10)
         fhat = fit$pred
         fit_tilde = knn.reg(matrix(x_grid, ncol=1), train=matrix(x_grid, ncol=1), y=y_tilde, k=10)
         fhat_tilde = fit_tilde$pred
         
         #plot(x_grid, fhat)
         #plot(x_grid, fhat_tilde)
         
         # Calculate mse.
         rho = rho + var(y_tilde)/sigma
         mse = mse + mean((fhat - f_grid)^2)
         mse_tilde = mse_tilde + mean((fhat_tilde - f_grid)^2)
         mse_ratio = mse_ratio + (mse_tilde / mse)
      }
      
      # Calculate average rho, mse, across Monte Carlo draws.
      df$rho[i] = rho / NMC
      df$mse[i] = mse / NMC
      df$mse_tilde[i] = mse_tilde / NMC
      df$mse_ratio[i] = mse_ratio / NMC
   }
   
   return(df)
}



#######################################################################
# Run several scenarios.
#######################################################################

# Uniform x's.
df1 = mseRatio(cc=1, unif_x=T, k=10, NMC=1000)

# Non-uniform x's.
df2 = mseRatio(cc=1, unif_x=F, k=10, NMC=1000)

# Uniform x's with smoothing.
df3 = mseRatio(cc=1, unif_x=T, k=30, NMC=1000)

# Non-uniform x's with smoothing.
df4 = mseRatio(cc=1, unif_x=F, k=30, NMC=1000)

#######################################################################
# Plot results.
#######################################################################

df = rbind.data.frame(df1, df2, df3, df4)
df$scenario = rep(1:4, each=nrow(df1))
df$smoothed = ifelse(df$k==10,"Less smooth (k=10)", "Smooth (k=30)")

ggplot(df, aes(x=alpha, y=mse_ratio, linetype=smoothed, color=unif_x)) + 
   geom_line(size=.8) + 
   scale_colour_manual(name="", labels=c('Non-uniform x        ', 'Uniform x        '), 
                       values=c("#000000","#D55E00")) +
   scale_linetype_discrete(name="") + 
   labs(x='Scale of sigma to rounding width', y='MSE inflation ratio') +
   coord_cartesian(ylim=c(1,2.1)) + 
   theme_Publication(legend_position='bottom', legend_direction='vertical') +
   theme(legend.key.size= unit(1, "cm"))

ggsave('./output/01-sim-study-mse-inflation/mse-inflation-ratio.pdf', width=6, height=5)


#######################################################################
# Plotting to show what's going on behind the scenes - illustrate scenarios.
#######################################################################

# This is for a single monte carlo draw of MSE, to visualize what's going on.
# Set up grid of scenarios.
# Note:  This is in long format for ggplot, so y and y_tilde are each their own secnarios here.

scenarios = expand.grid(round_y = c(T,F), unif_x=c(T,F), cc=1, alpha=c(.3,3), k=c(3,10, 30))
scenarios$sigma = scenarios$cc * scenarios$alpha
scenarios$scenario = 1:nrow(scenarios)

df = data.frame(matrix(0, nrow=0, ncol=13))
colnames(df) = c(colnames(scenarios), 'obs','x_grid', 'f_grid', 'y', 'fhat', 'sq_err')

# Loop through scenarios.
myseed = 456

for(i in 1:nrow(scenarios)){
   
   temp = suppressWarnings(cbind.data.frame(scenarios[i,],'obs'= 1:N))
   
   # Generate data.
   if(scenarios$unif_x[i]==T){
     set.seed(myseed)
      temp$x_grid = runif(N,-5,5)
   } else{
     set.seed(myseed)
      temp$x_grid = rescale(rbeta(N, 10, 2), to=c(-5,5))
   }
   
   temp$f_grid = f_true(temp$x_grid)
   
   # Generate y's, with or without rounding, depending on scenario.
   if(scenarios$round_y[i]==F){
     set.seed(myseed)
      temp$y = rnorm(N, temp$f_grid, scenarios$sigma[i])
   } else{
     set.seed(myseed)
      temp$y = round(rnorm(N, temp$f_grid, scenarios$sigma[i]), n_digits)
   }
   
   # Fit model.
   fit = knn.reg(matrix(temp$x_grid, ncol=1), 
                 train=matrix(temp$x_grid, ncol=1), y=temp$y, k=scenarios$k[i])
   temp$fhat = fit$pred
   
   # Squared error.
   temp$sq_err = (temp$fhat - temp$f_grid)^2	
   
   # Add temp data for scenario to results df.
   df = rbind.data.frame(df, temp)
}

# Plot results.
df$facet1 = paste0('scale ', df$alpha)

panel_a = ggplot(df %>% filter(unif_x==T, k==10), aes(x=x_grid, y=y)) + 
   geom_point(alpha=.15, colour='grey50') +
   geom_line(aes(y=f_grid), size=.6, colour='black') + 
   geom_line(aes(y=fhat), size=.6, colour="#D55E00") +
   facet_grid(facet1~ifelse(round_y==1, "y rounded", "y not rounded"), scales='free') + 
   labs(x='x', y='y', subtitle='A) Less smooth fit (k=10).') +
   theme_Publication()

panel_b = ggplot(df %>% filter(unif_x==T, k==30), aes(x=x_grid, y=y)) + 
   geom_point(alpha=.15, colour='grey50') +
   geom_line(aes(y=f_grid), size=.6, colour='black') + 
   geom_line(aes(y=fhat), size=.6, colour="#D55E00") +
   facet_grid(facet1~ifelse(round_y==1, "y rounded", "y not rounded"), scales='free') + 
   labs(x='x', y='y', subtitle='B) Smooth fit (k=30).') +
   theme_Publication()

plt = grid.arrange(panel_a, panel_b, ncol=2)
ggsave('./output/01-sim-study-mse-inflation/mse-sim-illustration.pdf', plot(plt), height=8, width=12)


