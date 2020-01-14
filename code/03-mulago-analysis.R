#####################################################################################
# Mulago analysis using tsBART with monotonicity and treatment for
# heaping bias due to rounding of some birth weights ot nearest 100g.
#####################################################################################

###########################################################################
# 1. Workspace prep and data setup.
###########################################################################

rm(list=ls())

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
library(ggrepel)

source('./code/helper-functions/ggtheme-publication.R')
source('./code/helper-functions/subgroupPosteriorInfo.R')
source('./code/helper-functions/testtrain-split.R')

install_github('jestarling/tsbart')

#======================================================================
# Read dataset and rename variables.
#     Variable renaming:
#        pec: newpet
#        gestational age: gest_age
#        birthweight: bw   
#        baby's sex: sex
#     We exclude fever since it's only available for the wave 3 data.  Previous analysis shows
#     it is not a significant predictor of birth weight.
#======================================================================s

# Read dataset and rename variables.
df = fread('data/mulago_data.csv', stringsAsFactors=T) %>%
   rename(pec = newpet, bw_kg = bw, infantsex = sex, 
          momage = mat_age, momjob = mat_jobtype, hiv=HIV, stillbirth=sb, protein=protein_level) %>%
   # Add a patient ID.
   mutate(id = 1:n(), bw_g = bw_kg * 1000, gest_age28 = gest_age-27) %>%
   # Convert infantsex into easy indicator variable.
   mutate(femaleinfant = ifelse(infantsex=="female",1,0)) %>%
   # Indicator for observations which are potentially rounded to nearest 100 g.
   mutate(maybe_rounded = ifelse(bw_g %% 100 == 0, 1, 0)) %>%
   # Select variables to keep.
   select(id, bw_kg, bw_g, maybe_rounded, pec, gest_age, gest_age28, femaleinfant, momage, momjob, ganda, hiv, bpsys_high,
          bpdias_high, hiv, protein, primip, testtrain) %>%
   # Add decade factor for momage.
   mutate(momage_decade = ifelse(momage<=20,"(13,20]",
                                 ifelse(momage<30,"(20,30]",
                                        ifelse(momage<40,"(30,40]","(40,46]"))))
# Convert ganda and hiv to indicators.
df$ganda = ifelse(df$ganda=="yes",1,0)
df$hiv = ifelse(df$hiv=="yes",1,0)

# Keep complete cases.
df = na.omit(df)

# Scale and center covariates which need it (momage, bpsys_high, bpdias_high).
scale_info = rbind.data.frame(cbind.data.frame('mean'=mean(df$momage), 'sd'=sd(df$momage)),
                              cbind.data.frame('mean'=mean(df$bpsys_high), 'sd'=sd(df$bpsys_high)),
                              cbind.data.frame('mean'=mean(df$bpdias_high), 'sd'=sd(df$bpdias_high)))
row.names(scale_info)=c('momage', 'bpsys_high', 'bpdias_high')

df = df %>% mutate(momage_sc = as.numeric(scale(momage))) %>%
   mutate(bpsys_high_sc = as.numeric(scale(bpsys_high))) %>% 
   mutate(bpdias_high_sc = as.numeric(scale(bpdias_high)))

# Separate training and test observations. 
# Stratified 80/20 split, by pec, gest_age, maybe_rounded.  
# Saved in dataset for reproducibility; was created by using df$testtrain = factor(testtrain(df)).
train = df %>% filter(testtrain=="train")
test = df %>% filter(testtrain=="test")

# Choose x covariates for analysis.
xcols = c('pec','femaleinfant','momage_sc','momjob','ganda','hiv','primip')

#======================================================================
# Investigate and plot binning.
#======================================================================

head(df)
summary(df$bw_kg); summary(df$bw_g)
table(df$maybe_rounded); round(prop.table(table(df$maybe_rounded)),3)

# Visualize binning in histogram.
ggplot(df, aes(x=bw_kg)) + 
   geom_histogram(bins=250) + 
   theme_Publication()  +
   scale_x_continuous(breaks=seq(0,5,by=1)) +
   labs(x='Birth weight (kg)', y='Count')
ggsave('./output/figure-01.pdf', width=6, height=4, dpi=800)
ggsave('./output/figure-01.tiff', width=6, height=4, dpi=800)

# Investigate rough size of sig2 compared to rounding (to nearest 100 grams).
mylm = lm(bw_kg ~ pec + gest_age28 + femaleinfant + momage + momjob + ganda + hiv +
             protein + primip, df)
sigma(mylm)
sd(df$bw_kg)


#======================================================================
# Calculate x-star, the "centroid x" patient, where patient takes on the 
# average or most common value of each covariate.  Set this up over
# all gest_ages, for pec=1 and pec=0, and add to test grid.
#======================================================================

xstar = df[1,]
xstar$id = 0
xstar$bw_kg = xstar$bw_g = xstar$maybe_rounded = NA
xstar$femaleinfant = median(df$femaleinfant)
xstar$momage = mean(df$momage)
xstar$momjob = names(which.max(table(df$momjob)))
xstar$ganda = median(df$ganda)
xstar$hiv = median(df$hiv)
xstar$momage_decade = names(which.max(table(df$momage_decade)))
xstar$bpsys_high = mean(df$bpsys_high)
xstar$bpdias_high = mean(df$bpdias_high)
xstar$protein = median(df$protein)
xstar$testtrain = 'centroid'

xstar = xstar[rep(rownames(xstar), each=2*length(unique(df$gest_age))), ]
xstar$pec = rep(c(0,1), each=nrow(xstar)/2)
xstar$gest_age = rep(sort(unique(df$gest_age)), times=2)
xstar$gest_age28 = xstar$gest_age - 27

#======================================================================
# Create grid of hypothetical test patients, to add to test set, for oos plotting.
# This will consist of the entire dataset, expanded for all time points.
#======================================================================

ptpanel = df[rep(row.names(df), each=length(unique(df$gest_age))),]
ptpanel$testtrain="ptpanel"
ptpanel$gest_age = rep(sort(unique(df$gest_age)), times=nrow(df))
ptpanel$gest_age28 = ptpanel$gest_age-27
ptpanel$bw_kg = NA
ptpanel$bw_g = NA
ptpanel$maybe_rounded = NA

# Set pec=1 for all obs.
ptpanel$pec=1

#  Make a second pt panel where we set pec=0 for all obs.
ptpanel_pec0 = ptpanel
ptpanel_pec0$pec = 0 

test = smartbind(test, xstar, ptpanel, ptpanel_pec0)
rm(ptpanel, ptpanel_pec0, xstar)


###########################################################################
# 2. Table 1, cohort description by pec/no pec.
###########################################################################

make_table1 = function(df, colvar){
   #--------------------------------------------------------------
   # FUNCTION: Creates a Table 1 of cohort characteristics, using
   #           the package qwraps2.  Note: Need to manually update
   #           the summaries object for each new data frame.
   #--------------------------------------------------------------
   # INPUTS:   df = the data frame
   #           colvar = the column variable
   #--------------------------------------------------------------
   # OUTPUT:   a data frame containing a cohort chars table
   #--------------------------------------------------------------
   require(qwraps2)
   
   # Set up column variables.
   df$trt = df$pec
   df$t = df$gest_age
   
   # Then proceed as usual.
   table = data.frame(matrix(NA,nrow=0,ncol=4))
   colnames(table) = c('Characteristic', 
                       paste0('Cohort (N = ',nrow(df),')'), 
                       paste0('No pecpsia (N = ', sum(df$trt==0),')'),
                       paste0('pecpsia (N = ', sum(df$trt==1),')')
   )
   
   require(rlist)
   require(tidyverse)
   require(qwraps2)
   options(qwraps2_markup = 'markdown')
   
   #--------------------------------------------------------------
   # Set up p-values vector.
   #--------------------------------------------------------------
   p = c(
      chisq.test(df[paste(colvar)], df$t)$p.value,
      t.test(df[paste(colvar)], df$momage)$p.value,  ### t-test
      chisq.test(df[paste(colvar)], df$momjob)$p.value,
      chisq.test(df[paste(colvar)], df$infantsex)$p.value,
      chisq.test(df[paste(colvar)], df$ganda)$p.value,
      chisq.test(df[paste(colvar)], df$hiv)$p.value,
      t.test(df[paste(colvar)], df$bpsys_high)$p.value, ### t-test
      t.test(df[paste(colvar)], df$bpdias_high)$p.value, ### t-test
      chisq.test(df[paste(colvar)], df$protein)$p.value,
      chisq.test(df[paste(colvar)], df$primip)$p.value
   )
   
   p = frmtp(p)
   p = gsub("\\*","",p)
   p = gsub("P ","",p)
   
   #--------------------------------------------------------------
   # Set up summaries info.
   #--------------------------------------------------------------
   summaries <-
      list("Gestational age in weeks, n (%)" = 
              list("28-31"       = ~ qwraps2::n_perc0(t %in% c(28:31), digits=2),
                   "32-37"         = ~ qwraps2::n_perc0(t %in% c(32:37), digits=2),
                   "38-39"         = ~ qwraps2::n_perc0(t %in% c(38,39), digits=2),
                   "40-42"         = ~ qwraps2::n_perc0(t %in% c(40:42), digits=2),
                   "43"         = ~ qwraps2::n_perc0(t %in% c(43), digits=2)
              ),
           "Maternal age in years, n (%)" =
              list("(13,20]"   = ~ qwraps2::n_perc0(momage_decade=='(13,20]', digits=2),
                   "(20,30]"   = ~ qwraps2::n_perc0(momage_decade=='(20,30]', digits=2),
                   "(30,40]"   = ~ qwraps2::n_perc0(momage_decade=='(30,40]', digits=2),
                   "(40,46]"   = ~ qwraps2::n_perc0(momage_decade=='(40,46]', digits=2)
              ),
           "Maternal job type, n (%)" = 
              list("Skilled"      = ~ qwraps2::n_perc0(momjob=="high", digits=2),
                   "Unskilled"    = ~ qwraps2::n_perc0(momjob=="med", digits=2),
                   "Unemployed"   = ~ qwraps2::n_perc0(momjob=="low", digits=2)
              ),
           "Infant sex, n (%)" = 
              list("Female"         = ~ qwraps2::n_perc0(femaleinfant==1, digits=2),
                   "Male"         = ~ qwraps2::n_perc0(femaleinfant==0, digits=2)
              ),
           "Ganda ethnicity, n (%)" = 
              list("No"         = ~ qwraps2::n_perc0(ganda==0, digits=2),
                   "Yes"         = ~ qwraps2::n_perc0(ganda==1, digits=2)
              ),
           "Maternal HIV, n (%)" = 
              list("Not Positive"         = ~ qwraps2::n_perc0(hiv==0, digits=2),
                   "Positive"         = ~ qwraps2::n_perc0(hiv==1, digits=2)
              ),
           "Max Systolic BP Delivery BP" = 
              list("mean pm sd"      = ~ qwraps2::mean_sd(bpsys_high, digits=2)
              ),
           "Max Diastolic Delivery BP" = 
              list("mean pm sd"      = ~ qwraps2::mean_sd(bpdias_high, digits=2)
              ),
           "Max Proteinuria, n (%)" = 
              list("0"         = ~ qwraps2::n_perc0(protein==0, digits=2),
                   "1"         = ~ qwraps2::n_perc0(protein==1, digits=2),
                   "2"         = ~ qwraps2::n_perc0(protein==2, digits=2),
                   "3"         = ~ qwraps2::n_perc0(protein==3, digits=2),
                   "4"         = ~ qwraps2::n_perc0(protein==4, digits=2)
              ),
           "Parity, n (%)" = 
              list("Multiparous"  = ~ qwraps2::n_perc0(primip==0, digits=2),
                   "Primiparous" = ~ qwraps2::n_perc0(primip==1, digits=2)
              )
   
      )
   
   #------------------------
   ### Create table.
   #------------------------
   
   # Set up basic table data.
   table = cbind(
      summary_table(df, summaries), # For entire df
      summary_table(dplyr::group_by(df, get(colvar)), summaries)) # For newpet columns
   
   # Extract elements needed to assemble table.
   grps = names(attr(unclass(table),'rgroups'))
   grp_nrows = attr(unclass(table),'rgroups')
   row_names = row.names(table)
   
   temp = unclass(table); row.names(temp) = NULL
   temp = as.data.frame(temp)
   
   # Assemble table. Begin with empty table1 dataframe. Loop through groups.
   table1 = cbind.data.frame('Characteristic' = as.character(row_names[0]), temp[0,])
   
   # Cast all rows as characters.
   for(i in 1:ncol(table1)){
      table1[,i] = as.character(table1[,i])
   }
   
   # Keep track of which rows we are in, as we loop through groups.
   row_marker = 1
   
   for(i in 1:length(grps)){
      
      # Add header row for group i.
      table1[nrow(table1)+1,] = ""
      table1[nrow(table1),"Characteristic"] = as.character(grps[i])
      
      # Add values for group i.
      idx = row_marker:(row_marker+grp_nrows[i]-1)
      table1 = rbind.data.frame(table1,
                                cbind.data.frame("Characteristic"=row_names[idx],
                                                 temp[idx,]))
      # Update row marker.
      row_marker = row_marker + grp_nrows[i] 
   }
   
   #-------------------------------
   # Format row and column names.
   #-------------------------------
   
   # Row names.
   row.names(table1) = NULL
   
   # The rest of this is all for column names.
   
   require(stringr)
   colnms = c('Characteristic',colnames(table))
   colnms = str_replace(colnms,'df ','')
   
   # Because str_replace not working on "get(colvar)" for some reason.
   for(i in 3:length(colnms)){
      colnms[i] = c(paste0(colvar,":",substr(colnms[i],13,nchar(colnms[i])) ))
   }
   
   colnames(table1) = colnms
   
   colnames(table1)[3] = gsub("pec: 0", "No pecpsia",colnames(table1)[3])
   colnames(table1)[4] = gsub("pec: 1", "pecpsia",colnames(table1)[4])
   
   table1[26,2] =  gsub("&plusmn;","&pm",table1[26,2])
   table1[26,3] =  gsub("&plusmn;","&pm",table1[26,3])
   table1[26,4] =  gsub("&plusmn;","&pm",table1[26,4])
   table1[28,2] =  gsub("&plusmn;","&pm",table1[28,2])
   table1[28,3] =  gsub("&plusmn;","&pm",table1[28,3])
   table1[28,4] =  gsub("&plusmn;","&pm",table1[28,4])


   #-------------------------------
   # Add p-values to table.
   #-------------------------------
   table1["P-value"] = ""
   table1[which(table1[,2]==""),"P-value"] = p
   
   #-------------------------------
   # Create vector of rows for hline argument when using xtable to print.
   #-------------------------------
   hlines = c(sort(c(which(table1[,2]=="")-1,which(table1[,2]==""))),nrow(table1))
   
   # Return table and hlines.
   return(list('table1'=table1,'hlines'=hlines))
}

# Create table.
tab1 = make_table1(df,colvar='pec')
captn = "Descriptive characteristics of dataset."

# Print and save table to .txt file.
print(xtable(tab1$table1, caption=captn,align=c('llrrrr')),
     include.rownames=FALSE,
     scalebox=.7,
     caption.placement='bottom',
     hline.after = tab1$hlines,
     file='./output/table-01.txt')


###########################################################################
# 3. Model Fitting.
###########################################################################

#======================================================================
# Fit model with tsbart, using monotonicity and treatment for rounded bwts.
#======================================================================

# Set up model matrices and cutpoints.
xx = tsbart::makeModelMatrix(train[,xcols, drop=F])
xxtest = tsbart::makeModelMatrix(test[,xcols, drop=F])

nsim = 10000; nburn = 100; ntree = 200;

# Fit using regular tsbart.
fit = tsbart(y=train$bw_kg, tgt=train$gest_age28, x=xx, 
             tpred=test$gest_age28, xpred=xxtest, 
             nsim=nsim, nburn=nburn, ntree=ntree)

# Fit with monotonicity.
fit_pava = tsbart(y=train$bw_kg, tgt=train$gest_age28, x=xx, 
                  tpred=test$gest_age28,xpred=xxtest, 
                  nsim=nsim, nburn=nburn, ntree=ntree, monotone="incr")

#======================================================================
# Save posterior means and credible intervals for each model, for in-sample and out-of-sample.
#======================================================================

# Non-monotone, in-sample.
train$yhat = colMeans(fit$mcmcdraws)
train$lb = apply(fit$mcmcdraws, 2, function(x) quantile(x,0.025))
train$ub = apply(fit$mcmcdraws, 2, function(x) quantile(x,0.975))

# Non-monotone, out-of-sample.
test$yhat = colMeans(fit$mcmcdraws_oos)
test$lb = apply(fit$mcmcdraws_oos, 2, function(x) quantile(x,0.025))
test$ub = apply(fit$mcmcdraws_oos, 2, function(x) quantile(x,0.975))

# Monotone, in-sample.
train$yhat_pava = colMeans(fit_pava$mcmcdraws)
train$lb_pava = apply(fit_pava$mcmcdraws, 2, function(x) quantile(x,0.025))
train$ub_pava = apply(fit_pava$mcmcdraws, 2, function(x) quantile(x,0.975))

# Monotone, out-of-sample.
test$yhat_pava = colMeans(fit_pava$mcmcdraws_oos)
test$lb_pava = apply(fit_pava$mcmcdraws_oos, 2, function(x) quantile(x,0.025))
test$ub_pava = apply(fit_pava$mcmcdraws_oos, 2, function(x) quantile(x,0.975))

#======================================================================
# Plot results: Figure 1, in-sample results and posterior means.
#======================================================================

#-------------------------------------------------------------------
# Mulago Figure 1: Panels A, B. Plot of in-sample results - scatter plots with posterior means for centroid.
#-------------------------------------------------------------------

# Arrange data in long form for plotting.
ggdf = train %>% 
   select(id, gest_age, pec, yhat_pava) %>% 
   mutate(pec = ifelse(pec==1,"PE","No PE"))

# Posterior mean for 'centroid' patient, for pec/no pec cases.
pmeans = test %>% 
   filter(testtrain=='centroid') %>%
   select(id, gest_age, pec, yhat_pava) %>%
   mutate(pec = ifelse(pec==1,"PE","No PE"))

# Add posterior mean for 'centroid' info to ggdf.
ggdf = dplyr::left_join(ggdf, pmeans %>% filter(pec=="PE") %>% select(pec,gest_age,yhat_pava), by=c('gest_age')) 
ggdf = dplyr::left_join(ggdf, pmeans %>% filter(pec=="No PE") %>% select(pec,gest_age,yhat_pava), by=c('gest_age')) 
ggdf = ggdf %>% select(id, gest_age, pec.x, yhat_pava.x, yhat_pava.y, yhat_pava) %>%
   rename(pec=pec.x, yhat_pava=yhat_pava.x, centroid_pec1 = yhat_pava.y, centroid_pec0 = yhat_pava)

### Plot panel of in-sample results.

p1 = ggplot(ggdf %>% filter(pec=="No PE"), aes(x=gest_age, y=yhat_pava)) + 
   geom_point(alpha=.5, colour='grey60') + 
   geom_line(aes(y=centroid_pec1), colour='blue', size=1) + 
   geom_line(aes(y=centroid_pec0), colour='forestgreen', size=1) + 
   labs(x=' ', y='Birth weight (kg)', subtitle = "A) No PE") +
   theme_Publication()

p2 = ggplot(ggdf %>% filter(pec=="PE"), aes(x=gest_age, y=yhat_pava)) + 
   geom_point(alpha=.5, colour='grey60') + 
   geom_line(aes(y=centroid_pec1), colour='blue', size=1) + 
   geom_line(aes(y=centroid_pec0), colour='forestgreen', size=1) + 
   labs(x='Gestational age (Wks)', y=' ', subtitle = "B) PE") +
   theme_Publication()

#-------------------------------------------------------------------
# Mulago Figure 1: Panel C. Plot of in-sample results - estimated PEC differences for everybody in dataset.
#-------------------------------------------------------------------d

# Calculate differences over gestation for each trajectory in pt panel (all pts in df, regardless of pec status).
ggdf_pec1 = test %>% filter(testtrain=='ptpanel', pec==1)
ggdf_pec0 = test %>% filter(testtrain=='ptpanel', pec==0)

ggdf = ggdf_pec1 %>% select(id, gest_age, yhat_pava) %>% rename(yhat_pec1=yhat_pava)
ggdf$yhat_pec0 = ggdf_pec0$yhat_pava
ggdf$pecdiff = ggdf$yhat_pec0 - ggdf$yhat_pec1
 
# Mean difference over gestation for centroid.
pmeans_1 = pmeans %>% filter(pec=="PE")
pmeans_0 = pmeans %>% filter(pec=="No PE")
pmeans = pmeans_1; pmeans$pecdiff = pmeans_0$yhat_pava - pmeans_1$yhat_pava
ggdf$pecdiff_mean = left_join(ggdf, pmeans, by='gest_age')$pecdiff.y

# Create plot.
p3 = ggplot(ggdf, aes(x=gest_age, y=pecdiff, colour=factor(id))) + 
   geom_line(alpha=.05) + 
   geom_line(aes(y=pecdiff_mean), colour='black',size=1) +
   scale_colour_manual(name='', values=rep('grey60', times=length(unique(ggdf$id)))) +
   guides(colour=F) + 
   coord_cartesian(ylim=c(0,1)) +
   labs(x=' ', y=' ', subtitle = "C) PE Decrease") +
   theme_Publication()

#-------------------------------------------------------------------
# Put panel together and clean up plot-specific helper data frames.
#-------------------------------------------------------------------
plt = grid.arrange(p1, p2, p3, ncol=3)
ggsave('./output/figure-05.pdf', plot(plt), width=9, height=3, dpi=800)
ggsave('./output/figure-05.tiff', plot(plt), width=9, height=3, dpi=800)

rm(ggdf, p1, p2, p3, pmeans, pmeans_0, pmeans_1, ggdf_pec0, ggdf_pec1)

#======================================================================
# Mulago Figure 2: Spaghetti plot of some trajectories, for monotone and non-monotone.
#======================================================================

# # Find pts with interesting differences for plotting.
# max(abs(test$yhat_pava - test$yhat))
# test$id[which(abs(test$yhat_pava - test$yhat)>.075 & test$gest_age<36)]

# Small set of patients, with monotone and non-monotone functions overlaid, to illustrate benefits.
idx = c(82, 2, 566)
ggdf = test %>% filter(id %in% idx) %>%
   mutate(fac_lab = paste0('Patient ',id))

# Plot trajectories in panel, for individual patients, with error bars, for mono and non-mono overlaid.
ggplot(ggdf %>% filter(pec==0), aes(x=gest_age, y=yhat_pava, colour=factor(id))) + 
   geom_ribbon(aes(ymin=lb_pava - mean(fit_pava$sigma), ymax=ub_pava + mean(fit_pava$sigma)), alpha=.8, colour='white', fill='grey80') + 
   geom_line(colour='black') +
   geom_line(aes(y=yhat), colour='black', linetype=2) + 
   facet_wrap(~fac_lab) + 
   coord_cartesian(expand=F) + 
   labs(x='Gestational age (Wks)', y='Birth weight (kg)') +
   theme_Publication()

ggsave('./output/figure-06.pdf', width=8, height=4, dpi=800)
ggsave('./output/figure-06.tiff', width=8, height=4, dpi=800)

rm(ggdf, idx)


###########################################################################
# 4. Posterior subgroup analysis.
###########################################################################

# # PEC is diagnosed based on sys_bp > 140 or dias_bp > 90, and protein>=1.  Create some covariates.
# train$sys_hi = ifelse(train$bpsys_high>=140, 1, 0)
# train$dias_hi = ifelse( train$bpdias_high>90, 1,0)
# train$either_bp_high = ifelse(train$bpsys_high>=140 | train$bpdias_high>90, 1,0)

### First, we know that pec --> early delivery.  We want to deconfound this by fitting
### the posterior CART model to the "residual" of each yhat_t estimate from its ybar_t, 
### ie how high or low an estimate is for its gestational age, instead of across the board.

### Heteroskedasticity matters here, too - to handle this, instead of modeling
### e_i = yhat(x_i,t_i) - ybar(t_i), we can model
### e_i = log[yhat(x_i,t_i) / ybar(t_i)]
### which has the effect of modeling log percent differences conditional on time.

ybar_t = train %>% dplyr::group_by(gest_age) %>%
   summarize('ybar_t' = mean(bw_kg))

# Create fit residuals as described above: r_i = log[yhat(x_i,t_i) / ybar(t_i)]
train = dplyr::left_join(train, ybar_t, by='gest_age')
train$resid = log(train$yhat_pava / train$ybar_t)

#======================================================================
# Fit CART model to posterior rel risk estimates to look at subgroups.
#======================================================================

# Fit CART tree to time-conditional residuals.
cart = rpart(resid ~ pec + femaleinfant + momage + momjob + ganda + hiv + primip,
             data=train, method='anova',model=T, control=list(maxdepth=3, cp=0))

tree = rpart.plot(cart,type=2)
rpart.plot(cart,type=4, nn.col='blue')

# Plot subgroups tree.
labs = paste(c(1:8), collapse="           ")

pdf('./output/figure-07.pdf', height=4, width=6)
   par(mfrow=c(2,1), mar=c(0,0,0,0)); layout(c(1,2), widths=c(1,1), heights=c(2,.125))
   rpart.plot(cart,type=2)
   
   plot(c(0, 1), c(0, .2), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
   text(x = 0.51, y = .15, labs,
        cex = 1.1, col = "black", family="sans", font=2, adj=0.5)
dev.off()

tiff('./output/figure-07.tiff', height=4, width=6, units='in', res=800)
   par(mfrow=c(2,1), mar=c(0,0,0,0)); layout(c(1,2), widths=c(1,1), heights=c(2,.125))
   rpart.plot(cart,type=2)
   
   plot(c(0, 1), c(0, .2), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
   text(x = 0.51, y = .15, labs,
        cex = 1.1, col = "black", family="sans", font=2, adj=0.5)
dev.off()s

#======================================================================
# Violin plots for CART model posterior exploration across tree levels.
#======================================================================

#---------------------------------------------------------
# Top level.
#---------------------------------------------------------

cart = rpart(resid ~ pec + femaleinfant + momage + momjob + ganda + hiv + primip,
             data=train, method='anova',model=T, control=list(maxdepth=1))
labels = factor(ifelse(cart$where==2, "1","2"), levels=c("1","2"))
grp_summary = getSubgroupInfo(fit_pava, tgt = train$gest_age, subgroups = labels)$post_mean_draws_subgrps

p1 = ggplot(grp_summary %>% gather(value=postmean, key='subgroup'), 
       aes(x=subgroup, y=postmean)) + 
   geom_violin(fill='grey80') + 
   labs(x='Node', y='Birth weight (kg)', subtitle='A) First level.') +
   theme_Publication()
p1

#---------------------------------------------------------
# Middle level.
#---------------------------------------------------------

cart = rpart(resid ~ pec + femaleinfant + momage + momjob + ganda + hiv + primip,
             data=train, method='anova',model=T, control=list(maxdepth=2))
labels = ifelse(cart$where==3, "1", ifelse(cart$where==4,"2", ifelse(cart$where==6,"3", "4")))
grp_summary = getSubgroupInfo(fit_pava, tgt = train$gest_age, subgroups = labels)$post_mean_draws_subgrps

p2 = ggplot(grp_summary %>% gather(value=postmean, key='subgroup'), 
            aes(x=subgroup, y=postmean)) + 
   geom_violin(fill='grey80') + 
   labs(x='Node', y='Birth weight (kg)', subtitle='B) Second level.') +
   theme_Publication()
p2

#---------------------------------------------------------
# Bottom level.
#---------------------------------------------------------

cart = rpart(resid ~ pec + femaleinfant + momage + momjob + ganda + hiv + primip,
             data=train, method='anova',model=T, control=list(maxdepth=3, cp=0))
label_dict = cbind.data.frame('leaf'=as.numeric(names(table(cart$where))), 'leafno'=rank(as.numeric(names(table(cart$where)))))
labels = factor(dplyr::left_join(data.frame('leaf'=cart$where), label_dict, by='leaf')$leafno)
grp_summary = getSubgroupInfo(fit_pava, tgt = train$gest_age, subgroups = labels)$post_mean_draws_subgrps

p3 = ggplot(grp_summary %>% gather(value=postmean, key='subgroup'), 
            aes(x=subgroup, y=postmean)) + 
   geom_violin(fill='grey80') + 
   labs(x='Node', y='Birth weight (kg)', subtitle='B) Third level.') +
   theme_Publication()
p3

#---------------------------------------------------------
# Add plots together to make panel.
#---------------------------------------------------------
plt = grid.arrange(p1, p2, p3, ncol=3)
ggsave('./output/figure-08.pdf', plot(plt), width=10, height=4, dpi=800)
ggsave('./output/figure-08.tiff', plot(plt), width=10, height=4, dpi=800)



