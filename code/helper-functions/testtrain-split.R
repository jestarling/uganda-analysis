########################################################################
# Split test/train data based on stratified sample
########################################################################

testtrain = function(data, test_pct = .3){
   
   ### Stratified sampling based on all combos of fd/nd and diab_htn.
   
   # Library for left_join.
   require(dplyr)
   require(sampling)
   
   # Combine nd and fd into one temp variable for convenience.
   #grps = expand.grid('fd' = c(0,1), 
   #                   'diab_htn' = factor(c('Neither','Diabetes','Htn','Both'), 
   #                                       levels = c('Neither', 'Diabetes', 'Htn','Both')))
   
   grps = expand.grid('pec' = c(0,1),
                      'maybe_rounded'=c(0,1),
                      'gest_age'=28:43)
   
   grps$strata = 1:nrow(grps)
   
   # Match two columns to get strata group (out of possible 12 groups).
   # data$strata = suppressWarnings(left_join(data, grps, by=c('fd','diab_htn'))$strata)
   data$strata = suppressWarnings(left_join(data, grps, by=c('pec','maybe_rounded','gest_age'))$strata)
   
   # Vectors to hold row indices of test and train datasets.
   test_rowidx = NULL
   train_rowidx = NULL
   
   # Loop through strata, and assign to train/test set.
   for(i in 1:length(grps$strata)){
      
      # Unique ids in strata i.
      ids_in_strata = unique(data$id[which(data$strata==i)])
      
      # Sample size of test set.
      sampsize = round(length(ids_in_strata) * test_pct)
      
      test_ids = sample(ids_in_strata, sampsize, replace=F)
      train_ids = ids_in_strata[-which(ids_in_strata %in% test_ids)]
      
      test_rowidx = c(test_rowidx, which(data$id %in% test_ids))
      train_rowidx = c(train_rowidx, which(data$id %in% train_ids))
      
   }
   
   # Vector of test/train indicators.
   test_train = rep('train',nrow(data))
   test_train[test_rowidx] = 'test'
   
   # Return output
   return(test_train)
}
