########################################################################
# getSubgroupInfo: Function to get posterior draws, means, and bounds 
#                  from subgroups of data.
#-----------------------------------------------------------------------
# INPUTS: my_output: A list; the result of the tsbart() function from tsbart package.
#         subgroups: An optional n-length vector of subgroups.
#         subgroups_pred: An optional out-of-sample n-length vector of subgroups for prediction.
#-----------------------------------------------------------------------
# OUTPUTS:
#   post_mean_draws_subgrps: a (nsim x ngroups) matrix of posterior mean draws for each subgroup.  (No tgt info.)
#   post_mean_t_draws: a (nsim x (ngroups x tlen)) matrix of posterior mean draws for each subgroup/tgt combo.
#   post_mean_t: A data frame where each row is a subgroup/tgt combo, and columns are:
#       tgt, subgroup, yhat (posterior mean), lb, ub (95% credible intervals)

########################################################################

getSubgroupInfo = function(my_output, tgt, subgroups){
   
   # Set up grids of times and groups.
   tgrid = sort(unique(tgt))
   grps = sort(unique(subgroups))

   # Number of times and groups.
   nt = length(tgrid)
   ng = length(unique(subgroups))

   # Empty list for output.
   out = list()
   
   # Obtain input matrix.
   my_matrix = my_output$mcmcdraws

   #-----------------------------------------------------------------------
   # a. Posterior draws for overall CATE (using subgroups).
   #-----------------------------------------------------------------------
   
   # In-sample.
   cate = as.data.frame(matrix(0,nrow=nrow(my_matrix),ncol=ng))
   colnames(cate) = grps
   
   for(g in 1:ng){
      cate[,g] =  rowMeans(my_matrix[,which(subgroups==grps[g])])
   }
   
   out$post_mean_draws_subgrps = cate
   
   #-----------------------------------------------------------------------
   # . Posterior draws for conditional mean at each time (using subgroups.)
   #-----------------------------------------------------------------------
   
   # In-sample.
   cpost_mean = as.data.frame(matrix(0,nrow=nrow(my_matrix),ncol=ng*nt))
   
   groups_and_times = expand.grid(grps,tgrid)
   colnames(cpost_mean) = do.call(paste0, expand.grid(grps,'-',tgrid))
   
   
   
   for(tg in 1:(ng*nt)){
      temp = my_matrix[,which(subgroups==groups_and_times[tg,1]
                              & tgt==groups_and_times[tg,2]), drop=F]
      
      cpost_mean[,tg] =  apply(temp,1,function(x) mean(x,na.rm=T))
   }
   
   out$post_mean_t_draws = cpost_mean
   
   
   #-----------------------------------------------------------------------
   # c. Posterior means and credible intervals for conditional posterior means at each tgt.
   #-----------------------------------------------------------------------
   
   # In-sample.
   out$post_mean_t = cbind.data.frame('tgt'=rep(tgrid, each=ng),
                                      'subgroup'=rep(grps, times=nt),
                                      'yhat'=apply(out$post_mean_t_draws,2,function(x) mean(x,na.rm=T)),
                                      'lb'=apply(out$post_mean_t_draws,2,function(x) quantile(x,.025,na.rm=T)),
                                      'ub'=apply(out$post_mean_t_draws,2,function(x) quantile(x,.975,na.rm=T)))
   
   ########################################################################
   # Return function output.
   ########################################################################
   return(out)
}
