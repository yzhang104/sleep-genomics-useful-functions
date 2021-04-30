################################################
# The p-filter: multilayer FDR control for grouped hypotheses 
# source: Barber RF, Ramdas A. The p-filter: multi-layer FDR control for grouped hypotheses. arXiv:151203397 [stat] 2016; published online Oct 28.
pfilter = function(P,alphas,group){
  # P in [0,1]^n = vector of p-values
  # alphas in [0,1]^M = vector of target FDR levels
  # groups is a n-by-M matrix; 
  #	groups[i,m] = which group does P[i] belong to,
  #		for the m-th grouping
  
  # change groups to a matrix of col_1: number of rows, col_2: convert trait and metabolties to numeric value
  if (is.null(dim(group))) {
    groups<-as.matrix(data.frame(n_row=c(1:length(group)),n_group=as.numeric(as.factor(group))))
    alphas<-alphas[1:2]
  } else {
    groups<-as.matrix(data.frame(n_row=c(1:nrow(group)),n_group1=as.numeric(as.factor(group[,1])),n_group2=as.numeric(as.factor(group[,2]))))
  }
  
  n = length(P)
  M = length(alphas)
  G = apply(groups,2,max) # G[m] = # groups, for grouping m
  
  
  Simes = list()
  for(m in 1:M){
    Simes[[m]]=rep(0,G[m])
    for(g in 1:G[m]){
      group = which(groups[,m]==g)
      Simes[[m]][g]=min(sort(P[group])*length(group)/(1:length(group)))
    }
  }
  
  
  # initialize
  thresh = alphas
  Sh = 1:n
  for(m in 1:M){
    pass_Simes_m = which(is.element(groups[,m],which(Simes[[m]]<=thresh[m])))
    Sh = intersect(Sh,pass_Simes_m)
  }
  done = FALSE
  
  
  while(!done){
    thresh_old = thresh
    for(m in 1:M){
      # which groups, for the m-th grouping, 
      #	have any potential discoveries?
      Shm = sort(unique(groups[Sh,m]))
      
      # run BH on Simes[[m]], constraining to Shm
      Pvals_m = rep(1.01,G[m]); # >1 for groups not in Dm
      Pvals_m[Shm] = Simes[[m]][Shm]
      khatm = max(0,which(sort(Pvals_m)<=(1:G[m])/G[m]*alphas[m]))
      thresh[m] = khatm/G[m]*alphas[m]
      Sh = intersect(Sh,
                     which(is.element(groups[,m],which(Pvals_m<=thresh[m]))))
    }
    if(all(thresh_old==thresh)){done = TRUE}
  }
  
  Sh_temp = Sh;
  Sh = rep(0,n); Sh[Sh_temp] = 1
  Sh
  
}


# function to loop pfilter over multiple models and strata
pfilter_loop<-function(data, model_input, strata_input,alpha_threshold){
  output_data<-data.frame()
  for (i in seq_along(model_input)){
    strata_output<-data.frame()
    for (j in seq_along(strata_input)){
      temp_data<-data%>%
        dplyr::filter(model==model_input[i],strata==strata_input[j])
      temp_data$pfilter_1<-pfilter(P=temp_data$p_val,alphas=alpha_threshold,group=temp_data$trait)
      temp_data$pfilter_2<-pfilter(P=temp_data$p_val,alphas=alpha_threshold,group=temp_data[,c("trait","metabolite")])
      strata_output<-rbind(strata_output,temp_data)
    }
    output_data<-rbind(output_data,strata_output)
  }
  return(output_data)
}
