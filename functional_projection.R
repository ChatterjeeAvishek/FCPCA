functional_projection=function(fdobj1,basis,class_mean,mean_deduction){
  
  #fdoj1 is a functional object i.e, the class of the functional data
  # basis is the new direction where we want to project the data
  # class_mean is the mean of the class
  # mean_deduction is either TRUE or FALSE
  
  dime=dim(basis$coefs)
  k=dime[2]
  
  
  dim1=dim(fdobj1$coefs)
  
  proj1=t(rep(1, times=k))
  
  
  tempfd1_projected=fdobj1
  basis=basis
  time1_array=0
  
 
    # 
    # for(i in 1:dim1[2]){
    #   dmean_class1=fdobj1[i]-class_mean
    #   start_time <- Sys.time()
    #   proj1=rbind(proj1,inprod(dmean_class1,basis))
    #   end_time <- Sys.time()
    #   time1_array[i]=(end_time - start_time)
    # }
    # ####
  
  dmean_class1=fdobj1[1]-class_mean
  jin=2
  start_time <- Sys.time()
  for(i in 2:dim1[2]){
    dmean_class1$coefs= cbind(dmean_class1$coefs,(fdobj1[i]-class_mean)$coefs)
    
  }
  #dmean_class1$coefs=A_temp
  
  
  proj1=inprod(dmean_class1,basis)
  end_time <- Sys.time()
  time1=(end_time - start_time)
  
  
 
   #####
  
  #proj1=proj1[-1,]
  proj1=as.matrix(proj1)
  proj1_S1=proj1
  start_time <- Sys.time()
  tempfd1_projected_mat=basis$coefs%*%t(proj1)
  end_time <- Sys.time()
  time1=time1+(end_time - start_time)
  
  tempfd1_projected$coefs=tempfd1_projected_mat
  result=list('coef_class'=proj1,'tempfd1_projected'=tempfd1_projected,'time_projection'=time1)
  return(result)
  #return(tempfd1_projected)
  
}