FCPCA=function(Class_train, Class_test,ClassInfo_column_no,nbasis_bspline,norder,No_fpca_projections){
  
  #Class_train: Train data matrix in which any column represents train class labels
  #Class_test: Test data matrix in which any column represents train class labels
  #ClassInfo_column_no: Mention the column number in which both train and test data have the class label informations (Test and Train data can't have distinct column numbers)
  #nbasis_bspline: No of basis elements you want to represent the functional data
  #norder: Smoothing Parameter (preferably fix it at 6 but it may depend on the problem)
  #No_fpca_projections: Vector containg number of eigenfunctions to use for subspace projections : example: c(1,1)
  #Author: Avishek Chatterjee mail id: ac17rs012@iiserkol.ac.in
  
  # Finding total no of class presents in the training data set  
  Class_lebels=unique(Class_train[,ClassInfo_column_no])
  No_of_classes=length(Class_lebels)
  No_of_element_per_class_train=0
  No_of_element_per_class_test=0
  
  Class_time_points=ncol(Class_train[,-1])
  Time = 1:Class_time_points;
  Rng = c(1,Class_time_points);
  
  
  
  Class_lebels_test=unique(Class_test[,ClassInfo_column_no])
  for (fi in 1:length(Class_lebels_test)) {
     if( length(which(Class_test[,ClassInfo_column_no]==fi))== 1){
       pos=which(Class_test[,ClassInfo_column_no]==fi)
       Class_test=rbind(Class_test,Class_test[pos,])
     }
  }
  
  
  #Basis (Bspline) representation of the functional data
  
  Basis = create.bspline.basis(Rng, nbasis_bspline, norder)
  D2fdPar = fdPar(Basis, lambda=1e-7)
  mat=as.matrix(rep(1,nbasis_bspline))
  tic()
  for (i in 1:No_of_classes) {
    Class_Name=paste("Class_train", as.character(i),sep = "")
    s_train=which(Class_train[,ClassInfo_column_no]==Class_lebels[i])
    Class1_train=Class_train[s_train,-1]
    assign(Class_Name,Class1_train)
    No_of_element_per_class_train[i]=length(s_train)
    
    #Functional object for train data
    Class_Name_fn=paste("tempfd_train", as.character(i),sep = "")
    val=smooth.basis(Time, t(Class1_train), D2fdPar)$fd
    assign(Class_Name_fn,val)
    
    meanlogprec = mean(val)
    mvec=meanlogprec$coefs
    mat=matrix(c(mat,mvec),ncol=i+1)
    
    #FPCA on different classes
    temppca = pca.fd(val,nharm=No_fpca_projections[i])
    harmfd = temppca$harmonics
    FPCA_Name_fn=paste("harmfd", as.character(i),sep = "")
    assign(FPCA_Name_fn,harmfd)
    
    #Cumulative variance explained
    total_variance=sum(temppca$varprop)
    print(c("total_variance for Subspace",as.character(i),as.character(total_variance)))
    
    ########
    
    Class_Name=paste("Class_test", as.character(i),sep = "")
    s_test=which(Class_test[,ClassInfo_column_no]==Class_lebels[i])
    Class1_test=Class_test[s_test,-1]
    assign(Class_Name,Class1_test)
    No_of_element_per_class_test[i]=length(s_test)
    
    
    if(No_of_element_per_class_test[i]!=0){
    #Functional object for test data
    Class_Name_fn=paste("tempfd_test", as.character(i),sep = "")
    val=smooth.basis(Time, t(Class1_test), D2fdPar)$fd
    assign(Class_Name_fn,val)
    }
  }
  # Grand Mean Calculation
  mat=mat[,-1]
  mean_fdoj=tempfd_train1
  mean_fdoj$coefs=mat
  grand_mean=mean(mean_fdoj)
  grand_mean_vec=grand_mean$coefs
  posterior_prob=rep(1,nrow(Class_test))
  
  for (i in 1:No_of_classes) {
    FPCA_Name_fn=paste("harmfd", as.character(i),sep = "")
    harmfd=eval(parse(text = FPCA_Name_fn))
    newbasis_S=Gram_Schmidt_orthonormalization(grand_mean_fd=grand_mean,mean_mat=mat,old_basis=harmfd,nbasis=No_fpca_projections[i],No_of_classes)
    New_Basis_Name_fn=paste("newbasis_S", as.character(i),sep = "")
    assign(New_Basis_Name_fn,newbasis_S)
    
    Subspace_name=paste("S", as.character(i),sep = "")
    Class_name_fd_mean=paste("tempfd_train", as.character(i),sep = "")
    val1=eval(parse(text = Class_name_fd_mean))
    meanlogprec = mean(val1)
    
    feature_extracted_dim=No_fpca_projections[i]+(No_of_classes-1)
    proj_train=rep(1,feature_extracted_dim)
    proj_test=rep(1,feature_extracted_dim)
    
    y_train=1
    y_test=1
    
    
    for (j in 1:No_of_classes) {
      Class_name_fd=paste("tempfd_train", as.character(j),sep = "")
      tempfd=eval(parse(text = Class_name_fd))
      S_tempfd_projected=functional_projection(tempfd,newbasis_S,meanlogprec,TRUE)
      Class_name_fd_projected=paste(Subspace_name,"_tempfd", as.character(j),"_projected",sep = "")
      assign(Class_name_fd_projected,S_tempfd_projected)
      proj_train=rbind(proj_train,S_tempfd_projected$coef_class)
      
      if(No_of_element_per_class_test[j]!=0){
      
      Class_name_fd_test=paste("tempfd_test", as.character(j),sep = "")
      tempfd_test=eval(parse(text = Class_name_fd_test))
      S_tempfd_projected_test=functional_projection(tempfd_test,newbasis_S,meanlogprec,TRUE)
      Class_name_fd_projected_test=paste(Subspace_name,"_tempfd", as.character(j),"_projected_test",sep = "")
      assign(Class_name_fd_projected_test,S_tempfd_projected_test)
      proj_test=rbind(proj_test,S_tempfd_projected_test$coef_class)
      
      
      
      y_test=c(y_test,rep(j,No_of_element_per_class_test[j]))
      }
      y_train=c(y_train,rep(j,No_of_element_per_class_train[j]))
    }
    proj_train=proj_train[-1,]
    proj_test=proj_test[-1,]
    y_train=y_train[-1]
    y_test=y_test[-1]
    
    #Performing LDA
    data_train=data.frame(proj_train)
    data_train$response=y_train
    data_test=data.frame(proj_test)
    data_test$response=y_test 
    # Include library(MASS)
    model <- lda(response~., data = data_train)
    # Make predictions
    predictions <- predict(model,data_test)
    posterior_prob=cbind(posterior_prob,predictions$posterior)
    
  }
  
  posterior_prob=posterior_prob[,-1]
  colnames(posterior_prob)<-NULL
  rownames(posterior_prob)<-NULL
  
  index_mat=0
  correct_classification=0
  start_point=1
  sum=0
  
  for (i in 1:No_of_classes) {
    index_mat=seq(i,(No_of_classes^2), by = No_of_classes)
    sum=sum+No_of_element_per_class_test[i]
    if(No_of_element_per_class_test[i]!=0){
    posterior_mat=posterior_prob[start_point:sum,]
    for (j in 1:No_of_element_per_class_test[i]) {
      comp=which(is.element(index_mat,which.max(as.vector(posterior_mat[j,]))))
      if(length(comp)==1){
        correct_classification=correct_classification+1}
    }
    start_point=No_of_element_per_class_test[i]+start_point
    }
  }
  toc()
  cr=correct_classification/(nrow(Class_test))
  print(cr)
  return(cr)
}
