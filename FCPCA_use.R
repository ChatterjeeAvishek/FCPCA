library(fda)
library(MASS)
library(tictoc)
datasetname="Beef"
source('F:/code final/FCPCA.R')
source('F:/code final/functional_projection.R')
source('F:/code final/Gram_Schmidt_orthonormalization.R')
classification_rate_array_fcpca=0
for (i in 1:100) {
  
  string_train=paste(datasetname,"_TRAIN_",as.character(i),".csv",sep="")
  string_test=paste(datasetname,"_TEST_",as.character(i),".csv",sep="")
  C_train=read.csv(string_train,header = FALSE)
  C_test=read.csv(string_test,header = FALSE)
  
  C_train=data.matrix(C_train)
  C_test=data.matrix(C_test)
  
  C_train=C_train[,-2]
  C_test=C_test[,-2]
  ntrain=nrow(C_train)
  ntest=nrow(C_test)
 
  cr1=FCPCA(C_train, C_test,ClassInfo_column_no=1,nbasis_bspline=33,norder=6,No_fpca_projections=c(3,3,3,3,3))
  classification_rate_array_fcpca[i]=cr1
  
}

mean(classification_rate_array_fcpca)
sd(classification_rate_array_fcpca)
