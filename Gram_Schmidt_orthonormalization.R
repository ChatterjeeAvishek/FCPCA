Gram_Schmidt_orthonormalization=function(grand_mean_fd,mean_mat,old_basis,nbasis,No_of_classes){
  
  # meanlogprec1=mean(fdobj1)
  # meanlogprec2=mean(fdobj2)
  # meanlogprec3=mean(fdobj3)
  
  # Grand mean Calculation 
  m1=grand_mean_fd
  m2=grand_mean_fd
  grand_mean_vec=m2$coefs
  k=nbasis
for (j in (1:(No_of_classes-1))) {
  m1$coefs=mean_mat[,j]
  v=m1-m2
  newbasis=old_basis
  A1=inprod(v,old_basis)
  A2=inprod(old_basis,old_basis)
  Icomp=0
  sum=0
  for (i in 1:k) {
    
    Icomp[i]=A1[i]/A2[i,i]
    
  }
  V_n=(Icomp*old_basis)
  for (i in 1:k){
    sum=sum+V_n$coefs[,k]
  }
  k=k+1
  sum=mean_mat[,j]-grand_mean_vec-sum
  #sum=sum/norm(sum)
  newbasis$coefs = t(rbind(t(old_basis$coefs),t(sum))) 
  old_basis=newbasis
}
#orthonormalization
A_orthononormal=inprod(newbasis,newbasis)
for (i in 1:ncol(A_orthononormal)){
newbasis$coefs[,i]=newbasis$coefs[,i]/sqrt(A_orthononormal[i,i])
}
return(newbasis)

} 
  

