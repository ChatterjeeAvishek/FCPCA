# FCPCA
## Instructions to run functional class wise principal component analysis: 

The codes are developed in R (Version: 3.5.1)
The following three functions are attached in the repository: 
1.	FCPCA.R
2.	Gram_Schmidt_orthonormalization.R
3.	functional_projection.R
Required Packages: “fda”, “MASS”, “tictoc”

### Description: 
The function FCPCA.R performs functional classwise principal component analysis. It uses the above mentioned two functions Gram_Schmidt_orthonormalization.R and functional_projection.R. The Gram_Schmidt_orthonormalization.R function takes the functional basis object as input and performs Gram Schmidt Orthonormalization process for functional data. The functional_projection.R function takes functional data matrix and the basis functions and projects the functional data with respect to the basis functions. Usage along with arguments of the FCPCA function is given below. 

### Use: 

FCPCA(Train_data_matrix, Test_data_matrix, ClassInfo_column_no, nbasis, norder, No_fpca_projections) 

### Inputs:

Train_data_matrix: Train data matrix in which any one column represents train class labels. 

Test_data_matrix: Test data matrix in which any one column represents train class labels. 

ClassInfo_column_no: Mention the column number in which both train and test data have the class label informations.

nbasis: Number of b-spline basis functions to consider.

norder: An integer that spacifies the order of b-splines

No_fpca_projections: Vector containg number of eigenfunctions to use for subspace projections : example: c(1,1)
