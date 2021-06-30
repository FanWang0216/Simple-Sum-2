##This Rscript contains functions for calculating the covariance matrix of Z-scores 
getcov_Zmeta<-function(weightmatrix,varz_subgroup){
  num_study=dim(weightmatrix)[2]
  weightld_group=list()
  for (i in c(1:num_study)){
    W_snp_cohort=weightmatrix[,i]
    sqrtW_snp_cohort=sqrt(W_snp_cohort)
    matrix_weight=sqrtW_snp_cohort%*%t(sqrtW_snp_cohort)
    weightld_group[[i]]=(matrix_weight)*(varz_subgroup[[i]])
  }
  sum_ld=Reduce('+',weightld_group)
  vec_allstudy=rowSums(weightmatrix)
  sqrt_sumW=sqrt(vec_allstudy)
  factor_mat=1/(sqrt_sumW%*%t(sqrt_sumW))
  final_ld=factor_mat*sum_ld
  return(final_ld)
}

##Example:
setwd('/Users/macbook/work/error-in-variable/SS2_papercode/sampledata/')
load('13subgroup_ld.Rdata')
weightmatrix=read.table('weightmatrix.txt',header=F)
ld_meta=getcov_Zmeta(weightmatrix,varz_subgroup)
############################################################################

get_subgroup_var<-function(G,X,Sigma){
  m=dim(Sigma)[1];m
  #to make the Sigma invertible
  Sigma_est=Sigma+diag(rep(6.610696e-05,m))
  P=solve(Sigma_est)-solve(Sigma_est)%*%X%*%solve(t(X)%*%solve(Sigma_est)%*%X)%*%t(X)%*%solve(Sigma_est)
  num=t(G)%*%P%*%G
  deno1=solve(diag(sqrt(diag(t(G)%*%P%*%G))))
  covZ=deno1%*%num%*%deno1
  return(covZ)
}

##Example:

##simulate genotypes
get_genotype<-function(p,number_indi){
  genotype_col=rbinom(number_indi,2,p)
  return(genotype_col)
}
number_SNP=663
maf=runif(number_SNP,0.05,0.5)
set.seed(6666)
genotype_mat1=as.matrix(sapply(maf,FUN=get_genotype,70))
genotype_mat1_stand=scale(genotype_mat1,center=TRUE,scale=TRUE)
library('MASS')
G=genotype_mat1_stand
##Simulate additional covariates
X=mvrnorm(n=70,mu=rep(0,3),Sigma=diag(rep(1,3)))
Sigma_relateness=read.table(file='Sigma_relateness.txt',header=F)

covariance_Z=get_subgroup_var(G,X,Sigma_relateness)
