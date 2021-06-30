##functions for performing SS2 analysis
#stage 1 test:
stage1<-function(eqtl_pval,ld_eqtl){
  doeigen_vareqtl=eigen(ld_eqtl, symmetric=TRUE, only.values = FALSE, EISPACK = FALSE)
  LD_eig<-doeigen_vareqtl$values
  statistic=sum(qnorm(eqtl_pval/2,lower.tail=F)^2)
  m=dim(ld_eqtl)[1]
  pv<-davies(statistic,LD_eig,h=rep(1,m),delta=rep(0,m),acc=0.0001)$Qq
  if (pv<=0.001){
    pv=abs(imhof(statistic,LD_eig,h=rep(1,m),delta=rep(0,m))$Qq)
  }
  return(pv)
}

#stage 2 test: 
#use_evidence='Tsquare' means using the square of eQTL Z-scores as eQTL evidence measure
#use_evidence='log10p' means using the -log10 (eQTL p-value) as eQTL evidence measure
get_A<-function(Tissue_vec,m){
  A=matrix(0,m,m)
  a_diag<-rep(0,m)
  for (j in 1:m){
    if (sum(Tissue_vec)==0){
      a_diag[j]<--1/m
    }else if (sum(Tissue_vec)==m){
      a_diag[j]<-1/m
    }else{
      Tbar<-mean(Tissue_vec)
      a_diag[j]=(Tissue_vec[j]-Tbar)/(sum(Tissue_vec^2)-m*(Tbar^2))
    }
  }
  diag(A)=a_diag
  return(A)
}
get_eigenmatrix<-function(sd_cor,Tissue_vec,m){
  result=list()
  matrix_A<-get_A(Tissue_vec,m)
  matrix_mid<-(t(solve(sd_cor)))%*%matrix_A%*%(solve(sd_cor))
  eigenvalues<-eigen(matrix_mid)$values
  eigenmatrix<-eigen(matrix_mid)$vectors
  result[[1]]=eigenvalues
  result[[2]]=eigenmatrix
  return(result)
}
get_simplesumstats<-function(wald,covariate,m){
  if (sum(covariate)==0){
    SS<--mean(wald)
  }else if (sum(covariate)==m){
    SS<-mean(wald)
  }else{
    reg<-lm(wald~covariate)
    SS<-summary(reg)$coefficients[2,1]
  }
  return(SS)
}
get_p<-function(m,eigen_value,teststats,method='imhof'){
  switch(method,imhof={
    pv<-imhof(teststats,eigen_value,h=rep(1,m),delta=rep(0,m))$Qq
  },
  davies={
    pv<-davies(teststats,eigen_value,h=rep(1,m),delta=rep(0,m))$Qq
  })
  return(pv)
}

stage2<-function(eqtl_pval,gwas_pval,ld_gwas,use_evidence='Tsquare'){
  switch(use_evidence,Tsquare={
    eqtl_evid=(qnorm(eqtl_pval/2,lower.tail=F))^2
  }, log10p={
    eqtl_evid=-log10(eqtl_pval)
  })
  chol_varz<-chol(solve(ld_gwas))
  wald=(qnorm(gwas_pval/2,lower.tail=F))^2
  m=length(wald)
  eig_org<-get_eigenmatrix(chol_varz,eqtl_evid,m)[[1]]
  ss_org<-get_simplesumstats(wald,eqtl_evid,m)
  pv_SS<-get_p(m,eig_org,ss_org,method='davies')
  if (pv_SS<=0.001){
    pv_SS<-abs(get_p(m,eig_org,ss_org,method='imhof'))
  }
  return(pv_SS)
}


