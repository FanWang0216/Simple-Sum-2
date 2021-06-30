##This Rscript will provide sample code for conducting the SS2 analysis
setwd('/Users/macbook/work/error-in-variable/SS2_papercode/')
library('CompQuadForm')
##sample gwas file  
gwas_file=read.table(file='sampledata/GWAS_file_gene_MUC4_tissue_Human_nasal_epithelial.txt',header=T)
gwas_pvalue=gwas_file$ALL.FIXED.PVAL
##sample eqtl file
eqtl_file=read.table(file='sampledata/eqtl_file_gene_MUC4_tissue_Human_nasal_epithelial.txt',header=T)
eqtl_pvalue=eqtl_file$pval
##covariance of the GWAS summary statistics
ld_gwas=as.matrix(read.table(file='sampledata/ldGWAS_file_gene_MUC4_tissue_Human_nasal_epithelial.txt',header=F))
##covariance of the eQTL summary statistics
ld_eqtl=as.matrix(read.table(file='sampledata/ldeqtl_file_gene_MUC4_tissue_Human_nasal_epithelial.txt',header=F))

source('SS2_Rfunctions.R')
##performing SS2 analysis by using sample dataset
stage1_pv=stage1(eqtl_pvalue,ld_eqtl);stage1_pv
##if significant, perform stage 2
stage2_pv=stage2(eqtl_pvalue,gwas_pvalue,ld_gwas,use_evidence='log10p');stage2_pv

