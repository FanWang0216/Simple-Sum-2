---
title: "Simple Sum 2 Colocalization Analysis"
output: 
  github_document:
        pandoc_args: --webtex
---
#Simple Sum 2 
Simple Sum 2 is a summary statistics-based colocalization tool that extends the frequently implemented **[Simple Sum method](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1008007)** (Gong et al, 2019) to a flexible procedure (SS2) that integrates GWAS summary statistics with eQTL summary statistics across any number of gene-by-tissue pairs, is applicable when there are overlapping participants in the two studies and can be applied to GWAS summary statistics computed through meta-analysis, even with related individuals. 

*   The **[sampledata](https://github.com/FanWang0216/SimpleSum2Colocalization/tree/main/sampledata)** folder provides all the sample datasets for conducting the Simple Sum 2 colocalization analysis, including the GWAS and eQTL files for two genes: _MUC20_ and _MUC4_ in tissue human nasal epithelial. 

* **[SS2_functions.R](https://github.com/FanWang0216/SimpleSum2Colocalization/blob/main/SS2_functions.R)** contains all the R functions for conducting the SS2 test which has been implemented in a web-based tool: **[LocusFocus](https://github.com/naim-panjwani/LocusFocus)**, with an additional R option for using the chi-square statsitics as eQTL evidence measure.

* **Commands** to run Simple Sum 2 for one GWAS/eQTL combination:  
 
    1. Obtain the **GWAS p-values** from the GWAS summary statistics file:      _gwas_file=read.table(file='sampledata/[GWAS_file_gene_MUC4_tissue_Human_nasal_epithelial.txt](https://github.com/FanWang0216/SimpleSum2Colocalization/blob/main/sampledata/GWAS_file_gene_MUC4_tissue_Human_nasal_epithelial.txt)  )',header=T)_  
    _gwas_pvalue=gwas_file$ALL.FIXED.PVAL_

    2. Obtain the **eQTL p-values** from the eQTL summary statistics file:      _eQTL_file=read.table(file='sampledata/[eqtl_file_gene_MUC4_tissue_Human_nasal_epithelial.txt](https://github.com/FanWang0216/SimpleSum2Colocalization/blob/main/sampledata/eqtl_file_gene_MUC4_tissue_Human_nasal_epithelial.txt)',header=T)_    
_eqtl_pvalue=eqtl_file$pval_
  
    3. Obtain the **covariance matrix** for the GWAS and eQTL summary statistics, respectively:  
  _cov_GWAS=as.matrix(read.table(file='sampledata/[ldGWAS_file_gene_MUC4_tissue_Human_nasal_epithelial.txt](https://github.com/FanWang0216/SimpleSum2Colocalization/blob/main/sampledata/ldGWAS_file_gene_MUC4_tissue_Human_nasal_epithelial.txt)',header=F))_  
  _cov_eQTL=as.matrix(read.table(file='sampledata/[ldeqtl_file_gene_MUC4_tissue_Human_nasal_epithelial.txt](https://github.com/FanWang0216/SimpleSum2Colocalization/blob/main/sampledata/eqtl_file_gene_MUC4_tissue_Human_nasal_epithelial.txt)',header=F))_

    5. Run **SS2 colocalization analysis**:  
    **_stage1(eqtl_pvalue,cov_eQTL)_** produces the p-value for the stage 1 test.  
    **_stage2(eqtl_pvalue,gwas_pvalue,cov_GWAS,use_evidence='log10p')_** produces the p-value for the stage 2 test by using -log10(eQTL p) as eQTL evidence measure. The p-value calculated by using the chi-square statistics as eQTL evidence measure can be done by specifying **_use_evidence='Tsquare'_**.

*   **[SS2_sampledataset_code.R](https://github.com/FanWang0216/SimpleSum2Colocalization/blob/main/SS2_sampledataset_code.R)** implements the above commands by calling functions in SS2_functions.R.

*   **[ld_function.R](https://github.com/FanWang0216/SimpleSum2Colocalization/blob/main/ld_function.R)** provides the R functions for calculating the covariance matrix of Z-scores obtained from **meta-analysis** with either independent or related samples.  

    * R function **_getcov_Zmeta_** produces the covariance of Z-scores from a **meta-analysis** according to the equation:
$$\begin{aligned}
\operatorname{cov}\left(Z_{m e t a, j}, Z_{m e t a, l}\right) &=\operatorname{cov}\left(\frac{\sum_{c=1}^{C} \sqrt{w_{c, j}} Z_{c, j}}{\sqrt{\sum_{c=1}^{C} w_{c, j}}}, \frac{\sum_{c=1}^{C} \sqrt{w_{c, l}} Z_{c, l}}{\sqrt{\sum_{c=1}^{C} w_{c, l}}}\right) \\
&=\frac{\sum_{c=1}^{C} \sqrt{w_{c, j}} \sqrt{w_{c, l}}}{\sqrt{\sum_{c=1}^{C} w_{c, j}} \sqrt{\sum_{c=1}^{C} w_{c, l}}} \operatorname{cov}\left(Z_{c, j}, Z_{c, l}\right);
\end{aligned}$$
 $Z_{meta,j}$: the Z-score for SNP $j$ obtained from a meta-analysis with $C$ sub-studies.  
 $Z_{c,j}$: the Z-score for SNP $j$ from study $c$.   
 $w_{c,j}$: the weight for SNP $j$ from study $c$.   
 $m$: the total number of SNPs at the locus.   
       * **Input variables** for _getcov_Zmeta_: 
          1. **_weightmatrix_**: a $m\times n_{c}$ matrix, where the $j,c$ th element is the weight $w_{c,j}$. An example is provided in the file **[weightmatrix.txt](https://github.com/FanWang0216/SimpleSum2Colocalization/blob/main/sampledata/weightmatrix.txt)**   
          2. **_var_subgroup_**: a list with $C$ elements, where each element is the covariance matrix of Z-scores for a particular study. An example is provided in the file **[13subgroup_ld.Rdata](https://github.com/FanWang0216/SimpleSum2Colocalization/blob/main/sampledata/13subgroup_ld.Rdata)**.   


    * R function **_get_subgroup_var_** calculates the covariance of Z-scores for studies that contain related individuals, according to equation:
$$\begin{aligned}
\operatorname{cov}\left(Z_{c, j}, Z_{c, l}\right)=\frac{G_{j}^{\top} P^{*} G_{l}}{\sqrt{G_{j}^{\top} P^{*} G_{j} \sqrt{G_{l}^{\top} P^{*} G_{l}}}}, \text { with } P^{*}=\Sigma^{-1}-\Sigma^{-1} X\left(X^{\top} \Sigma^{-1} X\right)^{-1} X^{\top} \Sigma^{-1}.
\end{aligned}$$
    * **Input variables** for _get_subgroup_var_:
         1. **_G_**: a $n_{c}\times m$ genotype matrix
	       2. **_X_**: a $n_{c}\times q$ matrix of covariates (i.e. sex and age), including the intercept.  
	       3. **_Sigma_**( $\Sigma$ ): a $n_{c}\times n_{c}$ matrix contains information for sample relatedness. When the Sigma is unknown, the user could obtain the estimated Sigma by using R package _nlme_ or _GMMAT_.
    * Examples for using R function _get_subgroup_var_ and _get_subgroup_var_ are provided in the Rscript **[ld_function.R](https://github.com/FanWang0216/SimpleSum2Colocalization/blob/main/ld_function.R)**.


