#Simple Sum 2 Colocalization Analysis
-   The **[sampledata](https://github.com/FanWang0216/SimpleSum2Colocalization/tree/main/sampledata)** folder provides all the sample datasets for conducting the Simple Sum 2 colocalization analysis, including the GWAS and eQTL files for two genes: *MUC20* and *MUC4* in tissue human nasal epithelial.

-   **[SS2\_functions.R](https://github.com/FanWang0216/SimpleSum2Colocalization/blob/main/SS2_functions.R)** contains all the R functions for conducting the SS2 test which has been implemented in the LocusFocus, with an additional R option for using the chi-square statsitics as eQTL evidence measure.

-   **Commands** to run Simple Sum 2 for one GWAS/eQTL combination:

    1.  Obtain the **GWAS p-values** from the GWAS summary statistics file: *gwas\_file=read.table(file='sampledata/[GWAS\_file\_gene\_MUC4\_tissue\_Human\_nasal\_epithelial.txt](https://github.com/FanWang0216/SimpleSum2Colocalization/blob/main/sampledata/GWAS_file_gene_MUC4_tissue_Human_nasal_epithelial.txt) )',header=T)*
        *gwas\_pvalue=gwas\_file$ALL.FIXED.PVAL*

    2.  Obtain the **eQTL p-values** from the eQTL summary statistics file: *eQTL\_file=read.table(file='sampledata/[eqtl\_file\_gene\_MUC4\_tissue\_Human\_nasal\_epithelial.txt](https://github.com/FanWang0216/SimpleSum2Colocalization/blob/main/sampledata/eqtl_file_gene_MUC4_tissue_Human_nasal_epithelial.txt)',header=T)*
        *eqtl\_pvalue=eqtl\_file$pval*

    3.  Obtain the **covariance matrix** for the GWAS and eQTL summary statistics, respectively:
        *cov\_GWAS=as.matrix(read.table(file='sampledata/[ldGWAS\_file\_gene\_MUC4\_tissue\_Human\_nasal\_epithelial.txt](https://github.com/FanWang0216/SimpleSum2Colocalization/blob/main/sampledata/ldGWAS_file_gene_MUC4_tissue_Human_nasal_epithelial.txt)',header=F))*
        *cov\_eQTL=as.matrix(read.table(file='sampledata/[ldeqtl\_file\_gene\_MUC4\_tissue\_Human\_nasal\_epithelial.txt](https://github.com/FanWang0216/SimpleSum2Colocalization/blob/main/sampledata/eqtl_file_gene_MUC4_tissue_Human_nasal_epithelial.txt)',header=F))*

    4.  Run **SS2 colocalization analysis**:
        ***stage1(eqtl\_pvalue,cov\_eQTL)*** produces the p-value for the stage 1 test.
        ***stage2(eqtl\_pvalue,gwas\_pvalue,cov\_GWAS,use\_evidence='log10p')*** produces the p-value for the stage 2 test by using -log10(eQTL p) as eQTL evidence measure. The p-value calculated by using the chi-square statistics as eQTL evidence measure can be done by specifying ***use\_evidence='Tsquare'***.

-   **[SS2\_sampledataset\_code.R](https://github.com/FanWang0216/SimpleSum2Colocalization/blob/main/SS2_sampledataset_code.R)** implements the above commands by calling functions in SS2\_functions.R.

-   **[ld\_function.R](https://github.com/FanWang0216/SimpleSum2Colocalization/blob/main/ld_function.R)** provides the R functions for calculating the covariance matrix of Z-scores obtained from **meta-analysis** with either independent or related samples.

    -   R function ***getcov\_Zmeta*** produces the covariance of Z-scores from a **meta-analysis** according to the equation:

        ![\\begin{aligned}
        \\operatorname{cov}\\left(Z\_{m e t a, j}, Z\_{m e t a, l}\\right) &=\\operatorname{cov}\\left(\\frac{\\sum\_{c=1}^{C} \\sqrt{w\_{c, j}} Z\_{c, j}}{\\sqrt{\\sum\_{c=1}^{C} w\_{c, j}}}, \\frac{\\sum\_{c=1}^{C} \\sqrt{w\_{c, l}} Z\_{c, l}}{\\sqrt{\\sum\_{c=1}^{C} w\_{c, l}}}\\right) \\\\
        &=\\frac{\\sum\_{c=1}^{C} \\sqrt{w\_{c, j}} \\sqrt{w\_{c, l}}}{\\sqrt{\\sum\_{c=1}^{C} w\_{c, j}} \\sqrt{\\sum\_{c=1}^{C} w\_{c, l}}} \\operatorname{cov}\\left(Z\_{c, j}, Z\_{c, l}\\right);
        \\end{aligned}](https://latex.codecogs.com/png.latex?%5Cbegin%7Baligned%7D%0A%5Coperatorname%7Bcov%7D%5Cleft%28Z_%7Bm%20e%20t%20a%2C%20j%7D%2C%20Z_%7Bm%20e%20t%20a%2C%20l%7D%5Cright%29%20%26%3D%5Coperatorname%7Bcov%7D%5Cleft%28%5Cfrac%7B%5Csum_%7Bc%3D1%7D%5E%7BC%7D%20%5Csqrt%7Bw_%7Bc%2C%20j%7D%7D%20Z_%7Bc%2C%20j%7D%7D%7B%5Csqrt%7B%5Csum_%7Bc%3D1%7D%5E%7BC%7D%20w_%7Bc%2C%20j%7D%7D%7D%2C%20%5Cfrac%7B%5Csum_%7Bc%3D1%7D%5E%7BC%7D%20%5Csqrt%7Bw_%7Bc%2C%20l%7D%7D%20Z_%7Bc%2C%20l%7D%7D%7B%5Csqrt%7B%5Csum_%7Bc%3D1%7D%5E%7BC%7D%20w_%7Bc%2C%20l%7D%7D%7D%5Cright%29%20%5C%5C%0A%26%3D%5Cfrac%7B%5Csum_%7Bc%3D1%7D%5E%7BC%7D%20%5Csqrt%7Bw_%7Bc%2C%20j%7D%7D%20%5Csqrt%7Bw_%7Bc%2C%20l%7D%7D%7D%7B%5Csqrt%7B%5Csum_%7Bc%3D1%7D%5E%7BC%7D%20w_%7Bc%2C%20j%7D%7D%20%5Csqrt%7B%5Csum_%7Bc%3D1%7D%5E%7BC%7D%20w_%7Bc%2C%20l%7D%7D%7D%20%5Coperatorname%7Bcov%7D%5Cleft%28Z_%7Bc%2C%20j%7D%2C%20Z_%7Bc%2C%20l%7D%5Cright%29%3B%0A%5Cend%7Baligned%7D "\begin{aligned}
        \operatorname{cov}\left(Z_{m e t a, j}, Z_{m e t a, l}\right) &=\operatorname{cov}\left(\frac{\sum_{c=1}^{C} \sqrt{w_{c, j}} Z_{c, j}}{\sqrt{\sum_{c=1}^{C} w_{c, j}}}, \frac{\sum_{c=1}^{C} \sqrt{w_{c, l}} Z_{c, l}}{\sqrt{\sum_{c=1}^{C} w_{c, l}}}\right) \\
        &=\frac{\sum_{c=1}^{C} \sqrt{w_{c, j}} \sqrt{w_{c, l}}}{\sqrt{\sum_{c=1}^{C} w_{c, j}} \sqrt{\sum_{c=1}^{C} w_{c, l}}} \operatorname{cov}\left(Z_{c, j}, Z_{c, l}\right);
        \end{aligned}")

         ![Z\_{meta,j}](https://latex.codecogs.com/png.latex?Z_%7Bmeta%2Cj%7D "Z_{meta,j}"): the Z-score for SNP ![j](https://latex.codecogs.com/png.latex?j "j") obtained from a meta-analysis with ![C](https://latex.codecogs.com/png.latex?C "C") sub-studies.  
        ![Z\_{c,j}](https://latex.codecogs.com/png.latex?Z_%7Bc%2Cj%7D "Z_{c,j}"): the Z-score for SNP ![j](https://latex.codecogs.com/png.latex?j "j") from study ![c](https://latex.codecogs.com/png.latex?c "c").  
        ![w\_{c,j}](https://latex.codecogs.com/png.latex?w_%7Bc%2Cj%7D "w_{c,j}"): the weight for SNP ![j](https://latex.codecogs.com/png.latex?j "j") from study ![c](https://latex.codecogs.com/png.latex?c "c").  
        ![m](https://latex.codecogs.com/png.latex?m "m"): the total number of SNPs at the locus.  
    -   **Input variables** for *getcov\_Zmeta*:
        1.  ***weightmatrix***: a ![m\\times n\_{c}](https://latex.codecogs.com/png.latex?m%5Ctimes%20n_%7Bc%7D "m\times n_{c}") matrix, where the ![j,c](https://latex.codecogs.com/png.latex?j%2Cc "j,c") th element is the weight ![w\_{c,j}](https://latex.codecogs.com/png.latex?w_%7Bc%2Cj%7D "w_{c,j}"). An example is provided in the file **[weightmatrix.txt](https://github.com/FanWang0216/SimpleSum2Colocalization/blob/main/sampledata/weightmatrix.txt)**
        2.  ***var\_subgroup***: a list with ![C](https://latex.codecogs.com/png.latex?C "C") elements, where each element is the covariance matrix of Z-scores for a particular study. An example is provided in the file **[13subgroup\_ld.Rdata](https://github.com/FanWang0216/SimpleSum2Colocalization/blob/main/sampledata/13subgroup_ld.Rdata)**.
    -   R function ***get\_subgroup\_var*** calculates the covariance of Z-scores for studies that contain related individuals, according to equation:

        ![\\begin{aligned}
        \\operatorname{cov}\\left(Z\_{c, j}, Z\_{c, l}\\right)=\\frac{G\_{j}^{\\top} P^{\*} G\_{l}}{\\sqrt{G\_{j}^{\\top} P^{\*} G\_{j} \\sqrt{G\_{l}^{\\top} P^{\*} G\_{l}}}}, \\text { with } P^{\*}=\\Sigma^{-1}-\\Sigma^{-1} X\\left(X^{\\top} \\Sigma^{-1} X\\right)^{-1} X^{\\top} \\Sigma^{-1}.
        \\end{aligned}](https://latex.codecogs.com/png.latex?%5Cbegin%7Baligned%7D%0A%5Coperatorname%7Bcov%7D%5Cleft%28Z_%7Bc%2C%20j%7D%2C%20Z_%7Bc%2C%20l%7D%5Cright%29%3D%5Cfrac%7BG_%7Bj%7D%5E%7B%5Ctop%7D%20P%5E%7B%2A%7D%20G_%7Bl%7D%7D%7B%5Csqrt%7BG_%7Bj%7D%5E%7B%5Ctop%7D%20P%5E%7B%2A%7D%20G_%7Bj%7D%20%5Csqrt%7BG_%7Bl%7D%5E%7B%5Ctop%7D%20P%5E%7B%2A%7D%20G_%7Bl%7D%7D%7D%7D%2C%20%5Ctext%20%7B%20with%20%7D%20P%5E%7B%2A%7D%3D%5CSigma%5E%7B-1%7D-%5CSigma%5E%7B-1%7D%20X%5Cleft%28X%5E%7B%5Ctop%7D%20%5CSigma%5E%7B-1%7D%20X%5Cright%29%5E%7B-1%7D%20X%5E%7B%5Ctop%7D%20%5CSigma%5E%7B-1%7D.%0A%5Cend%7Baligned%7D "\begin{aligned}
        \operatorname{cov}\left(Z_{c, j}, Z_{c, l}\right)=\frac{G_{j}^{\top} P^{*} G_{l}}{\sqrt{G_{j}^{\top} P^{*} G_{j} \sqrt{G_{l}^{\top} P^{*} G_{l}}}}, \text { with } P^{*}=\Sigma^{-1}-\Sigma^{-1} X\left(X^{\top} \Sigma^{-1} X\right)^{-1} X^{\top} \Sigma^{-1}.
        \end{aligned}")

    -   **Input variables** for *get\_subgroup\_var*:
        1.  ***G***: a ![n\_{c}\\times m](https://latex.codecogs.com/png.latex?n_%7Bc%7D%5Ctimes%20m "n_{c}\times m") genotype matrix
        2.  ***X***: a ![n\_{c}\\times q](https://latex.codecogs.com/png.latex?n_%7Bc%7D%5Ctimes%20q "n_{c}\times q") matrix of covariates (i.e. sex and age), including the intercept.
        3.  ***Sigma***( ![\\Sigma](https://latex.codecogs.com/png.latex?%5CSigma "\Sigma") ): a ![n\_{c}\\times n\_{c}](https://latex.codecogs.com/png.latex?n_%7Bc%7D%5Ctimes%20n_%7Bc%7D "n_{c}\times n_{c}") matrix contains information for sample relatedness. When the Sigma is unknown, the user could obtain the estimated Sigma by using R package *nlme* or *GMMAT*.
    -   Examples for using R function *get\_subgroup\_var* and *get\_subgroup\_var* are provided in the Rscript **[ld\_function.R](https://github.com/FanWang0216/SimpleSum2Colocalization/blob/main/ld_function.R)**.
