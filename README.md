# veQTL
Variable expression quantitative trait loci is a tool to detect variants (genotypes) that are associated with gene expression variability.
To run, genotypes need to be converted to a 012 matrix where genotypes are represented as counts of major (or minor) allele. Use -1 for no calls. 

The veQTL_engine.R script is setup to run on expression and genotype matrices with identical sample orders. The script will compute the robust Brown-Forsythe test for equal variance. The output is a dataframe with genes as rownames from expression matrix and snps as rownames from snp matrix with the test statistic and p value (unadjusted).
