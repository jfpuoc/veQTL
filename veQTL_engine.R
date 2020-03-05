## This code performs the Brown-Forsythe Levene-type test for x SNPs for all genes
library(data.table)
library(Biobase)

outfile = "veQTL_analysis" 

## Genotype matrix with the count of the major allele. Snp_id as rownames
snps = read.table("PATH/GenotypeMatrix.txt", header =T)
## Expression data, transcript ids as rownames
expr = read.table("PATH/ExpressionMatrix.txt", header =T, row.names =1)

## SNP and geno must be matrix with samples (in same order) as cols, genes/SNPs as rows and IDs as rownames
results<-apply(snps,1, function(group){

  group <- as.character(group)

  ## exclude uncalled genotype (Denoted with -1) and samples that are in a genotype with a least than 10 samples in total.
  gr<-names(table(group)[table(group) > 10])
  gr <- group %in% gr[gr!=-1]

  group <- group[gr]
  y <- expr[gr]
  #tmp = rownames(group)
  #tmp = c(tmp,rownames(y))

  N <- length(group) # total obs.
  reorder <- order(group)
  group <- group[reorder]
  y <- y[,reorder]
  group <- as.factor(group)

  k <- length(levels(group)) # number of obs
  n <- tapply(group,group, FUN = length) # number of obs. per group

  ## Calc Yi_bar (rowMedians) for each transcript. Equivalant to  Yi_bar <- apply(y, 1, function(x) tapply(x,group, median)) ## Caluclated group medians for every gene --> this is the slowest bit of code
  if (length(levels(group)) == 2){
  Yi_bar1 <- rowMedians(as.matrix(y[group==levels(group)[1]]))
  Yi_bar2 <- rowMedians(as.matrix(y[group==levels(group)[2]]))
  Yi_bar <- rbind(Yi_bar1, Yi_bar2)
  rownames(Yi_bar) <- levels(group)
  colnames(Yi_bar) <- rownames(y)
  } else {
  Yi_bar1 <- rowMedians(as.matrix(y[group==levels(group)[1]]))
  Yi_bar2 <- rowMedians(as.matrix(y[group==levels(group)[2]]))
  Yi_bar3 <- rowMedians(as.matrix(y[group==levels(group)[3]]))
  Yi_bar <- rbind(Yi_bar1, Yi_bar2,Yi_bar3)
  rownames(Yi_bar) <- levels(group)
  colnames(Yi_bar) <- rownames(y)
  }

  Zij <-abs(y - t(apply(Yi_bar,2, function(x) rep(x, n)))) # maybe imporve?

  ## Calc Zi. (rowMeans) for each abs deviation from medians for each transcript. Equivalant to Zi. = apply(Zij,1, function(x) tapply(x,group, mean))
  if (length(levels(group)) == 2){
    Zi.1 <- rowMeans(as.matrix(Zij[group==levels(group)[1]]))
    Zi.2 <- rowMeans(as.matrix(Zij[group==levels(group)[2]]))
    Zi. <- rbind(Zi.1, Zi.2)
    rownames(Zi.) <- levels(group)
    colnames(Zi.) <- rownames(Zij)
  } else {
    Zi.1 <- rowMeans(as.matrix(Zij[group==levels(group)[1]]))
    Zi.2 <- rowMeans(as.matrix(Zij[group==levels(group)[2]]))
    Zi.3 <- rowMeans(as.matrix(Zij[group==levels(group)[3]]))
    Zi. <- rbind(Zi.1, Zi.2,Zi.3)
    rownames(Zi.) <- levels(group)
    colnames(Zi.) <- rownames(y)
  }

  Z.. <- rowMeans(Zij)

  
  # calc test statistics for each transcript -- > W = ((N-k)/(k-1)) * (sum(n*(Zi. - Z..)^2)/sum((Zij-rep(Zi.,n))^2))
  x <- ((N-k)/(k-1))
  x2 <- rep(n, ncol(Zi.))* (Zi. - rep(Z.., each=k))^2
  x3 <-(Zij - t(apply(Zi. , 2, rep, n)))^2 ## maybe improve

  z1 <- colSums(x2)
  z2 <- rowSums(x3)

  W <-  x * z1/z2  
  # p-value, calculate for all test

  W_k2 <- W[k==2]
  k2=2
  p_k2 <- ifelse(W_k2 > 0, (1 - pf(W_k2, k2-1, N-k2)), NA) 
  W_k3 <- W[k ==3]
  k3=3
  p_k3 <- ifelse(W_k3 > 0, (1 - pf(W_k3, k3-1, N-k3)), NA) 
 
  k2_g <- !is.na(p_k2)
  k3_g <- !is.na(p_k3)

  p <- c(p_k2[k2_g], p_k3[k3_g])
  W <- c(W_k2[k2_g], W_k3[k3_g])

  return(c(statistic = W, pvalue = p))
})

gene <- rep(rownames(results[1:(nrow(results)/2),]), each = ncol(results))
SNP <- rep(colnames(results), nrow(results)/2)
gene<-gsub("statistic\\.", "", gene)

stat<-as.vector(t(results[1:(nrow(results)/2),]))
pvalue<-as.vector(t(results[((nrow(results)/2)+1):nrow(results),]))

veQTL_results <- data.frame(ensembl_gene_id=gene, snp_id=SNP, statistic = stat, p.value =pvalue)

save(veQTL_results, file = paste(outfile,".R", sep=""))
