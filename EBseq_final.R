I LIKE THIS ONE
1. remove low abundant genes (quite stringent)
2. MedianNorm normalization (any other will do I guess)
3. SVA on the whole dataset (~conditions only)
4. low cutoff. 0.5 PP
5. Consistent results. DNAJB9 always in top hits, other genes we know in the lab following the expected trajectories.

### EBseqHMM
# statistical analysis in ordered RNA-seq experiments
# identify genes and isoforms that have non-constant expression profile over the time points/positions, and cluster them into expression paths
library("EBSeqHMM")
library("EBSeq")
library("sva")
library("biomaRt")
library("edgeR")
library("sva")
library("RColorBrewer")

## load normalized, sva adjusted counts from edgeR
eset <- read.delim("~/aitor/subread-1.6.2-Linux-x86_64/bin/counts.txt", comment.char="#")
counts_short <- eset[,7:30]
row.names(counts_short) <- eset$Geneid

## metadata
samples <- as.data.frame(cbind(rep("a", 24), rep("b", 24), rep("c", 24)))
samples$V1 <- as.factor(c("D8", "D24", "D24", "D8", "D8", "D24", "M24", "M8", "M8", "M8", "M24", "M24", "T8", "T24", "T24", "T24", "T8", "T8", "TM8", "TM8", "TM24", "TM24", "TM8", "TM24"))
samples$V2 <- as.factor(c("R1", "R3", "R2", "R3", "R2", "R1", "R3", "R2", "R1", "R3", "R2", "R1", "R1", "R2", "R3", "R1", "R2", "R3", "R3", "R1", "R2", "R1", "R2", "R3"))
samples$V3 <- as.factor(c("8", "24", "24", "8", "8", "24", "24", "8", "8", "8", "24", "24", "8", "24", "24", "24", "8", "8", "8", "8", "24", "24", "8", "24"))
colnames(samples) <- c("condition","replicate", "times")
row.names(samples) <- substring(colnames(counts_short), 51, 56)

colnames(counts_short) <- rownames(samples)

## remove low abundant
table(rowSums(counts_short==0)==24)

keep.exprs <- filterByExpr(counts_short)
counts_short <- counts_short[keep.exprs,]

## Normalisation by the method of trimmed mean of M-values (TMM)
nsamples <- ncol(counts_short)
col <- brewer.pal(nsamples, "Paired")

boxplot(log2(counts_short), col=col, main="")
title(main="Log2 Raw data",ylab="Log-counts")

Sizes <- MedianNorm(counts_short)
norm_counts <- GetNormalizedMat(counts_short, Sizes)
boxplot(log2(norm_counts), las=2, col=col, main="")
title(main="MedianNorm Normalised data",ylab="Log2-counts")


plotMDS(log2(norm_counts), top = 1000, labels = NULL, col = as.numeric(samples$condition), cex = 2)
title(main="Treatment groups")
plotMDS(log2(norm_counts), top = 1000, labels = NULL, col = as.numeric(samples$replicate), cex = 2)
title(main="Replicate groups")


## SVA correction for batch effect
mod1 <- model.matrix(~0 + condition, samples)
mod0 <- model.matrix(~1, samples)

svobj <- svaseq(norm_counts, mod1, mod0) 

cleanY = function(y, mod, svs) {
  X = cbind(mod, svs)
  Hat = solve(t(X) %*% X) %*% t(X)
  beta = (Hat %*% t(y))
  rm(Hat)
  gc()
  P = ncol(mod)
  return(y - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),]))
}

cleaned_count <- cleanY(norm_counts, mod1, svobj$sv)

## check distribution of data + PCA
boxplot(log2(cleaned_count), las=2, col=col, main="")
title(main="Normalised+SVA data",ylab="Log2-counts")

plotMDS(log2(cleaned_count), top = 1000, labels = NULL, col = as.numeric(samples$condition), cex = 2)
title(main="Treatment groups")
plotMDS(log2(cleaned_count), top = 1000, labels = NULL, col = as.numeric(samples$replicate), cex = 2)
title(main="Replicate groups")

## Split dataset in the 2 time points wheree we want to find trends (independently)
# SKIP
# Sizes_8h <- Sizes[samples$times == 8]
# Sizes_24h <- Sizes[samples$times == 24]
# if norm + sva applied, set lib size factors to 1
Sizes_8h <- c(1,1,1,1,1,1,1,1,1,1,1,1)
Sizes_24h <- c(1,1,1,1,1,1,1,1,1,1,1,1)

eset_8h <- cleaned_count[,samples$times == 8]
eset_24h <- cleaned_count[,samples$times == 24]

Conditions_8h <- samples$condition[samples$times == 8]
Conditions_8h <- droplevels(Conditions_8h)
Conditions_24h <- samples$condition[samples$times == 24]
Conditions_24h <- droplevels(Conditions_24h)

### Running EBseqHMM. REPEAT FOR 24h
# FCV = , minimum abs(FC) value, default = 2
EBSeqHMMGeneOut_8h <- EBSeqHMMTest(Data=eset_8h, sizeFactors=Sizes_8h, Conditions=Conditions_8h, UpdateRd=5, FCV = 1.3)
GeneConfCalls <- GetConfidentCalls(EBSeqHMMGeneOut_8h, FDR=.05,cutoff=.5, OnlyDynamic=TRUE)
results_8h_TOP <- as.data.frame(GeneConfCalls$Overall) 
table(results_8h_TOP$Most_Likely_Path)

EBSeqHMMGeneOut_24h <- EBSeqHMMTest(Data=eset_24h, sizeFactors=Sizes_24h, Conditions=Conditions_24h, UpdateRd=5, FCV = 1.3)
GeneConfCalls <- GetConfidentCalls(EBSeqHMMGeneOut_24h, FDR=.05,cutoff=.5, OnlyDynamic=TRUE)
results_24h_TOP <- as.data.frame(GeneConfCalls$Overall) 
table(results_24h_TOP$Most_Likely_Path)

## visualize gene of interest
PlotExp(eset_24h, Conditions_24h, Name="ENSG00000128590")
PlotExp(eset_8h, Conditions_8h, Name="ENSG00000128590")

## Anotate results table with GENE SYMBOL, chromosome, etc.
listEnsembl(GRCh=37)
listEnsembl(version=91)
ensembl = useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
todo_genes <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','chromosome_name'), mart=ensembl)

results_8h_TOP <- cbind(todo_genes[match(row.names(results_8h_TOP), todo_genes$ensembl_gene_id),], results_8h_TOP)
results_24h_TOP <- cbind(todo_genes[match(row.names(results_24h_TOP), todo_genes$ensembl_gene_id),], results_24h_TOP)

## save annotated results
write.csv(results_8h_TOP, "results_8h.csv")
write.csv(results_24h_TOP, "results_24h.csv")

