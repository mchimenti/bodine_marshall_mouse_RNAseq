## Analysis of Marshall/Bodine Mouse RNA-Seq data; chow vs HFD, WT vs KO
## Date: 8.30.2018
## Author: Michael Chimenti
## Organism: mm10 / mouse
## Aligners: hisat2 / salmon
## Design: Diet + Genotype interactions by gender 
## Reps: 4

##########
## Imports
##########

source("https://bioconductor.org/biocLite.R")
biocLite("DEGreport")

#negative binomial GLM and related
library('DESeq2')
library('calibrate')
library('tximport')
library('readr')
#annotation
library('biomaRt')
library("AnnotationDbi")
library("org.Hs.eg.db")
#Exploratory analysis
library('tidyverse')
library('pcaExplorer')
#pathway and gene clusters
library('DEGreport')
#library(pathview)
#library(gage)
#library(gageData)
#library(ggplot2)

setwd("~/iihg/RNA_seq/bodine_lab/project_marshall_aug2018/") 

###########
##Function Defs
###########

get_annotation <- function(dds, biomart_dataset, idtype){
  if(is.null(biomart_dataset))
    stop("Select a species to generate the corresponding annotation.
         To obtain a list, type mart = useMart('ensembl'), followed by listDatasets(mart).")
  
  mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                  host="www.ensembl.org",
                  dataset=biomart_dataset)
  anns <- getBM(attributes = c(idtype, "external_gene_name", "description"),
                filters = idtype,
                values = rownames(dds),
                mart = mart)
  
  # keep and match with the ones that are actually there
  anns2 <- anns[match(rownames(dds), anns[, 1]), ]
  rownames(anns2) <- rownames(dds)
  # rename the columns rsp. add row names to be consistent with other function
  colnames(anns2) <- c("gene_id","gene_name","description")
  
  return(anns2)
}

## Volcano Plot function 
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="topright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=ext_gene, cex=textcx, offset=0.3, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}

#######################################
## tximport > DESeq2 > PCAExplorer
#######################################


samples <- read.table("samples.csv", sep=',', header=TRUE)
rownames(samples) <- samples$sample
samples$batch <- as.factor(samples$batch)

files <- file.path(getwd(), samples$sample, 'salmon', 'quant.sf')
names(files) <- samples$sample

tx2gene <- read_csv(file.path(getwd(), "tx2gene.csv"), col_names = FALSE)

#tx2gene$X1 <- tx2gene$X1 %>%
#  strsplit(split = '.', fixed = TRUE) %>%
#  sapply( "[", 1)  ## obtuse sapply statement needed b/c of annoying way strsplit returns list of lists

txi <- tximport(files, type="salmon", tx2gene=tx2gene)

ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ batch + sex + diet + geno)

ddsTxi <- ddsTxi[ rowSums(counts(ddsTxi)) > 5, ]
ddsTxi <- DESeq(ddsTxi)


#################
## Splitting out the males and females: Andrea wants everything contrasted within males only and females only 
## diet and geno interactions by male/female 
#################

#samples$genodiet <- paste0(samples$geno, samples$diet)
#ddsTxi4 <- DESeqDataSetFromTximport(txi,
#                                    colData = samples,
#                                    design = ~ batch + sex + genodiet)

ddsTxiF <- ddsTxi4[ , ddsTxi4$sex == 'F']
ddsTxiF$sex <- droplevels(ddsTxiF$sex)
design(ddsTxiF) <- ~batch + genodiet
ddsTxiF <- ddsTxiF[ rowSums(counts(ddsTxiF)) > 5, ]
ddsTxiF <- DESeq(ddsTxiF)

ddsTxiM <- ddsTxi4[ , ddsTxi4$sex == 'M']
ddsTxiM$sex <- droplevels(ddsTxiM$sex)
design(ddsTxiM) <- ~batch + genodiet
ddsTxiM <- ddsTxiM[ rowSums(counts(ddsTxiM)) > 5, ]
ddsTxiM <- DESeq(ddsTxiM)


##Females

# Wild-types High-Fat vs. Chow
res_wt_hfd <- results(ddsTxiF, contrast = c("genodiet", "WTHFD", "WTChow"))
res_wt_hfd <- na.omit(res_wt_hfd)  #drop NA rows
res_wt_hfd_sig <- res_wt_hfd[res_wt_hfd$padj < 0.1 & res_wt_hfd$baseMean > 5.0 & abs(res_wt_hfd$log2FoldChange) < 10, ]
res_wt_hfd_ord <- res_wt_hfd_sig[order(res_wt_hfd_sig$padj),]
res_wt_hfd_ord$ext_gene <- anno[row.names(res_wt_hfd_ord), "gene_name"]

png("volcano_ko_HFD_v_chow.png", 1200, 1500, pointsize=20, res=100)
volcanoplot(res_ko_hfd_ord, main = "DE genes Trim62 knockout HFD vs. Chow; log2FC > 0.75 padj <0.1", 
            lfcthresh=0.75, sigthresh=0.1, textcx=.5, xlim=c(-3,3), ylim = c(3,10))
dev.off()
write.csv(x = res_ko_hfd_ord[,my_cols], file = "DE_genes_knockout_HFD_vs_Chow_padj_0p1.csv")

# Knock outs (Trim63 -/-) high-fat vs. chow
res_ko_hfd <- results(ddsTxiF, contrast = c("genodiet", "KOHFD", "KOChow"))
res_ko_hfd <- na.omit(res_ko_hfd)  #drop NA rows
res_ko_hfd_sig <- res_ko_hfd[res_ko_hfd$padj < 0.1 & res_ko_hfd$baseMean > 5.0 & abs(res_ko_hfd$log2FoldChange) < 10, ]
res_ko_hfd_ord <- res_ko_hfd_sig[order(res_ko_hfd_sig$padj),]
res_ko_hfd_ord$ext_gene <- anno[row.names(res_ko_hfd_ord), "gene_name"]

png("volcano_ko_HFD_v_chow.png", 1200, 1500, pointsize=20, res=100)
volcanoplot(res_ko_hfd_ord, main = "DE genes Trim62 knockout HFD vs. Chow; log2FC > 0.75 padj <0.1", 
            lfcthresh=0.75, sigthresh=0.1, textcx=.5, xlim=c(-3,3), ylim = c(3,10))
dev.off()
write.csv(x = res_ko_hfd_ord[,my_cols], file = "DE_genes_knockout_HFD_vs_Chow_padj_0p1.csv")

# Chow baseline, Knockout vs. Chow
res_ko_wt <- results(ddsTxiF, contrast = c("genodiet", "KOChow", "WTChow"))

res_ko_wt <- na.omit(res_ko_wt)  #drop NA rows
res_ko_wt_sig <- res_ko_wt[res_ko_wt$padj < 0.1 & res_ko_wt$baseMean > 5.0 & abs(res_ko_wt$log2FoldChange) < 10, ]
res_ko_wt_ord <- res_ko_wt_sig[order(res_ko_wt_sig$padj),]
res_ko_wt_ord$ext_gene <- anno[row.names(res_ko_wt_ord), "gene_name"]

png("volcano_ko_HFD_v_chow.png", 1200, 1500, pointsize=20, res=100)
volcanoplot(res_ko_hfd_ord, main = "DE genes Trim62 knockout HFD vs. Chow; log2FC > 0.75 padj <0.1", 
            lfcthresh=0.75, sigthresh=0.1, textcx=.5, xlim=c(-3,3), ylim = c(3,10))
dev.off()
write.csv(x = res_ko_hfd_ord[,my_cols], file = "DE_genes_knockout_HFD_vs_Chow_padj_0p1.csv")

##Males only 

# Wild-types High-Fat vs. Chow
res_wt_hfd <- results(ddsTxiM, contrast = c("genodiet", "WTHFD", "WTChow"))
res_wt_hfd <- na.omit(res_wt_hfd)  #drop NA rows
res_wt_hfd_sig <- res_wt_hfd[res_wt_hfd$padj < 0.1 & res_wt_hfd$baseMean > 5.0 & abs(res_wt_hfd$log2FoldChange) < 10, ]
res_wt_hfd_ord <- res_wt_hfd_sig[order(res_wt_hfd_sig$padj),]
res_wt_hfd_ord$ext_gene <- anno[row.names(res_wt_hfd_ord), "gene_name"]

png("volcano_ko_HFD_v_chow.png", 1200, 1500, pointsize=20, res=100)
volcanoplot(res_ko_hfd_ord, main = "DE genes Trim62 knockout HFD vs. Chow; log2FC > 0.75 padj <0.1", 
            lfcthresh=0.75, sigthresh=0.1, textcx=.5, xlim=c(-3,3), ylim = c(3,10))
dev.off()
write.csv(x = res_ko_hfd_ord[,my_cols], file = "DE_genes_knockout_HFD_vs_Chow_padj_0p1.csv")

# Knock outs (Trim63 -/-) high-fat vs. chow
res_ko_hfd <- results(ddsTxiM, contrast = c("genodiet", "KOHFD", "KOChow"))
res_ko_hfd <- na.omit(res_ko_hfd)  #drop NA rows
res_ko_hfd_sig <- res_ko_hfd[res_ko_hfd$padj < 0.1 & res_ko_hfd$baseMean > 5.0 & abs(res_ko_hfd$log2FoldChange) < 10, ]
res_ko_hfd_ord <- res_ko_hfd_sig[order(res_ko_hfd_sig$padj),]
res_ko_hfd_ord$ext_gene <- anno[row.names(res_ko_hfd_ord), "gene_name"]

png("volcano_ko_HFD_v_chow.png", 1200, 1500, pointsize=20, res=100)
volcanoplot(res_ko_hfd_ord, main = "DE genes Trim62 knockout HFD vs. Chow; log2FC > 0.75 padj <0.1", 
            lfcthresh=0.75, sigthresh=0.1, textcx=.5, xlim=c(-3,3), ylim = c(3,10))
dev.off()
write.csv(x = res_ko_hfd_ord[,my_cols], file = "DE_genes_knockout_HFD_vs_Chow_padj_0p1.csv")

# Chow baseline, Knockout vs. Chow
res_ko_wt <- results(ddsTxiM, contrast = c("genodiet", "KOChow", "WTChow"))

res_ko_wt <- na.omit(res_ko_wt)  #drop NA rows
res_ko_wt_sig <- res_ko_wt[res_ko_wt$padj < 0.1 & res_ko_wt$baseMean > 5.0 & abs(res_ko_wt$log2FoldChange) < 10, ]
res_ko_wt_ord <- res_ko_wt_sig[order(res_ko_wt_sig$padj),]
res_ko_wt_ord$ext_gene <- anno[row.names(res_ko_wt_ord), "gene_name"]

png("volcano_ko_HFD_v_chow.png", 1200, 1500, pointsize=20, res=100)
volcanoplot(res_ko_hfd_ord, main = "DE genes Trim62 knockout HFD vs. Chow; log2FC > 0.75 padj <0.1", 
            lfcthresh=0.75, sigthresh=0.1, textcx=.5, xlim=c(-3,3), ylim = c(3,10))
dev.off()
write.csv(x = res_ko_hfd_ord[,my_cols], file = "DE_genes_knockout_HFD_vs_Chow_padj_0p1.csv")


