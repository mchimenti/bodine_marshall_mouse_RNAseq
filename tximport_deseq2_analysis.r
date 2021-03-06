## Analysis of Marshall/Bodine Mouse RNA-Seq data; chow vs HFD, WT vs KO
## Date: 7.25.2018
## Author: Michael Chimenti
## Organism: mm10 / mouse
## Aligners: hisat2 / salmon
## Design: Case/control diet + case/control KO + gender + known batch effects
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

##---------------launch PCA Explorer on dds object 
anno <- get_annotation(ddsTxi, 'mmusculus_gene_ensembl','ensembl_gene_id')
anno <- na.omit(anno)
rldTxi <- rlog(ddsTxi, blind=FALSE)
pcaExplorer(dds=ddsTxi,annotation=anno,rlt=rldTxi)

## look at dispersion estimates 
plotDispEsts(ddsTxi)

## PCA analysis shows possible genotype interaction with sex; add to model and reperform

ddsTxi2 <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ batch + sex + diet + geno + sex:geno)

ddsTxi2 <- ddsTxi2[ rowSums(counts(ddsTxi2)) > 5, ]
ddsTxi2 <- DESeq(ddsTxi2)
resultsNames(ddsTxi2)

plotMA(object = ddsTxi, alpha = 0.05)

#[1] "Intercept"        "batch_2_vs_1"     "batch_3_vs_1"     "sex_M_vs_F"       "diet_HFD_vs_Chow" "geno_WT_vs_KO"   
#[7] "sexM.genoWT"

## Interaction term sig testing
ddsTxi2_lrt <- DESeq(ddsTxi2, test = "LRT", reduced = ~ batch + sex + diet + geno)
res_lrt <- results(ddsTxi2_lrt)

sig_res_int_LRT <-res_lrt %>%
              data.frame() %>%
              rownames_to_column(var="gene") %>% 
              as_tibble() %>% 
              filter(padj < 0.1)

sig_res_int_LRT$ext_gene <- anno[sig_res_int_LRT$gene, "gene_name"]

png("test.png", 1200, 1500, pointsize=20, res=100)
volcanoplot(sig_res_int_LRT, main = "Volcano Plot:", lfcthresh=1.0, sigthresh=0.1, textcx=.6, xlim=c(-10, 11), ylim = c(3,8))
dev.off()

### There are 32 genes with padj < 0.1 between the full model and reduced w/o interaction term 

##############
## DE analysis
##############

## Genotype
res_geno <- results(ddsTxi2, contrast = c("geno","KO","WT"))
res_geno <- na.omit(res_geno)  #drop NA rows
res_geno_sig <- res_geno[res_geno$padj < 0.05 & res_geno$baseMean > 5.0,]
res_geno_ord <- res_geno_sig[order(res_geno_sig$padj),]
res_geno_ord$ext_gene <- anno[row.names(res_geno_ord), "gene_name"]

png("geno_DE_volcano.png", 1200, 1500, pointsize=20, res=100)
volcanoplot(res_geno_ord, main = "Volcano Plot: DE genes across genotypes, KO vs. WT (Male responses)", lfcthresh=1.0, sigthresh=0.1, textcx=.6, xlim=c(-12, 12), ylim = c(3,25))
dev.off()

## Diet
res_diet <- results(ddsTxi2, contrast = c("diet","HFD","Chow"))
res_diet <- na.omit(res_diet)  #drop NA rows
res_diet_sig <- res_diet[res_diet$padj < 0.05 & res_diet$baseMean > 5.0,]
res_diet_ord <- res_diet_sig[order(res_diet_sig$padj),]
res_diet_ord$ext_gene <- anno[row.names(res_diet_ord), "gene_name"]

png("diet_DE_volcano.png", 1200, 1500, pointsize=20, res=100)
volcanoplot(res_diet_ord, main = "Volcano Plot: DE genes across diet, HDF vs. Chow", lfcthresh=0.7, sigthresh=0.1, textcx=.7, xlim=c(-3, 5), ylim = c(3.5,10))
dev.off()

## Sex 
res_sex <- results(ddsTxi2, contrast = c("sex","M","F"))
res_sex <- na.omit(res_sex)  #drop NA rows
res_sex_sig <- res_sex[res_sex$padj < 0.05 & res_sex$baseMean > 5.0,]
res_sex_ord <- res_sex_sig[order(res_sex_sig$padj),]
res_sex_ord$ext_gene <- anno[row.names(res_sex_ord), "gene_name"]

png("sex_DE_volcano.png", 1200, 1500, pointsize=20, res=100)
volcanoplot(res_sex_ord, main = "Volcano Plot: DE genes across sex (M vs. F)", lfcthresh=1.0, sigthresh=0.1, textcx=.6, xlim=c(-10, 11), ylim = c(3,50))
dev.off()

## Genotype interaction with sex 
res_geno_int <- results(ddsTxi2, name = "sexM.genoWT")
res_geno_int <- na.omit(res_geno_int)  #drop NA rows
res_geno_int_sig <- res_geno_int[res_geno_int$padj < 0.05 & res_geno_int$baseMean > 5.0,]
res_geno_int_ord <- res_geno_int_sig[order(res_geno_int_sig$padj),]
res_geno_int_ord$ext_gene <- anno[row.names(res_geno_int_ord), "gene_name"]

png("geno_sex_int.png", 1200, 1500, pointsize=20, res=100)
volcanoplot(res_geno_int_ord, main = "Volcano Plot: DE across geno w/ gender interaction", lfcthresh=0.7, sigthresh=0.2, textcx=.8, xlim=c(-6,6), ylim = c(4,8))
dev.off()

#### DEG plots
degPlot(dds = ddsTxi2, res = res_diet_ord, n = 9, xs = "diet", group = "sex", ann = "ext_gene")
degPlot(dds = ddsTxi2, res = res_geno_ord, n = 9, xs = "geno", group = "sex")
degPlot(dds = ddsTxi2, res = res_sex_ord, n = 9, xs = "sex")
degPlot(ddsTxi2, res_geno_int_ord, n = 9, xs = "geno", group = "sex")

##
res_geno_ashr <- lfcShrink(ddsTxi2, coef = 6, type = 'ashr')
res_diet_ashr <- lfcShrink(ddsTxi2, coef = 5, type = 'ashr')

plotMA(res_geno_ashr)
## write results 

my_cols <- c("baseMean","log2FoldChange","padj","ext_gene")
write.csv(x = res_sex_ord[,my_cols], file = "DE_genes_gender_padj_0p05.csv")
write.csv(x = res_diet_ord[,my_cols], file = "DE_genes_diet_padj_0p05.csv")
write.csv(x = res_geno_ord[,my_cols], file = "DE_genes_geno_padj_0p05.csv")
write.csv(x = res_geno_int_ord[,my_cols], file = "DE_genes_geno_interact_padj_0p05.csv")

##################
## Diet and Sex interaction using factor paste method (rec by DESeq2 Vignette)
##################

samples$sexdiet <- paste0(samples$sex, samples$diet)

ddsTxi3 <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ batch + geno + sexdiet)

ddsTxi3 <- ddsTxi3[ rowSums(counts(ddsTxi3)) > 5, ]
ddsTxi3 <- DESeq(ddsTxi3)

res_fem_hfd <- results(ddsTxi3, contrast = c("sexdiet", "FHFD", "FChow"))
res_fem_hfd <- na.omit(res_fem_hfd)  #drop NA rows
res_fem_hfd_sig <- res_fem_hfd[res_fem_hfd$padj < 0.05 & res_fem_hfd$baseMean > 5.0 & abs(res_fem_hfd$log2FoldChange < 10), ]
res_fem_hfd_ord <- res_fem_hfd_sig[order(res_fem_hfd_sig$padj),]
res_fem_hfd_ord$ext_gene <- anno[row.names(res_fem_hfd_ord), "gene_name"]

png("volcano_female_HFD_v_chow.png", 1200, 1500, pointsize=20, res=100)
volcanoplot(res_fem_hfd_ord, main = "DE genes female HFD vs. Chow, log2FC > 1.5, padj < 0.01", 
            lfcthresh=1.5, sigthresh=0.01, textcx=.5, xlim=c(-5,5), ylim = c(2,25))
dev.off()

res_male_hfd <- results(ddsTxi3, contrast = c("sexdiet", "MHFD", "MChow"))
res_male_hfd <- na.omit(res_male_hfd)  #drop NA rows
res_male_hfd_sig <- res_male_hfd[res_male_hfd$padj < 0.05 & res_male_hfd$baseMean > 5.0 & abs(res_male_hfd$log2FoldChange) < 10, ]
res_male_hfd_ord <- res_male_hfd_sig[order(res_male_hfd_sig$padj),]
res_male_hfd_ord$ext_gene <- anno[row.names(res_male_hfd_ord), "gene_name"]

png("volcano_male_HDF_v_chow.png", 1200, 1500, pointsize=20, res=100)
volcanoplot(res_male_hfd_ord, main = "DE genes male HFD vs. Chow; log2FC > 1.0, padj <0.05", 
            lfcthresh=1.0, sigthresh=0.05, textcx=.5, xlim=c(-4,4), ylim = c(2,15))
dev.off()


write.csv(x = res_fem_hfd_ord[,my_cols], file = 'DE_genes_fem_HFD_vs_Chow_padj_p05.csv')
write.csv(x = res_male_hfd_ord[,my_cols], file = 'DE_genes_male_HFD_vs_Chow_padj_p05.csv')

##
degPlot(ddsTxi3, res = res_fem_hfd, xs = 'sexdiet', n=6)
degPlot(ddsTxi3, res = res_male_hfd, xs = 'sexdiet', n=6)


##################
## Diet and Genotype interaction using factor paste method 
##################

samples$genodiet <- paste0(samples$geno, samples$diet)
ddsTxi4 <- DESeqDataSetFromTximport(txi,
                                    colData = samples,
                                    design = ~ batch + sex + genodiet)

ddsTxi4 <- ddsTxi4[ rowSums(counts(ddsTxi4)) > 5, ]
ddsTxi4 <- DESeq(ddsTxi4)

res_ko_hfd <- results(ddsTxi4, contrast = c("genodiet", "KOHFD", "KOChow"))
res_ko_hfd <- na.omit(res_ko_hfd)  #drop NA rows
res_ko_hfd_sig <- res_ko_hfd[res_ko_hfd$padj < 0.1 & res_ko_hfd$baseMean > 5.0 & abs(res_ko_hfd$log2FoldChange) < 10, ]
res_ko_hfd_ord <- res_ko_hfd_sig[order(res_ko_hfd_sig$padj),]
res_ko_hfd_ord$ext_gene <- anno[row.names(res_ko_hfd_ord), "gene_name"]

png("volcano_ko_HFD_v_chow.png", 1200, 1500, pointsize=20, res=100)
volcanoplot(res_ko_hfd_ord, main = "DE genes Trim62 knockout HFD vs. Chow; log2FC > 0.75 padj <0.1", 
            lfcthresh=0.75, sigthresh=0.1, textcx=.5, xlim=c(-3,3), ylim = c(3,10))
dev.off()

res_wt_hfd <- results(ddsTxi4, contrast = c("genodiet", "WTHFD", "WTChow"))
res_wt_hfd <- na.omit(res_wt_hfd)  #drop NA rows
res_wt_hfd_sig <- res_wt_hfd[res_wt_hfd$padj < 0.1 & res_wt_hfd$baseMean > 5.0 & abs(res_wt_hfd$log2FoldChange) < 10, ]
res_wt_hfd_ord <- res_wt_hfd_sig[order(res_wt_hfd_sig$padj),]
res_wt_hfd_ord$ext_gene <- anno[row.names(res_wt_hfd_ord), "gene_name"]

png("volcano_wt_HFD_v_chow.png", 1200, 1500, pointsize=20, res=100)
volcanoplot(res_wt_hfd_ord, main = "DE genes wild-type Trim62 HFD vs. Chow; log2FC > 0.75, padj <0.1", 
            lfcthresh=0.75, sigthresh=0.1, textcx=.5, xlim=c(-3,3), ylim = c(3,5))
dev.off()

write.csv(x = res_ko_hfd_ord[,my_cols], file = "DE_genes_knockout_HFD_vs_Chow_padj_0p1.csv")
write.csv(x = res_wt_hfd_ord[,my_cols], file = "DE_genes_wildtype_HFD_vs_Chow_padj_0p1.csv")


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

res_wt_hfd <- results(ddsTxiF, contrast = c("genodiet", "WTHFD", "WTChow"))
res_wt_hfd <- na.omit(res_wt_hfd)  #drop NA rows
res_wt_hfd_sig <- res_wt_hfd[res_wt_hfd$padj < 0.1 & res_wt_hfd$baseMean > 5.0 & abs(res_wt_hfd$log2FoldChange) < 10, ]
res_wt_hfd_ord <- res_wt_hfd_sig[order(res_wt_hfd_sig$padj),]
res_wt_hfd_ord$ext_gene <- anno[row.names(res_wt_hfd_ord), "gene_name"]

res_ko_hfd <- results(ddsTxiF, contrast = c("genodiet", "KOHFD", "KOChow"))
res_ko_hfd <- na.omit(res_ko_hfd)  #drop NA rows
res_ko_hfd_sig <- res_ko_hfd[res_ko_hfd$padj < 0.1 & res_ko_hfd$baseMean > 5.0 & abs(res_ko_hfd$log2FoldChange) < 10, ]
res_ko_hfd_ord <- res_ko_hfd_sig[order(res_ko_hfd_sig$padj),]
res_ko_hfd_ord$ext_gene <- anno[row.names(res_ko_hfd_ord), "gene_name"]

