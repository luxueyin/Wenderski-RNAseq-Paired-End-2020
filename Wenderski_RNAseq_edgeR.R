# data analyzed following tutorial from 
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4937821/.

# load R packages
library(dplyr, warn.conflicts = FALSE)
library(tidyr)
library(limma)
library(edgeR)
library(Mus.musculus)
library(RColorBrewer)
library(gplots)


# load files into a vector
files <- c("/Users/xlu/Desktop/Bioinformatics/Wenderski_RNAseq/Outputs/star_countmatrix/SRR11313882.sorted.bamReadsPerGene.out.tab",
           "/Users/xlu/Desktop/Bioinformatics/Wenderski_RNAseq/Outputs/star_countmatrix/SRR11313884.sorted.bamReadsPerGene.out.tab",
           "/Users/xlu/Desktop/Bioinformatics/Wenderski_RNAseq/Outputs/star_countmatrix/SRR11313886.sorted.bamReadsPerGene.out.tab",
           "/Users/xlu/Desktop/Bioinformatics/Wenderski_RNAseq/Outputs/star_countmatrix/SRR11313888.sorted.bamReadsPerGene.out.tab",
           "/Users/xlu/Desktop/Bioinformatics/Wenderski_RNAseq/Outputs/star_countmatrix/SRR11313890.sorted.bamReadsPerGene.out.tab",
           "/Users/xlu/Desktop/Bioinformatics/Wenderski_RNAseq/Outputs/star_countmatrix/SRR11313892.sorted.bamReadsPerGene.out.tab",
           "/Users/xlu/Desktop/Bioinformatics/Wenderski_RNAseq/Outputs/star_countmatrix/SRR11313894.sorted.bamReadsPerGene.out.tab",
           "/Users/xlu/Desktop/Bioinformatics/Wenderski_RNAseq/Outputs/star_countmatrix/SRR11313896.sorted.bamReadsPerGene.out.tab",
           "/Users/xlu/Desktop/Bioinformatics/Wenderski_RNAseq/Outputs/star_countmatrix/SRR11313898.sorted.bamReadsPerGene.out.tab",
           "/Users/xlu/Desktop/Bioinformatics/Wenderski_RNAseq/Outputs/star_countmatrix/SRR11313900.sorted.bamReadsPerGene.out.tab",
           "/Users/xlu/Desktop/Bioinformatics/Wenderski_RNAseq/Outputs/star_countmatrix/SRR11313902.sorted.bamReadsPerGene.out.tab",
           "/Users/xlu/Desktop/Bioinformatics/Wenderski_RNAseq/Outputs/star_countmatrix/SRR11313904.sorted.bamReadsPerGene.out.tab")
read.delim(files[1], nrow=10)

# readDGE function from edgeR allows us to combine all 9 text files in one step.
all_counts <- readDGE(files, path=NULL, columns=c(1,2)) 
all_counts <- all_counts[-c(1:3),]
class(all_counts) 
dim(all_counts)

# organizing sample information
# sample-level information needs to be associated with the columns of the counts 
# matrix information include: experimental variables, both biological and technical, 
# that could have an effect on expression levels.

samplenames <- substring(colnames(all_counts), 77, 87) 
samplenames

colnames(all_counts) <- samplenames
group <- as.factor(c("WT", "WT", "WT", "WT", "WT", 
                     "KO", "KO", "KO", "KO", "KO", "KO", "KO"))
all_counts$samples$group <- group
all_counts$samples

sample_treatment <- as.factor(c("SRR11313882_WT","SRR11313884_WT","SRR11313886_WT","SRR11313888_WT",
                                "SRR11313890_WT","SRR11313892_KO","SRR11313894_KO","SRR11313896_KO",
                                "SRR11313898_KO","SRR11313900_KO","SRR11313902_KO","SRR11313904_KO"))
all_counts$samples$sample_treatment <- sample_treatment
all_counts$samples

# organizing gene annotations
geneid <- rownames(all_counts)

# using Mus.Musculus to generate gene annotation
# gene <- select(Mus.musculus, keys=geneid, columns=c("SYMBOL", "TXCHROM"), 
#               keytype="ENSEMBL")

mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
#head(listAttributes(human),100)
#attributes <- listAttributes(mouse)
#attributes[grep("exon", attributes$name),]

# using biomart to generate gene annotations = more complete
genes <- getBM(attributes = c("ensembl_gene_id","mgi_symbol", "entrezgene_id", 
                              "chromosome_name","start_position","end_position",
                              "description"),
              filter= "ensembl_gene_id",
              values = geneid,
              mart = mouse)

genes$genelength <- abs(genes$end_position - genes$start_position)
#remove duplicates:keeps only first entry
genes <- genes[!duplicated(genes$ensembl_gene_id),]
#match order: df[match(target, df$name),]
genes <- genes[match(geneid,genes$ensembl_gene_id),]

################################################################################
################################################################################
# ONLY RUN ONCE for Mus.Musculus
# gene <- gene[-c(1:3),]
dim(genes)
head(genes)
nrow(genes)


# the data frame from gene annotations is added to the data object and packaged 
# in DGEList-object containing raw count data with associated sample information
# and gene annotations.
all_counts$genes <- genes
all_counts

#########################################################
# Data processing

# transform raw counts to a scale to account for library size difference
# converting raw counts to CPM and log-CPM values using cpm function in edgeR.
# rpkm values can be calculated using rpkm function in edgeR if gene lengths are 
# available.
# cpm value of 1 for a gene = 20 counts in the sample with the lowest sequencing
# depth or 76 counts in the sample with the greatest sequencing depth
cpm <- cpm(all_counts)
keep <- rowSums(cpm(all_counts)>1)>=32
all_counts_cpm <- all_counts[keep,]
dim(all_counts_cpm)
# all_counts_cpm$genes <- genes

lcpm <- cpm(all_counts, log=TRUE)
nrow(lcpm)

# Calculation of the mean and median of library size
L <- mean(all_counts$samples$lib.size) * 1e-6
M <- median(all_counts$samples$lib.size) * 1e-6
c(L,M)
summary(lcpm)
# average library size (L) is 3.7 and median is 3.65

# Plot 1: log-CPM values across 24 samples before and after filtering

# removing lowly expressed genes throughout all samples (~19% genes have 0 counts
# acorss all nine samples).
table(rowSums(all_counts$counts==0)==24)
# 13857 genes have non zero expression levels

# filterByExpr function in edgeR provides an automatic way to filter genes while
# keeping as many genes with worthwhile counts.
keep.exprs <- filterByExpr(all_counts, group=group)
all_counts <- all_counts[keep.exprs,, keep.lib.sizes=FALSE]
dim(all_counts)
# median library size ~3.65 million; 10/3.65 = 2.74. filterByExpr function keeps genes
# with a CPM of at least a 2.74 in at least 10 samples (all cell types have 10 replicates).

lcpm.cutoff <- log2(10/M + 2/L)
# number of genes above the cutoff - 13889
table(rowSums(all_counts$counts==lcpm.cutoff)==24)
nsamples <- ncol(all_counts)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
lcpm <- cpm(all_counts, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")

# Plot 2: comparing unnormalized and normalized data
# normalization factors calculated are automatically stored in x$samples$norm.factors

all_counts_cpm <- calcNormFactors(all_counts, method = "TMM")
all_counts_cpm$samples$norm.factors
# 0.9830853 0.9809906 1.0317399 1.0278990 0.9863788 0.9845002 0.9929199 
# 0.9970182 1.0179439 1.0205731 0.9772660 0.9789304 1.0318855 1.0310936 
# 0.9947554 0.9949212 1.0114429 1.0139570 1.0020946 1.0018065 1.0122644 
# 1.0148627 0.9576999 0.9593897

all_counts_cpm2 <- all_counts_cpm
all_counts_cpm2$samples$norm.factors <- 1
all_counts_cpm2$counts[,1] <- ceiling(all_counts_cpm2$counts[,1]*0.05)
all_counts_cpm2$counts[,2] <- all_counts_cpm2$counts[,2]*5

par(mfrow=c(1,2))
lcpm <- cpm(all_counts_cpm2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="A. Unnormalised data", ylab="Log-cpm")
all_counts_cpm2 <- calcNormFactors(all_counts_cpm2)
all_counts_cpm2$samples$norm.factors

lcpm <- cpm(all_counts_cpm2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="B. Normalised data", ylab="Log-cpm")

# Plot 3: unsupervised clustering of samples
# dimension 1 & 2 - experimental groups
dev.off()
lcpm <- cpm(all_counts_cpm, log=TRUE)
par(mfrow=c(1,1))
col.group <- group
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
plotMDS(lcpm, labels=group, col=col.group)
title(main="A. Sample groups")

# Plot 4: interactive 
# glMDSPlot function generates and html page (launch = TRUE will open it in a browser)
# with a MDS plot on the left and a bar plot showing the proportion of variation 
# explained by each dimension in the right panel.
# hovering on individual points will reveal the sample label.
library(Glimma)
glMDSPlot(lcpm, labels=paste(sample_treatment, sep="_"), 
          groups=all_counts$samples[,c(2)],
          launch=TRUE)
# file:///Users/xlu/Desktop/Bioinformatics/Wenderski_RNAseq/glimma-plots/MDS-Plot.html

# Differential Gene Expression:
# setting up a design matrix for cell population
design <- model.matrix(~0+group)
colnames(design) <- gsub("group", "", colnames(design))
design

# Contrasts for pairwise comparisons between cell populations are set up in 
# limma using the makeContrasts function.
contr.matrix <- makeContrasts(
  WT.vs.KO = KO-WT,
  levels = colnames(design))
contr.matrix

# Plot 5: Mean Variance Relationship
dev.off()
v <- voom(all_counts_cpm2, design, plot=TRUE)
v

# Fitting linear models for comparisons of interest
# Linear modelling in limma is carried out using the lmFit and contrasts.fit functions 
# originally written for application to microarrays.
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)

#summary(decideTests(fit2,adjust.method="fdr", method="global"))
#sum <- summary(decideTests(fit2,adjust.method="BH", method="separate"))
#more weight to fold-changes in the ranking
#treat computes empirical Bayes moderated-t p-values relative to a minimum meaningful fold-change threshold.
#Instead of testing for genes that have true log-fold-changes different from zero, it tests whether the true
#log2-fold-change is greater than lfc in absolute value (McCarthy and Smyth, 2009).
#In other words, it uses an interval null hypothesis, where the interval is [-lfc,lfc].
#When the number of DE genes is large, treat is often useful for giving preference to larger fold-changes and
#for prioritizing genes that are biologically important.
#treat is concerned with p-values rather than posterior odds, so it does not compute the B-statistic lods.
#The idea of thresholding doesn't apply to F-statistics in a straightforward way, so moderated F-statistics are also not computed.
#When lfc=0, treat is identical to eBayes, except that F-statistics and B-statistics are not computed.
#The lfc threshold is usually chosen relatively small, because significantly DE genes must all have fold changes substantially greater than the testing threshold.
#Typical values for lfc are log2(1.1), log2(2) or log2(2). The top genes chosen by treat can be examined using topTreat.
str(efit)
plotSA(efit, main = "Final Model: Mean-variance trend")

summary(decideTests(efit))

# for lfc 0.5 cut off tfit <- treat(vfit, lfc=0.5)
tfit <- treat(vfit)
dt <- decideTests(tfit, adjust.method="BH", method="separate")
sum <- summary(dt)
write.csv(sum,"SummaryCount_DEGs.csv")

###
comparisons=(coef(tfit))
comparisons=colnames(comparisons)
comp_out <- as.data.frame(rownames(v$E))
names(comp_out) <- c("GeneID")
nrowkeep <- nrow(comp_out)
SumTableOut <- NULL
for(i in 1:length(comparisons)){
  #comparison name
  comp=comparisons[i]
  print(comp)
  #make comparisons
  tmp=topTreat(tfit,coef=i,number=nrowkeep,adjust.method="BH")
  #nrow(tmp[(tmp$adj.P.Val<0.05),]) # number of DE genes
  #LogFC values:https://support.bioconductor.org/p/82478/
  tmp$direction <- c("none")
  # cutoff for lfc parameter
  tmp$direction[which(tmp$logFC > 0.5)] = c("Increase")
  tmp$direction[which(tmp$logFC < -0.5)] = c("Decrease")
  tmp$significance <- c("nonDE")
  tmp$significance[which(tmp$adj.P.Val <0.05)] = c("DE")
  #summary counts table based on Ensemble Gene ID counts:
  SumTable <- table(tmp$significance,tmp$direction)
  SumTable <- as.data.frame(SumTable)
  SumTable$comparison <- paste(comp)
  SumTableOut <- rbind(SumTable,SumTableOut)
  #get geneids
  tmp$GeneID <- rownames(tmp)
  #gene gene names and expression levels
  tmp2 <- tmp
  tmp2$comparison <- paste(comp)
  write.csv(tmp2,file = paste(comp,"_DEgenes.csv"))
  #save to output:
  #merge <- merge(comp_out,tmp2, by= "GeneID")
  merge2 <- tmp2 %>% dplyr::select(ensembl_gene_id,logFC,t,P.Value,adj.P.Val,direction,significance)
  colnames(merge2) <- paste(colnames(merge2),comp,sep=".")
  colnames(merge2)[1] <- c("GeneID")
  comp_out <- merge(comp_out, merge2, by="GeneID")
  #data for plot with gene names:
  genenames <- tmp2 %>% dplyr::select(adj.P.Val,logFC,mgi_symbol) %>% distinct()
  #names for plots
  plotname <- gsub("\\."," ",comp)
  plotname <- gsub("vs"," vs ",plotname)
  #volcano plot
  pdf(file = paste(comp,"_Volcano.pdf", sep=""), wi = 9, he = 6, useDingbats=F)
  with(genenames, plot(logFC, -log10(adj.P.Val), pch=20,col="gray", main=paste(plotname,"\nVolcano plot", sep=" "), ylab =c("-log10(adj.pvalue)"),xlab =c("Log Fold Change") ))
  #color points red when sig and log2 FC > 2 and blue if log2 FC < -2
  with(subset(genenames, logFC < -0.5 & -log10(adj.P.Val) > -log10(.05)), points(logFC, -log10(adj.P.Val), pch=20, col="blue"))
  with(subset(genenames, logFC > 0.5 & -log10(adj.P.Val) > -log10(.05)), points(logFC, -log10(adj.P.Val), pch=20, col="red"))
  #add lines
  abline(h = -log10(.05), col = c("black"), lty = 2, lwd = 1)
  abline(v = c(-0.5,0.5), col = "black", lty = 2, lwd = 1)
  #Label points with the textxy function from the calibrate plot
  library(calibrate)
  with(subset(genenames, adj.P.Val<0.05 & abs(logFC)>1.5), textxy(logFC, -log10(adj.P.Val), labs=mgi_symbol, cex=.4))
  dev.off()
}

#get log2CPM counts from voom and put in dataframe:
library(plotrix)
#library(xlsx)
#average log2 CPM and sem
all_counts_cpm2 <- as.data.frame(v$E)
all_counts_cpm2$GeneID<- rownames(v$E)

all_counts_cpm3 <- all_counts_cpm2 %>% group_by(GeneID) %>% gather(sample,log2CPM, 1:(ncol(all_counts_cpm2)-1))
all_counts_cpm3 <- as.data.frame(all_counts_cpm3)

GeneSummary <- all_counts_cpm3 %>% group_by(GeneID, sample) %>%
  summarize(meanlog2CPM = mean(log2CPM)) %>%
  ungroup()  %>%
  #dplyr::select(GeneID, sample, meanlog2CPM) %>%
  spread(sample, meanlog2CPM)

v$genes <- genes
geneinfo <- v$genes
GeneSummary2 <- merge(GeneSummary,geneinfo,by.x="GeneID",by.y="ensembl_gene_id")

write.csv(GeneSummary2,file= "MeanLog2CPM_pergroup_GeneExpression.csv")


# plot 6 ven diagram showing number of genes DE in the comparing WT vs KO
de.common <- which(dt[,1]!=0 & dt[,1]!=0)
length(de.common)

head(tfit$gene$mgi_symbol[de.common], n=20)

vennDiagram(dt[,1:1], circle.col=c( "salmon"))
# write.fit(tfit, dt, file="results.txt")

# temp <- as.data.frame(tfit)

# examining individual DE genes from top to bottom
WT.vs.KO <- topTreat(tfit, coef=1, n=Inf)
head(WT.vs.KO)

plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1], xlim=c(-2,13))
################################################################################

# get row ID that is above lcpm.cutoff threshold
id <- which(apply(lcpm, 1, function(x) any(x > lcpm.cutoff)))
# remove all rows containing values lower than lcpm.cutoff
nlcpm <- lcpm[id,]

# check the number of rows in nlcpm
nrow(nlcpm)
nrow(lcpm)

glMDPlot(tfit, coef=1, status=dt, main=colnames(tfit)[1],
         display.columns = c("mgi_symbol", "ensembl_gene_id", "entrezgene_id"),
         anno = all_counts_cpm2$genes,
         side.main="mgi_symbol", counts=lcpm, groups=group, launch=TRUE)
# file:///Users/xlu/Desktop/Bioinformatics/Wenderski_RNAseq/glimma-plots/MD-Plot.html 
################################################################################

# heatmap
WT.vs.KO.topgenes <- WT.vs.KO
i <- which(v$genes$mgi_symbol %in% WT.vs.KO.topgenes)
mycol <- colorpanel(1000,"blue","white","red")
library(gplots)
heatmap.2(lcpm, scale="row",
          labRow=v$genes$mgi_symbol[i], labCol=group,
          col=mycol, trace="none", density.info="none", 
          margin=c(8,6), lhei=c(2,10), dendrogram="column")
