set.seed("1234")
library(DESeq2)
library(dplyr)
#library(vidger)
library("BiocParallel")
args = commandArgs(trailingOnly=TRUE)
register(MulticoreParam(4))
countdata = read.table(args[1], header = T, row.names = NULL, check.names = F)
colnames(countdata)[1] = "Geneid";
countdata = aggregate(.~Geneid, countdata, max)
row.names(countdata) = countdata[,1]; countdata = countdata[,-1]
sampleInfo <- read.table(args[2], header = T,sep = "\t", row.names = 1)
countdata = countdata[,as.character(as.matrix(row.names(sampleInfo)))]
countdata  =  countdata[rowSums(countdata) >=  10,]
countdata = data.matrix(countdata)
ddsMat <- DESeqDataSetFromMatrix(countData = countdata, colData = sampleInfo, design = ~Replicate)
ddsMat<- estimateSizeFactors(ddsMat)
rld <- rlog(ddsMat, blind = FALSE)
system("mkdir -p data")
write.table(assay(rld), file = "data/readCountRlogNorm.xls", sep = "\t",col.names = NA, quote=FALSE)
diffExp <- DESeq(ddsMat)
system("mkdir -p figures")
pdf(file = "figures/Dispersion_plots.pdf", height = 5, width = 5.5)
plotDispEsts(diffExp)
dev.off()
#####Extra information#####
## counts per sample
sink("data/size_factors.txt")
total_counts = apply(assay(diffExp), 2, sum)
sizeFactors(diffExp)
sink()

###############Samples-Relationships##########
library(RColorBrewer)
library(gplots)
library(circlize)
library(ComplexHeatmap)
library(Hmisc)
#distance matrix

col4 <- colorRampPalette(c("darkblue","darkgreen", "yellow", "darkviolet", "darkmagenta"))
sampleDists <- as.matrix(dist( t(assay(rld)) ))
write.table(sampleDists,file = "data/rlog_Normalized_Values.txt",sep = "\t", quote=FALSE)
ht = Heatmap(sampleDists, col = col4(100), heatmap_legend_param  =  list(title = NULL))
pdf(file = "figures/SampleDistanceRlog.pdf", height = 5, width = 5.5)
draw(ht,heatmap_legend_side  =  "left")
dev.off()

#cor matrix
ht = Heatmap(cor(assay(rld)), col = col4(100), heatmap_legend_param  =  list(title = NULL))
pdf(file = "figures/SampleCorrelationRlog.pdf", height = 5, width = 5.5)
draw(ht,heatmap_legend_side  =  "left")
dev.off()

#PCA
pdf("figures/PCA_var1_var2.pdf", width=7, height=6)
data <- plotPCA(rld, intgroup="Replicate", returnData=TRUE)
percentVar = round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=Replicate, shape=name)) +
  geom_hline(aes(yintercept=0), colour="grey") +
  geom_vline(aes(xintercept=0), colour="grey") +
  geom_point(size=5)+
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw(base_size = 14) +
  ggtitle("PCA\n") + labs(color="Groups", shape="Sample Names")+
  scale_shape_manual(values=c(0:18,33:17))
dev.off()

################ Contrasts to be performed ##########
normReadCount = counts(diffExp, normalized = TRUE)
write.table(normReadCount, file = "data/readCountNorm.xls", sep = "\t",col.names = NA, quote=FALSE)
mCountdata =  data.frame(factor(sampleInfo$Replicate),t(normReadCount), check.names  =  FALSE)
colnames(mCountdata)[1] = "Sample"
mCountdata = aggregate(.~Sample, mCountdata, mean)
row.names(mCountdata) = mCountdata[,1];
mCountdata  =  mCountdata[,-1];mCountdata = t(mCountdata)
ind = match(colnames(mCountdata),unique(sampleInfo$Replicate))
mCountdata = mCountdata[,ind]
write.table(mCountdata, file = "data/readCountMatrixMergedRepli.xls", sep = "\t",col.names = NA, quote=FALSE)
mCountData=mCountdata
comp = read.table(args[3], header = T, row.names = NULL)
outs<-within(comp, contrasts <- paste(comp$Contrast,comp$Reference,sep='_vs_'))[3]
outs
comp = data.frame("Replicate",comp[1:(length(comp)/2)],comp[((length(comp)/2)+1):length(comp)])
listDiffMat_fdr = apply(comp, 1, function(x){  results(diffExp, contrast = x, alpha=0.05)})

####### Results, Volcano and MA ##############

#sapply(1:length(listDiffMat_fdr),function(x)write.table(listDiffMat_fdr[[x]], file=sprintf("data/%s%d.All.xls", outs$contrasts[x], x), quote=FALSE, sep="\t"))
for (i in 1:length(comp$Reference)) {
	#output_comps<- paste("data/",comp$Contrast[i],"_vs_",comp$Reference[i],".all.xls", sep="")
	output_comps<- paste("data/",outs$contrasts[i],".all.xls", sep="")
	mat<- as.data.frame(listDiffMat_fdr[i])
	mat$diffexpressed <- "No_change"
	mat$diffexpressed[mat$log2FoldChange >= 1 & mat$pvalue < 0.05] <- "Up"
	mat$diffexpressed[mat$log2FoldChange <= -1 & mat$pvalue < 0.05] <- "Down"
	mat$diffexpressed[mat$log2FoldChange >= 1 & mat$pvalue >= 0.05] <- "Non_significant"
	mat$diffexpressed[mat$log2FoldChange <= -1 & mat$pvalue >= 0.05] <- "Non_significant"
    write.table(file=output_comps, mat, sep="\t", quote=FALSE)
	
	#output_filt<- paste("data/",comp$Contrast[i],"_vs_",comp$Reference[i],".DEGs.xls", sep="")
	output_filt<- paste("data/",outs$contrasts[i],".DEGs.xls", sep="")
	filter_mat<- mat %>% filter(padj < 0.05) %>% filter(log2FoldChange >= 1 | log2FoldChange <= -1 )
    write.table(file=output_filt, filter_mat, sep="\t", quote=FALSE)
	
    mycolors <- c("blue", "red", "green", "black")
    names(mycolors) <- c("Up", "Down", "Non_significant", "No_change")
    #output_volcano<- paste("figures/", comp$Contrast[i],"_vs_",comp$Reference[i],".Volcano.pdf", sep="")
    output_volcano<- paste("figures/", outs$contrasts[i],".Volcano.pdf", sep="")
    pdf(file=output_volcano, width=7, height=6)
    p<- ggplot(data=mat, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed)) + geom_point(size=0.5) + theme_minimal() +  geom_vline(xintercept=c(-1, 1), col="red") + geom_hline(yintercept=-log10(0.05), col="red") + scale_colour_manual(values = mycolors)
    print(p)
    dev.off()
    
    #output_maplot<- paste("figures/", comp$Contrast[i],"_vs_",comp$Reference[i],".MA.pdf", sep="")
    output_maplot<- paste("figures/", outs$contrasts[i],".MA.pdf", sep="")
    pdf(file=output_maplot, width=7, height=6)
    p<- ggplot(data=mat, aes(x=log10(baseMean), y=log2FoldChange, col=diffexpressed)) + geom_point(size=0.5) + theme_minimal() + geom_hline(yintercept=c(-1,1), col="red") + scale_colour_manual(values = mycolors)
    print(p)
    dev.off()
}

####### Expression density plot #########
toplot = data.frame(counts(diffExp, normalized=T))
toplot = stack(toplot, select=colnames(toplot))
pdf("figures/Density_plot.sample_read_counts.pdf", width=7, height=6)
p = ggplot( toplot, aes(values, colour=ind, alpha=0.5))
p + geom_line(aes(color=ind), stat="density", alpha=0.5) +
  scale_x_log10(name="\nnormalized counts", breaks=c(0.1,1,10,100,1000,10000,100000), limits=c(0.1,100000) ) +
  scale_y_continuous(name="density\n") +
  scale_colour_discrete(name="Samples") +
  geom_vline(xintercept=10, colour="grey", linetype = "dashed") +
  theme_minimal() +
  ggtitle("Density plot\n") +
  theme()
dev.off()
rm(ddsMat); gc()

#####Comparative heatmap#####
library( "genefilter" )

sideCols=brewer.pal(12, "Set3")[colData(rld)$Replicate]

topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 100)
gene<- rownames(assay(rld)[ topVarGenes,])
fcs<-(sapply(1:length(listDiffMat_fdr),function(x)listDiffMat_fdr[[x]][gene,]$log2FoldChange))
rownames(fcs)<- gene
colnames(fcs)= outs$contrasts
pdf("figures/Heatmap_100_variable_genes_FC.pdf")
heatmap.2(fcs, scale="row", trace="none", dendrogram="row", col = colorRampPalette(rev(brewer.pal(6, "BrBG")) ), cexRow=0.4, cexCol=0.75, key=TRUE, margins=c(10,10),srtCol=45)
dev.off()
pdf("figures/Heatmap_100_variable_genes.pdf")
heatmap.2(assay(rld)[ topVarGenes,], scale="row", trace="none", dendrogram="row", col = colorRampPalette(rev(brewer.pal(9, "RdBu")) ), ColSideColors=sideCols, cexRow=0.4, cexCol=0.75, key=TRUE, margins=c(10,10),srtCol=45)
dev.off()

topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 50)
gene<- rownames(assay(rld)[ topVarGenes,])
fcs<-(sapply(1:length(listDiffMat_fdr),function(x)listDiffMat_fdr[[x]][gene,]$log2FoldChange))
rownames(fcs)<- gene
colnames(fcs)= outs$contrasts
pdf("figures/Heatmap_50_variable_genes_FC.pdf")
heatmap.2(fcs, scale="row", trace="none", dendrogram="row", col = colorRampPalette(rev(brewer.pal(6, "BrBG")) ), cexRow=0.4, cexCol=0.75, key=TRUE, margins=c(10,10),srtCol=45)
dev.off()
pdf("figures/Heatmap_50_variable_genes.pdf")
heatmap.2(assay(rld)[ topVarGenes,], scale="row", trace="none", dendrogram="row", col = colorRampPalette(rev(brewer.pal(9, "RdBu")) ), ColSideColors=sideCols, cexRow=0.4, cexCol=0.75, key=TRUE, margins=c(10,10),srtCol=45)
dev.off()

topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 25)
gene<- rownames(assay(rld)[ topVarGenes,])
fcs<-(sapply(1:length(listDiffMat_fdr),function(x)listDiffMat_fdr[[x]][gene,]$log2FoldChange))
rownames(fcs)<- gene
colnames(fcs)= outs$contrasts
pdf("figures/Heatmap_25_variable_genes_FC.pdf")
heatmap.2(fcs, scale="row", trace="none", dendrogram="row", col = colorRampPalette(rev(brewer.pal(6, "BrBG")) ), cexRow=0.4, cexCol=0.75, key=TRUE, margins=c(10,10),srtCol=45)
dev.off()
pdf("figures/Heatmap_25_variable_genes.pdf")
heatmap.2(assay(rld)[ topVarGenes,], scale="row", trace="none", dendrogram="row", col = colorRampPalette(rev(brewer.pal(9, "RdBu")) ), ColSideColors=sideCols, cexRow=0.4, cexCol=0.75, key=TRUE, margins=c(10,10),srtCol=45)
dev.off()





####Session information#######
sink("data/DESeq2.session_info.txt")
print("seed is 1234")
sessionInfo()
sink()

