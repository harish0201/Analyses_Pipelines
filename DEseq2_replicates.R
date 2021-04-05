set.seed("1234")
library(DESeq2)
library("BiocParallel")
args = commandArgs(trailingOnly=TRUE)
register(MulticoreParam(4))
countdata = read.table(args[1], header = T, row.names = NULL, check.names = F)
colnames(countdata)[1] = "Geneid";
countdata = aggregate(.~Geneid, countdata, max)
row.names(countdata) = countdata[,1]; countdata = countdata[,-1]
sampleInfo <- read.table(args[2], header = T,sep = "\t", row.names = 1)
countdata = countdata[,as.character(as.matrix(row.names(sampleInfo)))]
countdata  =  countdata[rowSums(countdata) >=  5,]
countdata = data.matrix(countdata)
ddsMat <- DESeqDataSetFromMatrix(countData = round(countdata), colData = sampleInfo, design = ~Replicate)
rld <- rlog(ddsMat, blind = FALSE)
system("mkdir -p data")
write.table(assay(rld), file = "data/readCountRlogNorm.xls", sep = "\t",col.names = NA, quote=FALSE)
diffExp <- DESeq(ddsMat)

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
system("mkdir -p figures")
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
ggsave(file=sprintf("figures/PCA_var1_var2.pdf"), width=7, height=6)


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
comp = data.frame("Replicate",comp[1:(length(comp)/2)],comp[((length(comp)/2)+1):length(comp)])
listDiffMat_nofdr = apply(comp, 1, function(x){  results(diffExp, contrast = x)})
#colNames = apply(comp, 1, function(x){ paste(x[-1], collapse  =  " Vs ") })

#export out tables:
sapply(1:length(listDiffMat_nofdr),function(x)write.table(listDiffMat_nofdr[[x]],file=sprintf('data/Comparison%d.xls', x), quote=FALSE, sep="\t"))

####### Expression density plot #########
toplot = data.frame(counts(diffExp, normalized=T))
toplot = stack(toplot, select=colnames(toplot))
p = ggplot( toplot, aes(values, colour=ind, alpha=0.5))
p + geom_line(aes(color=ind), stat="density", alpha=0.5) +
  scale_x_log10(name="\nnormalized counts", breaks=c(0.1,1,10,100,1000,10000,100000), limits=c(0.1,100000) ) +
  scale_y_continuous(name="density\n") +
  scale_colour_discrete(name="Samples") +
  geom_vline(xintercept=10, colour="grey", linetype = "dashed") +
  theme_minimal() +
  ggtitle("Density plot\n") +
  theme()
ggsave(file=sprintf("figures/Density_plot.sample_read_counts.pdf"), width=7, height=6)
rm(ddsMat); gc()

#####Comparative heatmap#####
library( "genefilter" )
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 100)
pdf("figures/100_Heatmap_genes.pdf")
dev.off()
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 50)
pdf("figures/50_Heatmap_genes.pdf")
heatmap.2(assay(rld)[ topVarGenes,], scale="row", trace="none", dendrogram="row", col = colorRampPalette(rev(brewer.pal(9, "RdBu")) )(255), ColSideColors=c( Control="gray", DPN="darkgreen", OH
T="orange" )[colData(rld)$Replicate], cexRow=0.3, cexCol=0.75, key=TRUE, margins=c(15,15),srtCol=45)
dev.off()
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 25)
pdf("figures/25_Heatmap_genes.pdf")
heatmap.2(assay(rld)[ topVarGenes,], scale="row", trace="none", dendrogram="row", col = colorRampPalette(rev(brewer.pal(9, "RdBu")) )(255), ColSideColors=c( Control="gray", DPN="darkgreen", OH
T="orange" )[colData(rld)$Replicate], cexRow=0.4, cexCol=0.75, key=TRUE, margins=c(15,15),srtCol=45)
dev.off()


####Session information#######
sink("DESeq2.session_info.txt")
sessionInfo()
sink()
system("rm Rplots.pdf")
