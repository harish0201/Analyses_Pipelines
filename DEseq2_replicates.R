library(DESeq2)
library("BiocParallel")
register(MulticoreParam(30))
countdata = read.table("counts.matrix", header = T, row.names = NULL, check.names = F)
colnames(countdata)[1] = "Geneid";
countdata = aggregate(.~Geneid, countdata, max)
row.names(countdata) = countdata[,1]; countdata = countdata[,-1]
#name = gsub("featureCounts.","", colnames(countdata))
#colnames(countdata) = name
sampleInfo <- read.table("SampleInfo_Replicates",header = T,sep = "\t",row.names = 1)
countdata = countdata[,as.character(as.matrix(row.names(sampleInfo)))]
countdata  =  countdata[rowSums(countdata) >=  1,]
countdata = data.matrix(countdata)
ddsMat <- DESeqDataSetFromMatrix(countData = round(countdata), colData = sampleInfo, design = ~Replicate)
rld <- rlog(ddsMat, blind = FALSE)
system("mkdir -p data")
write.table(assay(rld), file = "data/readCountRlogNorm.xls", sep = "\t",col.names = NA, quote=FALSE)
diffExp <- DESeq(ddsMat)
rm(ddsMat); gc()

###############Samples-Relationships##########
library(RColorBrewer)
library(gplots)
library(circlize)
library(ComplexHeatmap)
library(Hmisc)
#distance matrix
system("mkdir -p figures/absPairVenn")
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
png(file = "figures/PCA.png", height = 10, width = 10, units="in", res=600)
plotPCA(rld, intgroup="Replicate")
dev.off()

################ COMPS##########
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
comp = read.table("contrast", header = T, row.names = NULL)
comp = data.frame("Replicate",comp[1:(length(comp)/2)],comp[((length(comp)/2)+1):length(comp)])
listDiffMat = apply(comp, 1, function(x){  results(diffExp, ,contrast = x)})
#colNames = apply(comp, 1, function(x){ paste(x[-1], collapse  =  " Vs ") })

#export out tables:
sapply(1:length(listDiffMat),function(x)write.table(listDiffMat[[x]],file=sprintf('data/Comparison%d.txt', x), quote=FALSE, sep="\t"))

###############VENNS############################
library(VennDiagram)
cols = c("red","green","blue","orange","magenta","turquoise1","pink","purple","olivedrab1")
system("mkdir -p figures/absPairVenn")
#pairwise comparision on absolute value on normalized counts
vennComb  =  combn(colnames(mCountData), 2, simplify  =  FALSE)
lapply(vennComb, function(x){
GenesList = lapply(x, function(y) { row.names(mCountData)[mCountData[,y] > 1];    })
names(GenesList) = x
name = paste(x,  collapse = "_join_");
name = paste(name, ".pdf", sep = "")
name = paste("figures/absPairVenn/", name , sep = "")
pdf(file = name)
v0 <-venn.diagram(GenesList, filename = NULL, col  =  cols[1:length(GenesList)], fill = cols[1:length(GenesList)], euler.d = FALSE, scaled = F, cat.pos = c(0,0))
grid.newpage()
grid.draw(v0)
dev.off()
})

#####Comparative heatmap#####
library( "genefilter" )
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 50)
pdf("figures/50_Heatmap_genes.pdf")
heatmap.2(assay(rld)[ topVarGenes,], scale="row", trace="none", dendrogram="row", col = colorRampPalette(rev(brewer.pal(9, "RdBu")) )(255), ColSideColors=c( Control="gray", DPN="darkgreen", OHT="orange" )[colData(rld)$Replicate], cexRow=0.3, cexCol=0.75, key=TRUE, margins=c(15,15),srtCol=45)
dev.off()
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 25)
pdf("figures/25_Heatmap_genes.pdf")
heatmap.2(assay(rld)[ topVarGenes,], scale="row", trace="none", dendrogram="row", col = colorRampPalette(rev(brewer.pal(9, "RdBu")) )(255), ColSideColors=c( Control="gray", DPN="darkgreen", OHT="orange" )[colData(rld)$Replicate], cexRow=0.4, cexCol=0.75, key=TRUE, margins=c(15,15),srtCol=45)
dev.off()
