#### Runs of Homozygosity
#### Assumes variations have been called from Resequencing analyses and sex chromosomes have been identified per PSMC analyses

#1. Convert into compatible format for RZooRoH
bcftools convert -g outfile –-chrom –-tag GT variants.final.vcf
sed -e ’s/-nan/0/g’ outfile.gen > newfile.gen

#2. Read the same in R
library(RZooRoH)
A <- "newfile.gen"
IDS<- "sample.ids"
B <- zoodata(genofile = A, zformat = "gp", samplefile = IDs)
#model with pre-defined rates for 5 classes (4 HBD and 1 non-HBD# class)
mix5R <- zoomodel(predefined = TRUE, K=6, krates = c(4,8,16,32,64,64)) #(>22000 SNPs means we can go for k>6)
mix5R
ZR1 <- zoorun(mix5R,B,localhbd = TRUE)

ZR1@mixc
#will give the results for each sample in a row for each K in column.
```
[,1]         [,2]         [,3]         [,4]      [,5]      [,6]
[1,] 1.612398e-05 1.774601e-05 2.209367e-05 4.539151e-05 0.3072168 0.6926819
[2,] 1.610590e-05 1.770317e-05 2.201947e-05 4.563731e-05 0.3156164 0.6842821
```

ZR1@krates 
```
[,1] [,2] [,3] [,4] [,5] [,6]
[1,]    4    8   16   32   64   64
[2,]    4    8   16   32   64   64
```
###Rate class inbreeding
x <- 1-ZR1@realized[,6]
hist(x,main="",xlab="Inbreeding coefficient",col='tomato', nc=30, xlim=c(0,0.15))

###Boxplot: lower the box, => inbreeding has occured
zooplot_prophbd(ZR1, cols = 'red4', style = 'boxplot')

#3. Regions with plink:
##1. Create .ped and .map files
plink --vcf 6.vcf.gz --const-fid --allow-extra-chr  --recode12

##2. Runs of Homozygosity #new_cmds
plink --file 9 --allow-extra-chr --homozyg-window-snp 50 --homozyg-window-threshold 0.05 --homozyg-kb 500 --homozyg-window-missing 2 --homozyg-window-het 1 --homozyg-gap 1000 --homozyg-density 50 --out 9.out
