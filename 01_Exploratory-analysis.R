#load
  library(gplots)
  library(devtools)
  library(Biobase)
  library(dplyr)
  library(RSkittleBrewer)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")
  library(org.Hs.eg.db)
  library(AnnotationDbi)

tropical=  c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')
palette(tropical)
par(pch=19)

#load_data
con = url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
pdata=pData(bm)
edata=exprs(bm)
fdata = fData(bm)
ls()

#tables
table(pdata$gender)
table(pdata$gender,pdata$race)

#summary
summary(edata)

# To see missing values
# Use option useNA to include NA's in table
table(pdata$age,useNA="ifany")

# is.na check for NA values
table(is.na(pdata$age))

# Checking for other common missing names (blank)
sum(pdata$age==" ")
# to remove NA values
sum(pdata$age==" ",na.rm=TRUE)

# Check genomic data for NAs
is.na(edata)[1,]
sum(is.na(edata))

# Make the distribution of NA's by genes
gene_na = rowSums(is.na(edata))
table(gene_na)

# Make the distribution of NA's by samples
sample_na = colSums(is.na(edata))
table(sample_na)


#dimensions (row-pdata=col-edata, row-edata=row-fdata)
dim(pdata)
dim(edata)
dim(fdata)

#For making boxplot
boxplot(edata[,1]) #not good
boxplot(log2(edata[,1]+1)) #not good either
boxplot(log2(edata+1),col=2,range=0)

#Histograms
par(mfrow=c(1,2))
hist(log2(edata[,1]+1),col=2)
hist(log2(edata[,2]+1),col=2) #sample by sample

#Density plot
plot(density(log2(edata[,1]+1)),col=2)
lines(density(log2(edata[,2]+1)),col=3) #overlays plot 2 over 1

qqplot(log2(edata[,1]+1), log2(edata[,2]+1),col=3)
abline(c(0,1))

#Bland-altman
mm = log2(edata[,1]+1) - log2(edata[,2]+1)
aa = log2(edata[,1]+1) + log2(edata[,2]+1)
plot(aa,mm,col=2)

#remove low-count data
edata = as.data.frame(edata)
filt_edata = filter(edata,rowMeans(edata) > 1)
dim(filt_edata)
boxplot(as.matrix(log2(filt_edata+1)),col=2)

#Checking for consistency (dataset vs meta-data)
aeid = as.character(fdata[,1])
aeid[1:5] #ensembl ids
chr = AnnotationDbi::select(org.Hs.eg.db,keys=aeid,keytype="ENSEMBL",columns="CHR")
head(chr)

dim(chr)
dim(edata)
# Take non-duplicated chromsomes
chr = chr[!duplicated(chr[,1]),]

# Confirm that the annotation still is in the right order
all(chr[,1] == rownames(edata))

# Select the chromosome Y samples
edata=as.data.frame(edata)
edatay = dplyr::filter(edata,chr$CHR=="Y")
dim(edatay)
# Males have Y chromsome expression as expected
boxplot(colSums(edatay) ~ pdata$gender)
points(colSums(edatay) ~ jitter(as.numeric(pdata$gender)),
        col=as.numeric(pdata$gender),
        pch=19)


#Multi-variate
ematrix = as.matrix(edata)[rowMeans(edata) > 10000,]
dim(ematrix)
heatmap(ematrix)
colramp = colorRampPalette(c(3,"white",2))(9)
colramp
heatmap(ematrix,col=colramp)
heatmap(ematrix,col=colramp,Rowv=NA,Colv=NA)
heatmap.2(ematrix,col=colramp,Rowv=NA,Colv=NA,
          dendrogram="none", scale="row",trace="none")
