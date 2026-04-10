tropical=  c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')
palette(tropical)

  library(devtools)
  library(Biobase)

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
pdata=pData(bm)
edata=as.data.frame(exprs(bm))
fdata = fData(bm)
ls()

hist(rnorm(1000),col=2)
hist(edata[,1],col=2,breaks=100) #mostly zeroes in dataset
hist(log(edata[,1]),col=2,breaks=100)
min(log(edata))
quantile(log(edata[,1]))

#adding 1 to obs doesn't change already large obs but helps with zeroes
min(log(edata[,1] + 1))
hist(log(edata[,1] + 1),breaks=100,col=2)

#for better visualisation
hist(log2(edata[,1] + 1),breaks=100,col=2,xlim=c(1,15),ylim=c(0,400))

hist(rowSums(edata==0),col=2)

#remove low values (lowly expressed genes)
low_genes = rowMeans(edata) < 5
table(low_genes)
filt_edata = filter(edata,!low_genes)
dim(filt_edata)

low_genes2 = rowMedians(as.matrix(edata)) < 5
table(low_genes2,low_genes)
filt_edata2 = filter(edata,!low_genes2)
dim(filt_edata2)

hist(log2(filt_edata[,1] + 1),col=2)
