tropical=  c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')
palette(tropical)
par(pch=19)

#load
  library(devtools)
  library(Biobase)

#loading data
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)
ls()

#reducing size of dataset (filter)
edata = edata[rowMeans(edata) > 100, ]
dim(edata)
edata = log2(edata + 1)
edata_centered = edata - rowMeans(edata) #otherwise 1st singular value of vector will always be mean level
                                         #but we want to see variation between samples or genes
svd1 = svd(edata_centered)
names(svd1)
svd1$d #has 129 values as 129 columns in edata
dim(edata)
dim(svd1$v) #variation across genes
dim(svd1$u) #variation across samples

#plotting
plot(svd1$d,ylab="Singular value",col=2)

#plotting variance
plot(svd1$d^2/sum(svd1$d^2),ylab="Percent Variance Explained",col=2)

#for first two pc/sv
par(mfrow=c(1,2))
plot(svd1$v[,1],col=2,ylab="1st PC")
plot(svd1$v[,2],col=2,ylab="2nd PC")
plot(svd1$v[,1],svd1$v[,2],col=2,ylab="2nd PC",xlab="1st PC")

#for comparing both studies
plot(svd1$v[,1],svd1$v[,2],ylab="2nd PC",
     xlab="1st PC",col=as.numeric(pdata$study))

boxplot(svd1$v[,1] ~ pdata$study,border=c(1,2))
points(svd1$v[,1] ~ jitter(as.numeric(pdata$study)),col=as.numeric(pdata$study))

#
pc1 = prcomp(edata)
plot(pc1$rotation[,1],svd1$v[,1])

edata_centered2 = t(t(edata) - colMeans(edata))
svd2 = svd(edata_centered2)
plot(pc1$rotation[,1],svd2$v[,1],col=2) 

#
edata_outlier = edata_centered
edata_outlier[6,] = edata_centered[6,] * 10000
svd3 = svd(edata_outlier)
plot(svd1$v[,1],svd3$v[,1],xlab="Without outlier",ylab="With outlier")

#
plot(svd3$v[,1],edata_outlier[6,],col=4)