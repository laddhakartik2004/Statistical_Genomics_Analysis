tropical=  c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')
palette(tropical)
par(pch=19)

#load
  library(devtools)
  library(Biobase)
  library(limma)
  library(edge)
  library(genefilter)

#load dataset
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file=con)
close(con)
bot = bottomly.eset
pdata=pData(bot)
edata=as.matrix(exprs(bot))
fdata = fData(bot)
ls()

#Transforming data
edata = log2(as.matrix(edata) + 1)
edata = edata[rowMeans(edata) > 10, ] #removing lowly expressed genes

#T-test (for multigroup comparison use- genefilter or rowftests)
tstats_obj = rowttests(edata,pdata$strain) # works for 2 gp comparison
names(tstats_obj)
hist(tstats_obj$statistic,col=2)

#permutation
set.seed(135) #cuz got to have a reproducible study
strain = pdata$strain
strain0 = sample(strain) #sample command for permutation
tstats_obj0 = rowttests(edata,strain0) 
hist(tstats_obj0$statistic,col=2,xlim=c(-5,2)) #repeat to see if u get anything diff
quantile(tstats_obj0$statistic) 
quantile(tstats_obj$statistic) #more symm- batch effect?

#F-statistics (lane number has 8 diff levels)
table(pdata$lane.number)
fstats_obj = rowFtests(edata,as.factor(pdata$lane.number))
names(fstats_obj)
hist(fstats_obj$statistic,col=2)

#moderated
mod = model.matrix(~ pdata$strain)
fit_limma = lmFit(edata,mod)
ebayes_limma = eBayes(fit_limma) #shrink stats
head(ebayes_limma$t)

#plotting
plot(ebayes_limma$t[,2],-tstats_obj$statistic,col=4, #tstats multiplied by minus
      xlab="Moderated T-stat",ylab="Regular T-stat")  #same direction
abline(c(0,1),col="darkgrey",lwd=3)

#for all 8 levels of lane.number
mod_adj = model.matrix(~ pdata$strain + as.factor(pdata$lane.number))
fit_limma_adj = lmFit(edata,mod_adj)
ebayes_limma_adj = eBayes(fit_limma_adj)
head(ebayes_limma_adj$t)

plot(ebayes_limma_adj$t[,2],-tstats_obj$statistic,col=3,
     xlab="Adjusted T-stat",ylab="T-stat")
abline(c(0,1),lwd=3,col="darkgrey")

#Factor
mod_lane = model.matrix(~ as.factor(pdata$lane.number))
fit_limma_lane = lmFit(edata,mod_lane)
ebayes_limma_lane = eBayes(fit_limma_lane) 
head(ebayes_limma_lane$t)

#comparing diff lanes
top_lane = topTable(ebayes_limma_lane, coef=2:7,
                    number=dim(edata)[1],sort.by="none")
head(top_lane)
plot(top_lane$F,fstats_obj$statistic,
     xlab="Moderated F-statistic",ylab="F-statistic",col=3)

#Using edge
edge_study = build_study(edata, grp = as.factor(pdata$lane.number))
de_obj = lrt(edge_study)
qval = qvalueObj(de_obj)
plot(qval$stat,fstats_obj$statistic,col=4,
      xlab="F-stat from edge",ylab="F-stat from genefilter")

#Adjust for strain
edge_study2 = build_study(edata, grp = as.factor(pdata$lane.number),
                        adj.var=pdata$strain)
de_obj2 = lrt(edge_study2)
qval2 = qvalueObj(de_obj2)
plot(qval2$stat,fstats_obj$statistic,col=4,
      xlab="Adj F-stat from edge",ylab="F-stat from genefilter")
