tropical=  c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')
palette(tropical)

#load
  library(devtools)
  library(Biobase)
  library(limma)
  library(edge)
  library(genefilter)
  library(qvalue)

#load dataset
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file=con)
close(con)
bot = bottomly.eset
pdata=pData(bot)
edata=as.matrix(exprs(bot))
fdata = fData(bot)
ls()

#transforming data
edata = log2(as.matrix(edata) + 1)
edata = edata[rowMeans(edata) > 10, ]

#p value from Fstats
fstats_obj = rowFtests(edata,as.factor(pdata$strain))
head(fstats_obj)
hist(fstats_obj$p.value,col=2)

#Using edge package
edge_study = build_study(edata, grp = pdata$strain, 
                         adj.var = as.factor(pdata$lane.number))
de_obj = lrt(edge_study) #likelihood ratio test
qval = qvalueObj(de_obj)
hist(qval$pvalues,col=3)
#strange p value (missing a variable or modeling strategy may not be right)

#Moderated p values
mod = model.matrix(~ pdata$strain + pdata$lane.number)
fit_limma = lmFit(edata,mod)
ebayes_limma = eBayes(fit_limma)
limma_pvals = topTable(ebayes_limma,number=dim(edata)[1])$P.Value
hist(limma_pvals,col=4)

#empirically calculating p values
set.seed(3333)
B = 1000
tstats_obj = rowttests(edata,pdata$strain)
tstat0 = matrix(NA,nrow=dim(edata)[1],ncol=B)
tstat = tstats_obj$statistic
strain = pdata$strain
for(i in 1:B){
  strain0 = sample(strain)
  tstat0[,i] = rowttests(edata,strain0)$statistic
} #repermutes strain data every time it loops

emp_pvals = empPvals(tstat,tstat0)
hist(emp_pvals,col=2)

#using Bonferroni correction (FWER, 1fp)
fp_bonf = p.adjust(fstats_obj$p.value,method="bonferroni")
hist(fp_bonf,col=3)
quantile(fp_bonf)
sum(fp_bonf <0.05) #no statistically sig results

#Benjamini-Hochberg correction (FDR)
fp_bh = p.adjust(fstats_obj$p.value,method="BH")
hist(fp_bh,col=3)
quantile(fp_bh)
sum(fp_bh <0.05) #none sig

#limma package
limma_pvals_adj = topTable(ebayes_limma,number=dim(edata)[1])$adj.P.Val
hist(limma_pvals_adj,col=2)
sum(limma_pvals_adj <0.05) #2 adj values

#qvalue package
qval_limma = qvalue(limma_pvals)
summary(qval_limma) #tells diff threshold levels
qval$pi0 #estimated fraction of null hypothesis

#qvalues for edge obj
qval = qvalueObj(de_obj)
summary(qval)