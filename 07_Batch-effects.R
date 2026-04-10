tropical=  c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')
palette(tropical)
par(pch=19)

#load
  library(devtools)
  library(Biobase)
  library(sva)
  library(bladderbatch)
  library(snpStats)

#loading data
data(bladderdata)
pheno = pData(bladderEset)
edata = exprs(bladderEset)
dim(edata)
dim(pheno)
head(pheno)

#as we have the batch variable, making a model matrix
mod = model.matrix(~as.factor(cancer) + as.factor(batch),data=pheno)
fit = lm.fit(mod,t(edata)) #coeffs adjusted for batch effects 
hist(fit$coefficients[2,],col=2,breaks=100) 

table(pheno$cancer,pheno$batch) #check correl b/w batch variable and outcome variable

#Another approach to adjust batch effects (batch-combat)
batch = pheno$batch
modcombat = model.matrix(~1, data=pheno)
modcancer = model.matrix(~cancer, data=pheno)
combat_edata = ComBat(dat=edata, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
combat_fit = lm.fit(modcancer,t(combat_edata))
hist(combat_fit$coefficients[2,],col=2,breaks=100)

#original fit v/s combat fit plot
plot(fit$coefficients[2,],combat_fit$coefficients[2,],col=2,
      xlab="Linear Model",ylab="Combat",xlim=c(-5,5),ylim=c(-5,5))
abline(c(0,1),col=1,lwd=3)

#sva package
mod = model.matrix(~cancer,data=pheno)
mod0 = model.matrix(~1, data=pheno)
sva1 = sva(edata,mod,mod0,n.sv=2) #surrogate variables=2
names(sva1) #new covariate sv created
summary(lm(sva1$sv ~ pheno$batch))

#2nd surrogate variable is correl with batch variable, so
boxplot(sva1$sv[,2] ~ pheno$batch)
points(sva1$sv[,2] ~ jitter(as.numeric(pheno$batch)),col=as.numeric(pheno$batch))

#sva identified new covariates
modsv = cbind(mod,sva1$sv)
modsv
fitsv = lm.fit(modsv,t(edata))

par(mfrow=c(1,2))

#comparing SVA to Combat
plot(fitsv$coefficients[2,],combat_fit$coefficients[2,],col=2,
      xlab="SVA",ylab="Combat",xlim=c(-5,5),ylim=c(-5,5))
abline(c(0,1),col=1,lwd=3)\

#comparing SVA to linear model (less shrunken,more correl)
plot(fitsv$coefficients[2,], fit$coefficients[2,],col=2,
      xlab="SVA",ylab="linear model",xlim=c(-5,5),ylim=c(-5,5))
abline(c(0,1),col=1,lwd=3)

#snpStats
data(for.exercise)
controls <- rownames(subject.support)[subject.support$cc==0]
controls #diff populations
use <- seq(1, ncol(snps.10), 10)
ctl.10 <- snps.10[controls,use]

#Calculate principal components
xxmat <- xxt(ctl.10, correct.for.missing=FALSE) #linear algebra calc
evv <- eigen(xxmat, symmetric=TRUE) #decomposition
pcs <- evv$vectors[,1:5]
pcs[1,]

#Comparing pcs
pop <- subject.support[controls,"stratum"]
plot(pcs[,1],pcs[,2],col=as.numeric(pop),
      xlab="PC1",ylab="PC2")
legend(0,0.15,legend=levels(pop),pch=19,col=1:2)