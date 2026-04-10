tropical=  c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')
palette(tropical)
par(pch=19)

#load
  library(devtools)
  library(Biobase)
  library(broom)

#dataset
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
pdata=pData(bm)
edata=as.data.frame(exprs(bm))
fdata = fData(bm)
ls()

#tidy data
edata = as.matrix(edata)
lm1 = lm(edata[1,] ~ pdata$age)
lm1
tidy(lm1)

#plotting
plot(pdata$age,edata[1,], col=1)
abline(lm1$coeff[1],lm1$coeff[2], col=2,lwd=3)


pdata$gender
table(pdata$gender)

#association between gender and gene expression (not very strong here)
boxplot(edata[1,] ~ pdata$gender)
points(edata[1,] ~ jitter(as.numeric(pdata$gender)),
       col=as.numeric(pdata$gender))

#dummy variables (equal 1 if True)-
dummy_m = pdata$gender=="M"
dummy_m
dummy_m*1
dummy_f = pdata$gender=="F"
dummy_f
dummy_f*1

lm2 = lm(edata[1,] ~ pdata$gender)
tidy(lm2)

mod2 = model.matrix(~pdata$gender)
mod2

#Tissue type
table(pdata$tissue.type)
pdata$tissue.type == "adipose"
pdata$tissue.type == "adrenal"

#intercept term represents tissue type adipose (doesn't appear in list of coeff)
tidy(lm(edata[1,] ~ pdata$tissue.type ))

#adjusting for variables
lm3 = lm(edata[1,] ~ pdata$age + pdata$gender)
tidy(lm3)

#interaction model(for females, genderM term is zero)
lm4 = lm(edata[1,] ~ pdata$age*pdata$gender)
tidy(lm4)

#Overlaying plots
lm4 = lm(edata[6,] ~ pdata$age)
plot(pdata$age,edata[6,],col=2) #see outlier
abline(lm4,col=1,lwd=3) #outlier doesn't affect the line mucn


index = 1:19
lm5 = lm(edata[6,] ~ index)
plot(index,edata[6,],col=2)
abline(lm5,col=1,lwd=3) #outlier affects the line

lm6 = lm(edata[6,-19] ~ index[-19]) #subtracting outlier
abline(lm6,col=3,lwd=3)

legend(5,1000,c("With outlier","Without outlier"),col=c(1,3),lwd=3)


#residuals
par(mfrow=c(1,2))
hist(lm6$residuals,col=2)
hist(lm5$residuals,col=3)

#plotting residuals colored by different variables
colramp = colorRampPalette(1:4)(17)
lm9 = lm(edata[2,] ~ pdata$age)
plot(lm9$residuals,col=colramp[as.numeric(pdata$tissue.type)])
