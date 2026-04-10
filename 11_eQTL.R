tropical=  c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')
palette(tropical)
par(pch=19)

install.packages("MatrixEQTL")
#load
  library(devtools)
  library(Biobase)
  library(MatrixEQTL)

#loading data from the package
base.dir = find.package("MatrixEQTL")
SNP_file_name = paste(base.dir, "/data/SNP.txt", sep="");
expression_file_name = paste(base.dir, "/data/GE.txt", sep="")
covariates_file_name = paste(base.dir, "/data/Covariates.txt", sep="")
output_file_name = tempfile()

#reading expression, snp and covariate data
expr = read.table(expression_file_name,sep="\t",
                  header=T,row.names=1)
expr[1,]

snps = read.table(SNP_file_name,sep="\t",
                  header=T,row.names=1)
snps[1,]

cvrt = read.table(covariates_file_name,sep="\t",
                  header=T,row.names=1)
cvrt[1,]
dim(cvrt) #2 cvrts exist (gender, age)
cvrt[2,]

#linear model for expression and snp
e1 = as.numeric(expr[1,]) 
s1 = as.numeric(snps[1,])
lm1 = lm(e1 ~ s1)
tidy(lm1) #p value large, not highly assoc

#Plots
plot(e1 ~ jitter(s1),
     col=(s1+1),xaxt="n",xlab="Genotype",ylab="Expression")
axis(1,at=c(0:2),labels=c("AA","Aa","aa"))
lines(lm1$fitted ~ s1,type="b",pch=15,col="darkgrey") #plot fitted values

dim(expr)

#parameters
pvOutputThreshold = 1e-2 #p value threshold
errorCovariance = numeric() 
useModel = modelLINEAR

#
snps = SlicedData$new()        #sliced data obj
snps$fileDelimiter = "\t"      # the TAB character
snps$fileOmitCharacters = "NA" # denote missing values;
snps$fileSkipRows = 1          # one row of column labels
snps$fileSkipColumns = 1       # one column of row labels
snps$fileSliceSize = 2000      # read file in chunks of 2,000
snps$LoadFile( SNP_file_name )

#
gene = SlicedData$new()
gene$fileDelimiter = "\t"       # the TAB character
gene$fileOmitCharacters = "NA"  # denote missing values;
gene$fileSkipRows = 1           # one row of column labels
gene$fileSkipColumns = 1        # one column of row labels
gene$fileSliceSize = 2000       # read file in pieces of 2,000 rows
gene$LoadFile(expression_file_name)

#empty obj, so no adjustment
cvrt = SlicedData$new()

#RUNNING MatrixEQTL
me = Matrix_eQTL_engine(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = NULL,
    pvOutputThreshold = pvOutputThreshold,
    useModel = useModel, 
    errorCovariance = errorCovariance, 
    verbose = TRUE,
    pvalue.hist = TRUE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE); #calculate FDR or not

plot(me)

names(me)
me$time.in.sec #time elapsed
me$param #all diff parameters
me$all$neqtls #no. of eqtls passed threshold=1
me$all$eqtls #tells eqtl snp, gene, p val, FDR
