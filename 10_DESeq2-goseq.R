tropical=  c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')
palette(tropical)
par(pch=19)

#load packages (DESeq2, goseq)
  library(devtools)
  library(Biobase)
  library(goseq)
  library(DESeq2)

head(supportedGenomes()) #details about first six databases
head(supportedGeneIDs()) 

#loading system file
temp_data =read.table(system.file("extdata","Li_sum.txt",
                                     package="goseq"),sep="\t",
                                     header=TRUE,
                                     stringsAsFactors=FALSE)

#first column with ensembl ids labelled as 'genes'
#removing the column and setting rownames as gene names

temp_data[1,] 
expr= temp_data[,-1]
rownames(expr) = temp_data[,1]

expr = expr[rowMeans(expr) > 5,] #removing lowly expressed genes
grp=factor(rep(c("Control","Treated"),times=c(4,3)))#set gp variables
grp
pdata  = data.frame(grp)

#DESeq
de = DESeqDataSetFromMatrix(expr, pdata, ~grp)
de_fit = DESeq(de) #identify differentially expressed genes
de_results = results(de_fit)
head(de_results) #p value and adj p value for every gene

#finding genes that are differentially expressed (FDiscoveryR)
genes = as.integer(de_results$padj < 0.05)
not_na = !is.na(genes) #removing those with not enough data
names(genes) = rownames(expr)
genes = genes[not_na]

## ------------------------------------------------------------------------
head(supportedGenomes(),n=12)[,1:4]

#probability weight function (hg19 genome)
pwf=nullp(genes,"hg19","ensGene")
head(pwf)

#using goseq for enrichment statistics (GO for gene ontology)
GO.wall=goseq(pwf,"hg19","ensGene")
head(GO.wall) #overrepr, underrepr p values, numDE

#for particular Molecular Function categ
GO.MF=goseq(pwf,"hg19","ensGene",test.cats=c("GO:MF"))
head(GO.MF)
