# Genomic Data Analysis
## Johns Hopkins University Coursework (Coursera) 
## Overview
This repository contains my work and learning outcomes from the Statistics for Genomic Data Science course offered by Johns Hopkins University (Prof Jeff Leek). The project focuses on applying computational and statistical methods to analyze high-throughput genomic datasets and extract biologically meaningful insights. Scripts are rewritten and annotated to reflect my understanding.
## Dataset
The datasets used in this project are derived from publicly available gene expression repositories.

Bodymap dataset- http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData

Montgomery-Pickrell dataset- http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData

Bottomly dataset- http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData

### for goseq and DESeq2
Install the packages from Bioconductor manager using

if (!requireNamespace("BiocManager", quietly = TRUE))

    install.packages("BiocManager")


BiocManager::install("goseq")

Similarly, replace the package name with DESeq2 to download it. To load system files from goseq,

temp_data =read.table(system.file("extdata","Li_sum.txt",

                                    package="goseq"),sep="\t",
                                    
                                    header=TRUE,
                                    
                                    stringsAsFactors=FALSE)

### for eQTL (Expression Quantitative Trait Loci)
Install MatrixEQTL using

install.packages("MatrixEQTL")

Then load sample data from the base directory using

base.dir = find.package("MatrixEQTL")

SNP_file_name = paste(base.dir, "/data/SNP.txt", sep="");

expression_file_name = paste(base.dir, "/data/GE.txt", sep="")

covariates_file_name = paste(base.dir, "/data/Covariates.txt", sep="")


## Tools
1. RStudio
2. Bioconductor
3. bowtie (datasets)

## R Packages required 
### (most from bioconductor)
gplots, devtools, Biobase, dplyr, RSkittleBrewer, org.Hs.eg.db, Annotationdbi, dendextend,
preprocessCore, broom, snpStats, sva, bladderbatch, limma, edge, genefilter, qvalue,
DESeq2, goseq, MatrixEQTL

## Clone my repository
git clone https://github.com/laddhakartik2004/Statistical_Genomics_Analysis.git
