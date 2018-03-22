BoneNotebook
================
Matthew Schmitz
March 12, 2018

This is an analysis for Candice where we're going to compare the RNAseq data she has to the seq data from the Ambrosi et al data.

Ambrosi et al. Data Analysis
============================

This data was taken from ENA with accession number ERP013883 (<http://> www.ebi.ac.uk/ena). Following quality control, reads were quantified using Salmon.

First I'll load the libraries and functions I'll need for the analysis (skip this part).

``` r
library(DESeq2)
library(EBSeq)
library(matrixStats)
library(biomaRt)
library(gplots)
library(clusterProfiler)
library(RColorBrewer)
library(tximport)
cols <-  colorRampPalette(rev(brewer.pal(11,"RdBu")))(50)

#see Github lengning
library(EACI)


median.normalize <- function(x){
  GetNormalizedMat(x,MedianNorm(x))
}

round.log <- function(s,base=2){
  round(log(s+1, base),digits = 1)
}

std.heatmap <- function(M,...){
  heatmap.2(M,Rowv = F,Colv = F,trace="none",col = cols,...)
}

rn.merge <- function(x,y,fill=0,simple.intersect=F){
  rn <- intersect(rownames(x),rownames(y))
  zerosx <- setdiff(rownames(x),rownames(y))
  zerosy <- setdiff(rownames(y),rownames(x))
  out <- cbind(x[rn,,drop=F],y[rn,,drop=F])
  if(simple.intersect){return(out)}
  else{
  if(length(zerosx)!=1  &length(zerosy)!=1){
    zx <- matrix(fill, nrow=length(zerosx), ncol =ncol(y), dimnames = list(zerosx,NULL))
    zy <- matrix(fill, nrow=length(zerosy), ncol =ncol(x), dimnames = list(zerosy,NULL))
    zx <- cbind(x[zerosx,],zx)
    zy <- cbind(zy,y[zerosy,])
  }else if(length(zerosx)==1){
    zx <- rep(fill, ncol(y))
    zy <- matrix(fill, nrow=length(zerosy), ncol =ncol(x), dimnames = list(zerosy,NULL))
    zx <- c(x[zerosx,],zx)
    zy <- cbind(zy,y[zerosy,])
  }else if(length(zerosy)==1){
    zx <- matrix(fill, nrow=length(zerosx), ncol =ncol(y), dimnames = list(zerosx,NULL))
    zy <- rep(fill, ncol(x))
    print(zx)
    print(zy)
    zx <- cbind(x[zerosx,],zx)
    zy <- c(zy,y[zerosy,])
  }
  out <- rbind(out,rbind(zx,zy))
  return(out)}
}

cor.compare <- function(x,y,min=0, varX=NULL ,interest.set =NULL, ...){
  d <- rn.compare(x,y)
  x <- d[[1]]
  y <- d[[2]]
  i = intersect(rownames(x), rownames(y))
  i = i[rowMaxs(x[i,],na.rm = T)>=min | rowMaxs(y[i,],na.rm = )>=min]
  if(!is.null(interest.set)){
    i = interest.set[interest.set%in%i]
  }
  if(!is.null(varX)){
    i = i[order(rowMeans(cbind(rowSds(x[i,]), rowSds(y[i,]))),decreasing = T)]
    i = i[1:ifelse(varX>length(i),length(i),varX)]
  }
  print("Num Genes:")
  print(length(i))
  return(cor(as.matrix(x[i,]),as.matrix(y[i,]), ...))
}

rn.compare <- function(x,y,fill=0){
  rn <- as.character(intersect(rownames(x),rownames(y)))
  zerosx <- setdiff(rownames(x),rownames(y))
  zerosy <- setdiff(rownames(y),rownames(x))
  zx <- matrix(fill, nrow=length(zerosx), ncol =ncol(y), dimnames = list(zerosx,NULL))
  zy <- matrix(fill, nrow=length(zerosy), ncol =ncol(x), dimnames = list(zerosy,NULL))
  nx <- rbind(x,zy)
  ny <- rbind(y,zx)
  return(list(nx[rownames(ny),,drop=F],ny[rownames(nx),,drop=F]))
}
```

``` r
datapath <- "~/code/data/GEOData/seq/AmbrosiBone/out/"
fileList <- dir(datapath)
fileList <- paste0(datapath,"/",fileList,"/","quant.sf")

dsList <- lapply(fileList,read.csv2, sep="\t",header=T,row.names=1,stringsAsFactors=F)
allRownames <- Reduce(union,lapply(dsList,rownames))

#Use Biomart to get the 
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host="www.ensembl.org")
rnSymbol <- getBM(attributes = c("ensembl_transcript_id_version","mgi_symbol","description"),filters = c("ensembl_transcript_id_version"),values =allRownames ,mart = mart) 
rnSymbolGenes <- rnSymbol[rnSymbol$mgi_symbol!="",]
rnSymbolGenes <- rnSymbol[rnSymbol$mgi_symbol!=""& !grepl("predicted gene", rnSymbol$description),]

#Load the table from ENA with the names of each sample
sampleMat <- read.table(file = "~/code/data/GEOData/seq/AmbrosiBone/PRJEB12408.txt",sep="\t",header = T,stringsAsFactors = F)
print(sampleMat)
```

    ##    study_accession sample_accession secondary_sample_accession
    ## 1       PRJEB12408     SAMEA3724741                 ERS1031890
    ## 2       PRJEB12408     SAMEA3724742                 ERS1031891
    ## 3       PRJEB12408     SAMEA3724743                 ERS1031892
    ## 4       PRJEB12408     SAMEA3724744                 ERS1031893
    ## 5       PRJEB12408     SAMEA3724745                 ERS1031894
    ## 6       PRJEB12408     SAMEA3724746                 ERS1031895
    ## 7       PRJEB12408     SAMEA3724747                 ERS1031896
    ## 8       PRJEB12408     SAMEA3724748                 ERS1031897
    ## 9       PRJEB12408     SAMEA3724749                 ERS1031898
    ## 10      PRJEB12408     SAMEA3724750                 ERS1031899
    ## 11      PRJEB12408     SAMEA3724751                 ERS1031900
    ## 12      PRJEB12408     SAMEA3724752                 ERS1031901
    ## 13      PRJEB12408     SAMEA3724741                 ERS1031890
    ## 14      PRJEB12408     SAMEA3724742                 ERS1031891
    ## 15      PRJEB12408     SAMEA3724743                 ERS1031892
    ## 16      PRJEB12408     SAMEA3724744                 ERS1031893
    ## 17      PRJEB12408     SAMEA3724745                 ERS1031894
    ## 18      PRJEB12408     SAMEA3724746                 ERS1031895
    ## 19      PRJEB12408     SAMEA3724747                 ERS1031896
    ## 20      PRJEB12408     SAMEA3724748                 ERS1031897
    ## 21      PRJEB12408     SAMEA3724749                 ERS1031898
    ## 22      PRJEB12408     SAMEA3724750                 ERS1031899
    ## 23      PRJEB12408     SAMEA3724751                 ERS1031900
    ## 24      PRJEB12408     SAMEA3724752                 ERS1031901
    ##    experiment_accession run_accession tax_id scientific_name
    ## 1            ERX1425553    ERR1354078  10090    Mus musculus
    ## 2            ERX1425554    ERR1354079  10090    Mus musculus
    ## 3            ERX1425555    ERR1354080  10090    Mus musculus
    ## 4            ERX1425556    ERR1354081  10090    Mus musculus
    ## 5            ERX1425557    ERR1354082  10090    Mus musculus
    ## 6            ERX1425558    ERR1354083  10090    Mus musculus
    ## 7            ERX1425559    ERR1354084  10090    Mus musculus
    ## 8            ERX1425560    ERR1354085  10090    Mus musculus
    ## 9            ERX1425561    ERR1354086  10090    Mus musculus
    ## 10           ERX1425562    ERR1354087  10090    Mus musculus
    ## 11           ERX1425563    ERR1354088  10090    Mus musculus
    ## 12           ERX1425564    ERR1354089  10090    Mus musculus
    ## 13           ERX1425565    ERR1354090  10090    Mus musculus
    ## 14           ERX1425566    ERR1354091  10090    Mus musculus
    ## 15           ERX1425567    ERR1354092  10090    Mus musculus
    ## 16           ERX1425568    ERR1354093  10090    Mus musculus
    ## 17           ERX1425569    ERR1354094  10090    Mus musculus
    ## 18           ERX1425570    ERR1354095  10090    Mus musculus
    ## 19           ERX1425571    ERR1354096  10090    Mus musculus
    ## 20           ERX1425572    ERR1354097  10090    Mus musculus
    ## 21           ERX1425573    ERR1354098  10090    Mus musculus
    ## 22           ERX1425574    ERR1354099  10090    Mus musculus
    ## 23           ERX1425575    ERR1354100  10090    Mus musculus
    ## 24           ERX1425576    ERR1354101  10090    Mus musculus
    ##       instrument_model library_name library_layout
    ## 1  Illumina HiSeq 2500     15587394         PAIRED
    ## 2  Illumina HiSeq 2500     15587395         PAIRED
    ## 3  Illumina HiSeq 2500     15587396         PAIRED
    ## 4  Illumina HiSeq 2500     15587397         PAIRED
    ## 5  Illumina HiSeq 2500     15587398         PAIRED
    ## 6  Illumina HiSeq 2500     15587399         PAIRED
    ## 7  Illumina HiSeq 2500     15587400         PAIRED
    ## 8  Illumina HiSeq 2500     15587401         PAIRED
    ## 9  Illumina HiSeq 2500     15587402         PAIRED
    ## 10 Illumina HiSeq 2500     15587403         PAIRED
    ## 11 Illumina HiSeq 2500     15587404         PAIRED
    ## 12 Illumina HiSeq 2500     15587405         PAIRED
    ## 13 Illumina HiSeq 2500     15587394         PAIRED
    ## 14 Illumina HiSeq 2500     15587395         PAIRED
    ## 15 Illumina HiSeq 2500     15587396         PAIRED
    ## 16 Illumina HiSeq 2500     15587397         PAIRED
    ## 17 Illumina HiSeq 2500     15587398         PAIRED
    ## 18 Illumina HiSeq 2500     15587399         PAIRED
    ## 19 Illumina HiSeq 2500     15587400         PAIRED
    ## 20 Illumina HiSeq 2500     15587401         PAIRED
    ## 21 Illumina HiSeq 2500     15587402         PAIRED
    ## 22 Illumina HiSeq 2500     15587403         PAIRED
    ## 23 Illumina HiSeq 2500     15587404         PAIRED
    ## 24 Illumina HiSeq 2500     15587405         PAIRED
    ##                                                                                                                                            fastq_ftp
    ## 1  ftp.sra.ebi.ac.uk/vol1/fastq/ERR135/008/ERR1354078/ERR1354078_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/ERR135/008/ERR1354078/ERR1354078_2.fastq.gz
    ## 2  ftp.sra.ebi.ac.uk/vol1/fastq/ERR135/009/ERR1354079/ERR1354079_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/ERR135/009/ERR1354079/ERR1354079_2.fastq.gz
    ## 3  ftp.sra.ebi.ac.uk/vol1/fastq/ERR135/000/ERR1354080/ERR1354080_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/ERR135/000/ERR1354080/ERR1354080_2.fastq.gz
    ## 4  ftp.sra.ebi.ac.uk/vol1/fastq/ERR135/001/ERR1354081/ERR1354081_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/ERR135/001/ERR1354081/ERR1354081_2.fastq.gz
    ## 5  ftp.sra.ebi.ac.uk/vol1/fastq/ERR135/002/ERR1354082/ERR1354082_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/ERR135/002/ERR1354082/ERR1354082_2.fastq.gz
    ## 6  ftp.sra.ebi.ac.uk/vol1/fastq/ERR135/003/ERR1354083/ERR1354083_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/ERR135/003/ERR1354083/ERR1354083_2.fastq.gz
    ## 7  ftp.sra.ebi.ac.uk/vol1/fastq/ERR135/004/ERR1354084/ERR1354084_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/ERR135/004/ERR1354084/ERR1354084_2.fastq.gz
    ## 8  ftp.sra.ebi.ac.uk/vol1/fastq/ERR135/005/ERR1354085/ERR1354085_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/ERR135/005/ERR1354085/ERR1354085_2.fastq.gz
    ## 9  ftp.sra.ebi.ac.uk/vol1/fastq/ERR135/006/ERR1354086/ERR1354086_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/ERR135/006/ERR1354086/ERR1354086_2.fastq.gz
    ## 10 ftp.sra.ebi.ac.uk/vol1/fastq/ERR135/007/ERR1354087/ERR1354087_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/ERR135/007/ERR1354087/ERR1354087_2.fastq.gz
    ## 11 ftp.sra.ebi.ac.uk/vol1/fastq/ERR135/008/ERR1354088/ERR1354088_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/ERR135/008/ERR1354088/ERR1354088_2.fastq.gz
    ## 12 ftp.sra.ebi.ac.uk/vol1/fastq/ERR135/009/ERR1354089/ERR1354089_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/ERR135/009/ERR1354089/ERR1354089_2.fastq.gz
    ## 13 ftp.sra.ebi.ac.uk/vol1/fastq/ERR135/000/ERR1354090/ERR1354090_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/ERR135/000/ERR1354090/ERR1354090_2.fastq.gz
    ## 14 ftp.sra.ebi.ac.uk/vol1/fastq/ERR135/001/ERR1354091/ERR1354091_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/ERR135/001/ERR1354091/ERR1354091_2.fastq.gz
    ## 15 ftp.sra.ebi.ac.uk/vol1/fastq/ERR135/002/ERR1354092/ERR1354092_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/ERR135/002/ERR1354092/ERR1354092_2.fastq.gz
    ## 16 ftp.sra.ebi.ac.uk/vol1/fastq/ERR135/003/ERR1354093/ERR1354093_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/ERR135/003/ERR1354093/ERR1354093_2.fastq.gz
    ## 17 ftp.sra.ebi.ac.uk/vol1/fastq/ERR135/004/ERR1354094/ERR1354094_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/ERR135/004/ERR1354094/ERR1354094_2.fastq.gz
    ## 18 ftp.sra.ebi.ac.uk/vol1/fastq/ERR135/005/ERR1354095/ERR1354095_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/ERR135/005/ERR1354095/ERR1354095_2.fastq.gz
    ## 19 ftp.sra.ebi.ac.uk/vol1/fastq/ERR135/006/ERR1354096/ERR1354096_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/ERR135/006/ERR1354096/ERR1354096_2.fastq.gz
    ## 20 ftp.sra.ebi.ac.uk/vol1/fastq/ERR135/007/ERR1354097/ERR1354097_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/ERR135/007/ERR1354097/ERR1354097_2.fastq.gz
    ## 21 ftp.sra.ebi.ac.uk/vol1/fastq/ERR135/008/ERR1354098/ERR1354098_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/ERR135/008/ERR1354098/ERR1354098_2.fastq.gz
    ## 22 ftp.sra.ebi.ac.uk/vol1/fastq/ERR135/009/ERR1354099/ERR1354099_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/ERR135/009/ERR1354099/ERR1354099_2.fastq.gz
    ## 23 ftp.sra.ebi.ac.uk/vol1/fastq/ERR135/000/ERR1354100/ERR1354100_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/ERR135/000/ERR1354100/ERR1354100_2.fastq.gz
    ## 24 ftp.sra.ebi.ac.uk/vol1/fastq/ERR135/001/ERR1354101/ERR1354101_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/ERR135/001/ERR1354101/ERR1354101_2.fastq.gz
    ##    sra_ftp sample_title
    ## 1       NA       1_ZFP+
    ## 2       NA      1_CD24-
    ## 3       NA      1_CD24+
    ## 4       NA      1_Sca1-
    ## 5       NA       2_ZFP+
    ## 6       NA      2_CD24-
    ## 7       NA      2_CD24+
    ## 8       NA      2_Sca1-
    ## 9       NA       3_ZFP+
    ## 10      NA      3_CD24-
    ## 11      NA      3_CD24+
    ## 12      NA      3_Sca1-
    ## 13      NA       1_ZFP+
    ## 14      NA      1_CD24-
    ## 15      NA      1_CD24+
    ## 16      NA      1_Sca1-
    ## 17      NA       2_ZFP+
    ## 18      NA      2_CD24-
    ## 19      NA      2_CD24+
    ## 20      NA      2_Sca1-
    ## 21      NA       3_ZFP+
    ## 22      NA      3_CD24-
    ## 23      NA      3_CD24+
    ## 24      NA      3_Sca1-

``` r
# dsAgList <-  lapply(dsList,function(x){
#   rnsgs <-  rnSymbolGenes[rnSymbolGenes$ensembl_transcript_id_version %in% rownames(x),]
#   x - x[rnsgs$ensembl_transcript_id_version,]
#   ret <- aggregate(as.integer(x$NumReads), by=list(rnsgs$mgi_symbol),sum)
#   rownames(ret) <- ret[,1]
#   ret[,-1,drop=F]
# })
#
#ambrosiMat <-  as.matrix(Reduce(rn.merge,dsAgList))
txList <-  tximport(fileList,type="salmon",txOut = T)
```

    ## reading in files with read_tsv

    ## 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24

``` r
ambrosiMat <-  summarizeToGene(txList,rnSymbolGenes)$counts[-1,]
```

    ## removing duplicated transcript rows from tx2gene
    ## transcripts missing from tx2gene: 12967
    ## summarizing abundance
    ## summarizing counts
    ## summarizing length

``` r
colnames(ambrosiMat) <- sampleMat$sample_title

#There's a failed sample in there
std.heatmap(cor(ambrosiMat,method="spearman"))
```

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Rowv is FALSE, while dendrogram is `both'. Omitting row
    ## dendogram.

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Colv is FALSE, while dendrogram is `column'. Omitting column
    ## dendogram.

![](BoneNotebook_files/figure-markdown_github/loadData-1.png)

``` r
#Correlation between tech replicates super high. Pool the reads from the technical replicates.
ambrosiMat <- ambrosiMat[,1:12]+ambrosiMat[,13:24]
#Get rid of failed sample
ambrosiMat <- ambrosiMat[,-which(colnames(ambrosiMat)=="2_CD24-")]
#Check failed sample removal
std.heatmap(cor(ambrosiMat,method = "spearman"))
```

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Rowv is FALSE, while dendrogram is `both'. Omitting row
    ## dendogram.

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Colv is FALSE, while dendrogram is `column'. Omitting column
    ## dendogram.

![](BoneNotebook_files/figure-markdown_github/loadData-2.png)

``` r
ambrosiMatNorm <-  median.normalize(ambrosiMat[rowMaxs(ambrosiMat)>2,])
condits <- sapply(strsplit(colnames(ambrosiMat),"_"),function(x)x[2])
ambrosiMat <- ambrosiMat[,order(condits)]
ambrosiMatNorm <- ambrosiMatNorm[,order(condits)]
condits <- condits[order(condits)]
condits <- gsub("\\+","plus",condits)
condits <- gsub("\\-","minus",condits)
```

Now we have normalized counts in "ambrosiMatNorm", the conditions in "condits", and the unnormalized counts for differential expression in "ambrosiMat."

### A few genes of interest

``` r
head(ambrosiMat)
```

    ##                 1_CD24-   3_CD24-    1_CD24+    2_CD24+    3_CD24+
    ## 0610009B22Rik 168.00000 291.00000  240.00000 499.000000  360.00000
    ## 0610009O20Rik 499.00000 784.00000 1293.00000 467.000000 1054.00000
    ## 0610010F05Rik   0.00000 379.00000    7.00000 132.000000    0.00000
    ## 0610010K14Rik 697.63399 117.73001  273.36182 561.736115  226.69839
    ## 0610012G03Rik  12.02473  27.97399   77.97179  39.534603   22.93296
    ## 0610030E20Rik 181.16982  83.18100   60.10451   5.043553  565.96365
    ##                 1_Sca1-    2_Sca1-  3_Sca1-    1_ZFP+   2_ZFP+    3_ZFP+
    ## 0610009B22Rik 174.00000 347.000000 301.5604 441.00000 437.0000 562.00000
    ## 0610009O20Rik 575.00000 663.000000 219.0000 520.00000 982.0000 744.00000
    ## 0610010F05Rik 284.00000 209.000001 107.0000 216.00000  11.0000   7.00000
    ## 0610010K14Rik 177.14379  99.160890 151.7189 598.96206 967.5398 448.18364
    ## 0610012G03Rik  69.56521   2.412995   0.0000  64.76743  71.5059  28.75185
    ## 0610030E20Rik  34.00000 873.346671  19.0000 368.30748 390.3008 438.34378

``` r
head(ambrosiMatNorm)
```

    ##                 1_CD24-   3_CD24-    1_CD24+    2_CD24+    3_CD24+
    ## 0610009B22Rik 155.31431 352.25089  212.88925 397.418205  347.91985
    ## 0610009O20Rik 461.32047 949.01959 1146.94084 371.932468 1018.63200
    ## 0610010F05Rik   0.00000 458.77350    6.20927 105.128663    0.00000
    ## 0610010K14Rik 644.95558 142.51031  242.48247 447.383083  219.09131
    ## 0610012G03Rik  11.11674  33.86208   69.16398  31.486515   22.16342
    ## 0610030E20Rik 167.48967 100.68929   53.31502   4.016833  546.97218
    ##                 1_Sca1-    2_Sca1-   3_Sca1-    1_ZFP+     2_ZFP+
    ## 0610009B22Rik 235.30292 329.986313 393.85481 303.75905 297.504949
    ## 0610009O20Rik 777.58148 630.492581 286.02629 358.17393 668.535147
    ## 0610010F05Rik 384.05764 198.752564 139.74801 148.77994   7.488683
    ## 0610010K14Rik 239.55431  94.298952 198.15342 412.56269 658.690778
    ## 0610012G03Rik  94.07412   2.294684   0.00000  44.61155  48.680458
    ## 0610030E20Rik  45.97873 830.525787  24.81507 253.68873 265.712629
    ##                   3_ZFP+
    ## 0610009B22Rik 426.661806
    ## 0610009O20Rik 564.833423
    ## 0610010F05Rik   5.314293
    ## 0610010K14Rik 340.254163
    ## 0610012G03Rik  21.827964
    ## 0610030E20Rik 332.783894

``` r
print(condits)
```

    ##  [1] "CD24minus" "CD24minus" "CD24plus"  "CD24plus"  "CD24plus" 
    ##  [6] "Sca1minus" "Sca1minus" "Sca1minus" "ZFPplus"   "ZFPplus"  
    ## [11] "ZFPplus"

``` r
barplot(ambrosiMatNorm["Esr1",],las=2,main="Esr1")
```

![](BoneNotebook_files/figure-markdown_github/checkData-1.png)

``` r
barplot(ambrosiMatNorm["Esr2",],las=2,main="Esr2")
```

![](BoneNotebook_files/figure-markdown_github/checkData-2.png)

``` r
barplot(ambrosiMatNorm["Gper1",],las=2,main="Gper1")
```

![](BoneNotebook_files/figure-markdown_github/checkData-3.png)

### PCA

Let's redo the Principal component analysis (or singular value decomposition, svd) to check against figure 5 of the paper.

``` r
ambrosiCLN <- round.log(ambrosiMatNorm+1,2)
ambrosiCLN <- ambrosiCLN[rowSds(ambrosiCLN)>1,]
svAmbrosi <- svd((ambrosiCLN-rowMeans(ambrosiCLN))/rowSds(ambrosiCLN))

std.heatmap(cor(ambrosiMatNorm,method = "spearman"))
```

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Rowv is FALSE, while dendrogram is `both'. Omitting row
    ## dendogram.

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Colv is FALSE, while dendrogram is `column'. Omitting column
    ## dendogram.

![](BoneNotebook_files/figure-markdown_github/svd-1.png)

``` r
conditNums <- sapply(condits,function(x)which(x==sort(unique(condits))))
#Sca1minus is osteo, ZFP (mature) and CD24- (less mature) are adipocytes, cd24+ is multipotent
plot(svAmbrosi$v[,1:2],col=conditNums,xlab="PC1",ylab="PC2")
legend("bottomright",legend = unique(names(conditNums)),col=1:4,fill = 1:4)
```

![](BoneNotebook_files/figure-markdown_github/svd-2.png)

``` r
eaciout <- list()
l <- 1
eacivector <- svAmbrosi$u[,l]
names(eacivector) <- rownames(ambrosiCLN)
eaciout[[l]] <- eacitest(eacivector,"org.Mm.eg","SYMBOL",sets = "GO")$setscores
```

    ## Loading necessary libraries...

    ## Loaded Package org.Mm.eg.db

    ## Converting annotations to data.frames ...

    ## iteration 1 done; time  9.51 sec 
    ## iteration 2 done; time  7.51 sec 
    ## iteration 3 done; time  7.78 sec 
    ## iteration 4 done; time  8.01 sec 
    ## iteration 5 done; time  7.24 sec 
    ## iteration 6 done; time  7.81 sec 
    ## iteration 7 done; time  7.58 sec 
    ## iteration 8 done; time  7.77 sec 
    ## iteration 9 done; time  7.83 sec 
    ## iteration 10 done; time  8.03 sec

    ## Labeling output ...

    ## Loaded Package GO.db

``` r
l <- 2
eacivector <- svAmbrosi$u[,l]
names(eacivector) <- rownames(ambrosiCLN)
eaciout[[l]] <- eacitest(eacivector,"org.Mm.eg","SYMBOL",sets = "GO")$setscores
```

    ## Loading necessary libraries...

    ## Loaded Package org.Mm.eg.db

    ## Converting annotations to data.frames ...

    ## iteration 1 done; time  5.94 sec 
    ## iteration 2 done; time  7.86 sec 
    ## iteration 3 done; time  9.61 sec 
    ## iteration 4 done; time  7.36 sec 
    ## iteration 5 done; time  7.45 sec 
    ## iteration 6 done; time  7.24 sec 
    ## iteration 7 done; time  5.73 sec 
    ## iteration 8 done; time  8.26 sec 
    ## iteration 9 done; time  9.45 sec 
    ## iteration 10 done; time  8.23 sec

    ## Labeling output ...

    ## Loaded Package GO.db

### GO enrichment of Principal Components

So that reiterates the Ambrosi analysis. Now let's look at the continuous GO enrichment in the genes that contribute to PC1 (separates the osteocyte and progenitors from adipocytes)

``` r
#PC1 Positive
a <- eaciout[[1]][eaciout[[1]]$set.mean>0,]
print(a[1:25,])
```

    ##                                                               Term
    ## GO:1904181          positive regulation of membrane depolarization
    ## GO:0010935            regulation of macrophage cytokine production
    ## GO:0035455                            response to interferon-alpha
    ## GO:0035641                         locomotory exploration behavior
    ## GO:0070006                          metalloaminopeptidase activity
    ## GO:0043034                                               costamere
    ## GO:0005161         platelet-derived growth factor receptor binding
    ## GO:0030742                           GTP-dependent protein binding
    ## GO:0043649                     dicarboxylic acid catabolic process
    ## GO:0051497            negative regulation of stress fiber assembly
    ## GO:0097440                                         apical dendrite
    ## GO:0038191                                      neuropilin binding
    ## GO:0060716             labyrinthine layer blood vessel development
    ## GO:0034312                               diol biosynthetic process
    ## GO:0052744 phosphatidylinositol monophosphate phosphatase activity
    ## GO:0097320                              plasma membrane tubulation
    ## GO:2001212                            regulation of vasculogenesis
    ## GO:0000062                                  fatty-acyl-CoA binding
    ## GO:0042581                                        specific granule
    ## GO:1990126       retrograde transport, endosome to plasma membrane
    ## GO:0051654             establishment of mitochondrion localization
    ## GO:0031579                              membrane raft organization
    ## GO:0004602                         glutathione peroxidase activity
    ## GO:0030574                              collagen catabolic process
    ## GO:0021846                         cell proliferation in forebrain
    ##            Ontology    set.mean      set.sd set.size         pval
    ## GO:1904181       BP 0.009550961 0.007362549       11 0.000000e+00
    ## GO:0010935       BP 0.008480568 0.004542898       10 0.000000e+00
    ## GO:0035455       BP 0.008176674 0.004452738       16 0.000000e+00
    ## GO:0035641       BP 0.008128965 0.003041400       12 0.000000e+00
    ## GO:0070006       MF 0.008062144 0.003396604        9 0.000000e+00
    ## GO:0043034       CC 0.007416668 0.003652416       12 0.000000e+00
    ## GO:0005161       MF 0.007303739 0.004115494       10 0.000000e+00
    ## GO:0030742       MF 0.007236354 0.005686460       16 0.000000e+00
    ## GO:0043649       BP 0.007048522 0.004348495        9 6.661338e-16
    ## GO:0051497       BP 0.006929356 0.004769006       19 1.998401e-15
    ## GO:0097440       CC 0.006927306 0.008652274       17 1.998401e-15
    ## GO:0038191       MF 0.006925863 0.009478934       11 1.998401e-15
    ## GO:0060716       BP 0.006921120 0.003255819       13 1.998401e-15
    ## GO:0034312       BP 0.006855821 0.007104589       10 3.774758e-15
    ## GO:0052744       MF 0.006738662 0.002955610        9 1.110223e-14
    ## GO:0097320       BP 0.006549348 0.004641722       11 5.995204e-14
    ## GO:2001212       BP 0.006467521 0.007580034       11 1.225686e-13
    ## GO:0000062       MF 0.006281751 0.003783333       10 5.995204e-13
    ## GO:0042581       CC 0.006236637 0.008351424       11 8.757439e-13
    ## GO:1990126       BP 0.006187573 0.004099407       10 1.318279e-12
    ## GO:0051654       BP 0.006159705 0.003921547       10 1.660672e-12
    ## GO:0031579       BP 0.006144897 0.006097797        9 1.876721e-12
    ## GO:0004602       MF 0.006067209 0.005298084        9 3.548717e-12
    ## GO:0030574       BP 0.005944741 0.003864840       10 9.534373e-12
    ## GO:0021846       BP 0.005903418 0.007465946       18 1.325007e-11

``` r
#PC1 Negative
a <- eaciout[[1]][eaciout[[1]]$set.mean<0,]
print(a[1:25,])
```

    ##                                                                              Term
    ## GO:0003417                                     growth plate cartilage development
    ## GO:0042555                                                            MCM complex
    ## GO:0006271                      DNA strand elongation involved in DNA replication
    ## GO:0050699                                                      WW domain binding
    ## GO:0000940                                 condensed chromosome outer kinetochore
    ## GO:0031643                                     positive regulation of myelination
    ## GO:0003688                                         DNA replication origin binding
    ## GO:0043142                          single-stranded DNA-dependent ATPase activity
    ## GO:0051084                            'de novo' posttranslational protein folding
    ## GO:0030206                               chondroitin sulfate biosynthetic process
    ## GO:1904666                        regulation of ubiquitin protein ligase activity
    ## GO:0050911 detection of chemical stimulus involved in sensory perception of smell
    ## GO:0000800                                                        lateral element
    ## GO:0005861                                                       troponin complex
    ## GO:0002076                                                 osteoblast development
    ## GO:0061436                                          establishment of skin barrier
    ## GO:0000796                                                      condensin complex
    ## GO:0006198                                                 cAMP catabolic process
    ## GO:0010369                                                           chromocenter
    ## GO:0034501                                    protein localization to kinetochore
    ## GO:0018279                          protein N-linked glycosylation via asparagine
    ## GO:0005251                           delayed rectifier potassium channel activity
    ## GO:0004936                                     alpha-adrenergic receptor activity
    ## GO:0030033                                                   microvillus assembly
    ## GO:0045120                                                             pronucleus
    ##            Ontology     set.mean      set.sd set.size         pval
    ## GO:0003417       BP -0.007913482 0.006350914        9 8.242433e-20
    ## GO:0042555       CC -0.007165374 0.005537788        9 1.579214e-16
    ## GO:0006271       BP -0.006697905 0.003346225       10 1.226746e-14
    ## GO:0050699       MF -0.006399908 0.003496455       16 1.695746e-13
    ## GO:0000940       CC -0.006360513 0.004464935       12 2.379024e-13
    ## GO:0031643       BP -0.005818542 0.003866918       12 2.044247e-11
    ## GO:0003688       MF -0.005554578 0.003062234       11 1.558571e-10
    ## GO:0043142       MF -0.005363276 0.003052580       14 6.421543e-10
    ## GO:0051084       BP -0.005140655 0.007101077        9 3.143476e-09
    ## GO:0030206       BP -0.004989020 0.006338731        7 8.940523e-09
    ## GO:1904666       BP -0.004803296 0.007399239       11 3.089380e-08
    ## GO:0050911       BP -0.004620843 0.002097693       13 1.000393e-07
    ## GO:0000800       CC -0.004529402 0.006367887       14 1.773971e-07
    ## GO:0005861       CC -0.004385660 0.004648931        7 4.271912e-07
    ## GO:0002076       BP -0.004371176 0.006836826       10 4.660641e-07
    ## GO:0061436       BP -0.004240701 0.003095505       15 1.009057e-06
    ## GO:0000796       CC -0.004222971 0.003853109        7 1.118854e-06
    ## GO:0006198       BP -0.004200151 0.003667048       11 1.277187e-06
    ## GO:0010369       CC -0.004170686 0.002705810       10 1.513719e-06
    ## GO:0034501       BP -0.004121958 0.002540034       17 1.999980e-06
    ## GO:0018279       BP -0.003973351 0.005387521       15 4.590650e-06
    ## GO:0005251       MF -0.003906662 0.004965269       11 6.604611e-06
    ## GO:0004936       MF -0.003844353 0.004192067        4 9.230585e-06
    ## GO:0030033       BP -0.003498827 0.009998889       10 5.403767e-05
    ## GO:0045120       CC -0.003392271 0.007256207       13 9.040247e-05

Interesting... Now PC2 (separates osteocytes from progenitors)

``` r
#PC2 Positive
a <- eaciout[[2]][eaciout[[2]]$set.mean>0,]
print(a[1:25,])
```

    ##                                                                           Term
    ## GO:0000076                                          DNA replication checkpoint
    ## GO:0042168                                              heme metabolic process
    ## GO:0019825                                                      oxygen binding
    ## GO:0042555                                                         MCM complex
    ## GO:0043034                                                           costamere
    ## GO:0060004                                                              reflex
    ## GO:0035641                                     locomotory exploration behavior
    ## GO:0048821                                             erythrocyte development
    ## GO:0052744             phosphatidylinositol monophosphate phosphatase activity
    ## GO:0000940                              condensed chromosome outer kinetochore
    ## GO:0006271                   DNA strand elongation involved in DNA replication
    ## GO:0007076                                     mitotic chromosome condensation
    ## GO:0002098                                    tRNA wobble uridine modification
    ## GO:0043142                       single-stranded DNA-dependent ATPase activity
    ## GO:0030539                                          male genitalia development
    ## GO:0043567 regulation of insulin-like growth factor receptor signaling pathway
    ## GO:0008139                               nuclear localization sequence binding
    ## GO:0007064                                   mitotic sister chromatid cohesion
    ## GO:0034508                                         centromere complex assembly
    ## GO:0005847        mRNA cleavage and polyadenylation specificity factor complex
    ## GO:0071108                                 protein K48-linked deubiquitination
    ## GO:0003688                                      DNA replication origin binding
    ## GO:0042588                                                     zymogen granule
    ## GO:0032769                       negative regulation of monooxygenase activity
    ## GO:0043240                                     Fanconi anaemia nuclear complex
    ##            Ontology    set.mean      set.sd set.size         pval
    ## GO:0000076       BP 0.009693135 0.004932242       10 0.000000e+00
    ## GO:0042168       BP 0.007914790 0.009950480       14 0.000000e+00
    ## GO:0019825       MF 0.007654870 0.004611019       13 0.000000e+00
    ## GO:0042555       CC 0.007441273 0.004535643        9 0.000000e+00
    ## GO:0043034       CC 0.006902043 0.003410221       12 0.000000e+00
    ## GO:0060004       BP 0.006565326 0.003320756       13 2.220446e-16
    ## GO:0035641       BP 0.006333117 0.002822265       12 1.776357e-15
    ## GO:0048821       BP 0.005632202 0.010486405       21 1.544986e-12
    ## GO:0052744       MF 0.005609784 0.003603374        9 1.891154e-12
    ## GO:0000940       CC 0.005590678 0.004938136       12 2.245759e-12
    ## GO:0006271       BP 0.005444022 0.002063442       10 8.237189e-12
    ## GO:0007076       BP 0.005298614 0.006001636        9 2.891776e-11
    ## GO:0002098       BP 0.005133174 0.003917011        9 1.160048e-10
    ## GO:0043142       MF 0.004823193 0.002418823       14 1.398381e-09
    ## GO:0030539       BP 0.004799560 0.003085027       12 1.680430e-09
    ## GO:0043567       BP 0.004590763 0.002749857       16 8.207849e-09
    ## GO:0008139       MF 0.004501030 0.002706107       10 1.589782e-08
    ## GO:0007064       BP 0.004474369 0.004951665       13 1.930226e-08
    ## GO:0034508       BP 0.004467597 0.003425169       17 2.027394e-08
    ## GO:0005847       CC 0.004274121 0.002753735        8 8.006887e-08
    ## GO:0071108       BP 0.004062983 0.002363101       11 3.358228e-07
    ## GO:0003688       MF 0.004058608 0.002412097       11 3.456995e-07
    ## GO:0042588       CC 0.004049460 0.002115837       12 3.672673e-07
    ## GO:0032769       BP 0.004039837 0.002196960       12 3.913548e-07
    ## GO:0043240       CC 0.003699106 0.006487035       12 3.388582e-06

``` r
#PC2 Negative
a <- eaciout[[2]][eaciout[[2]]$set.mean<0,]
print(a[1:25,])
```

    ##                                                                  Term
    ## GO:0042608                                    T cell receptor binding
    ## GO:0004745                             retinol dehydrogenase activity
    ## GO:0005779                 integral component of peroxisomal membrane
    ## GO:0070402                                              NADPH binding
    ## GO:0016755       transferase activity, transferring amino-acyl groups
    ## GO:0045948            positive regulation of translational initiation
    ## GO:0030687                       preribosome, large subunit precursor
    ## GO:0048845                          venous blood vessel morphogenesis
    ## GO:1904424                                  regulation of GTP binding
    ## GO:0048875                       chemical homeostasis within a tissue
    ## GO:0032823          regulation of natural killer cell differentiation
    ## GO:0035859                                    Seh1-associated complex
    ## GO:0006677                         glycosylceramide metabolic process
    ## GO:0090502        RNA phosphodiester bond hydrolysis, endonucleolytic
    ## GO:0004467                  long-chain fatty acid-CoA ligase activity
    ## GO:0008356                                   asymmetric cell division
    ## GO:0009931 calcium-dependent protein serine/threonine kinase activity
    ## GO:0031902                                     late endosome membrane
    ## GO:0007214                  gamma-aminobutyric acid signaling pathway
    ## GO:0004115                3',5'-cyclic-AMP phosphodiesterase activity
    ## GO:0046966                           thyroid hormone receptor binding
    ## GO:0097546                                               ciliary base
    ## GO:0036158                                  outer dynein arm assembly
    ## GO:0045236                            CXCR chemokine receptor binding
    ## GO:0050774              negative regulation of dendrite morphogenesis
    ##            Ontology     set.mean      set.sd set.size         pval
    ## GO:0042608       MF -0.007225827 0.002514962       12 1.338052e-19
    ## GO:0004745       MF -0.006967215 0.004736985       10 2.483191e-18
    ## GO:0005779       CC -0.006921376 0.003991866       11 4.122337e-18
    ## GO:0070402       MF -0.006843526 0.003357262       11 9.677517e-18
    ## GO:0016755       MF -0.006644089 0.006152538       13 8.252173e-17
    ## GO:0045948       BP -0.006520814 0.006213942       12 3.009541e-16
    ## GO:0030687       CC -0.006427968 0.002866980       13 7.851767e-16
    ## GO:0048845       BP -0.006133401 0.003886870       10 1.506047e-14
    ## GO:1904424       BP -0.006115307 0.002902614       10 1.797771e-14
    ## GO:0048875       BP -0.006060489 0.006105307       11 3.064503e-14
    ## GO:0032823       BP -0.006047297 0.006936837       12 3.481768e-14
    ## GO:0035859       CC -0.006046898 0.008022817       11 3.495241e-14
    ## GO:0006677       BP -0.005997180 0.005862817       11 5.640960e-14
    ## GO:0090502       BP -0.005861337 0.004638792       12 2.045740e-13
    ## GO:0004467       MF -0.005442840 0.004942386        8 9.049884e-12
    ## GO:0008356       BP -0.005246803 0.004761979       11 4.866223e-11
    ## GO:0009931       MF -0.005184949 0.006868003        8 8.172605e-11
    ## GO:0031902       CC -0.004992897 0.003392507       13 3.937092e-10
    ## GO:0007214       BP -0.004960321 0.002685419       12 5.111538e-10
    ## GO:0004115       MF -0.004908306 0.003028564       11 7.728744e-10
    ## GO:0046966       MF -0.004798493 0.005462973       12 1.824991e-09
    ## GO:0097546       CC -0.004753170 0.005377825       19 2.587775e-09
    ## GO:0036158       BP -0.004745586 0.003401511        9 2.742657e-09
    ## GO:0045236       MF -0.004722545 0.004456333       10 3.270645e-09
    ## GO:0050774       BP -0.004630942 0.002674785       10 6.533070e-09

Sweet.

Time for some differential expression

``` r
cond <- as.factor(condits)
colnames(ambrosiMat) <- make.names(condits,unique = T)
dds <- DESeqDataSetFromMatrix(round(ambrosiMat),colData = DataFrame(cond),design = formula(~cond+0))
```

    ## converting counts to integer mode

``` r
DESeqOutput <-  DESeq(dds)
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

``` r
resultsNames(DESeqOutput)
```

    ## [1] "condCD24minus" "condCD24plus"  "condSca1minus" "condZFPplus"

``` r
resList <-  list(results(DESeqOutput,contrast = c(1,-1/3,-1/3,-1/3 ),cooksCutoff=T),results(DESeqOutput,contrast = c(-1/3,1,-1/3,-1/3 ),cooksCutoff=T),results(DESeqOutput,contrast = c(-1/3,-1/3,1,-1/3 ),cooksCutoff=T),results(DESeqOutput,contrast = c(-1/3,-1/3,-1/3 ,1),cooksCutoff=T))

ambrosiUpDown <-lapply(resList,function(res){
  res <- res[!is.na(res$padj),]
  list(rownames(res[res$padj<.1&res$log2FoldChange>0,]),rownames(res[res$padj<.1&res$log2FoldChange<0,]))
})

print(str(ambrosiUpDown))
```

    ## List of 4
    ##  $ :List of 2
    ##   ..$ : chr [1:419] "1700066M21Rik" "5430403G16Rik" "6430548M08Rik" "8430408G22Rik" ...
    ##   ..$ : chr [1:591] "1110008L16Rik" "1810041L15Rik" "2610008E11Rik" "3110062M04Rik" ...
    ##  $ :List of 2
    ##   ..$ : chr [1:896] "1190002N15Rik" "1700028J19Rik" "1700047I17Rik2" "1810013L24Rik" ...
    ##   ..$ : chr [1:800] "1700019D03Rik" "2200002D01Rik" "2310009B15Rik" "2310061I04Rik" ...
    ##  $ :List of 2
    ##   ..$ : chr [1:894] "1500015O10Rik" "1700008O03Rik" "1700047I17Rik2" "1810041L15Rik" ...
    ##   ..$ : chr [1:853] "0610040J01Rik" "1700029J07Rik" "1810043G02Rik" "2510039O18Rik" ...
    ##  $ :List of 2
    ##   ..$ : chr [1:689] "2200002D01Rik" "2810021J22Rik" "3830406C13Rik" "4930432E11Rik" ...
    ##   ..$ : chr [1:765] "1500015O10Rik" "1700028J19Rik" "1700047I17Rik2" "1700066M21Rik" ...
    ## NULL

So there's a bunch of genes differentially expressed up or down for each cell type compared to all the others.

``` r
for(i in 1:length(resList)){
  print(resultsNames(DESeqOutput)[i])
  res <- resList[[i]]
  res <- res[!is.na(res$padj),]
  write.table(res[order(res$log2FoldChange,decreasing = T),],paste0("~/Desktop/DEG_",resultsNames(DESeqOutput)[i],".txt"),sep = "\t",quote = F)
  res <- res[res$padj<.1&res$log2FoldChange>0,]
  mat <- log(ambrosiMatNorm+1,2)-rowMeans(log(ambrosiMatNorm+1,2))
  std.heatmap(mat[rownames(res[order(res$padj,decreasing = F),])[1:25],],main=paste(resultsNames(DESeqOutput)[i],"Up vs All\nLogFC vs mean"))
}
```

    ## [1] "condCD24minus"

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Rowv is FALSE, while dendrogram is `both'. Omitting row
    ## dendogram.

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Colv is FALSE, while dendrogram is `column'. Omitting column
    ## dendogram.

    ## [1] "condCD24plus"

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Rowv is FALSE, while dendrogram is `both'. Omitting row
    ## dendogram.

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Colv is FALSE, while dendrogram is `column'. Omitting column
    ## dendogram.

![](BoneNotebook_files/figure-markdown_github/hmDEs-1.png)

    ## [1] "condSca1minus"

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Rowv is FALSE, while dendrogram is `both'. Omitting row
    ## dendogram.

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Colv is FALSE, while dendrogram is `column'. Omitting column
    ## dendogram.

![](BoneNotebook_files/figure-markdown_github/hmDEs-2.png)

    ## [1] "condZFPplus"

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Rowv is FALSE, while dendrogram is `both'. Omitting row
    ## dendogram.

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Colv is FALSE, while dendrogram is `column'. Omitting column
    ## dendogram.

![](BoneNotebook_files/figure-markdown_github/hmDEs-3.png)![](BoneNotebook_files/figure-markdown_github/hmDEs-4.png)

``` r
for(i in 1:length(resList)){
  print(resultsNames(DESeqOutput)[i])
  res <- resList[[i]]
  res <- res[!is.na(res$padj),]
  res <- res[res$padj<.1&res$log2FoldChange>0,]
  std.heatmap(log(ambrosiMatNorm[rownames(res[order(res$padj,decreasing = F),])[1:25],]+1,2),main=paste(resultsNames(DESeqOutput)[i],"Up vs All"))
}
```

    ## [1] "condCD24minus"

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Rowv is FALSE, while dendrogram is `both'. Omitting row
    ## dendogram.

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Colv is FALSE, while dendrogram is `column'. Omitting column
    ## dendogram.

    ## [1] "condCD24plus"

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Rowv is FALSE, while dendrogram is `both'. Omitting row
    ## dendogram.

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Colv is FALSE, while dendrogram is `column'. Omitting column
    ## dendogram.

![](BoneNotebook_files/figure-markdown_github/hmDEs-5.png)

    ## [1] "condSca1minus"

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Rowv is FALSE, while dendrogram is `both'. Omitting row
    ## dendogram.

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Colv is FALSE, while dendrogram is `column'. Omitting column
    ## dendogram.

![](BoneNotebook_files/figure-markdown_github/hmDEs-6.png)

    ## [1] "condZFPplus"

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Rowv is FALSE, while dendrogram is `both'. Omitting row
    ## dendogram.

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Colv is FALSE, while dendrogram is `column'. Omitting column
    ## dendogram.

![](BoneNotebook_files/figure-markdown_github/hmDEs-7.png)![](BoneNotebook_files/figure-markdown_github/hmDEs-8.png)

Candice's Data
--------------

The data for the Ingraham lab RNAseq was also passed to Salmon after fastqc and trimming with trimgalore as suggested by the NuGen Prep.

``` r
datapath <- "~/code/IngrahamLabData/BoneSalmonOutputs/"
fileList <- dir(datapath)
fileList <- fileList[!grepl("Gene|pdf",fileList)]
dsList <- lapply(paste0(datapath,fileList),read.csv2, sep="\t",header=T,row.names=1,stringsAsFactors=F)
allRownames <- Reduce(union,lapply(dsList,rownames))

mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host="www.ensembl.org")
rnSymbol <- getBM(attributes = c("ensembl_transcript_id_version","mgi_symbol",'description'),filters = c("ensembl_transcript_id_version"),values =allRownames ,mart = mart) 
rnSymbolGenes <- rnSymbol[rnSymbol$mgi_symbol!=""& !grepl("predicted gene", rnSymbol$description),]
```

``` r
# dsAgList <-  lapply(dsList,function(x){
#   rnsgs <-  rnSymbolGenes[rnSymbolGenes$ensembl_transcript_id_version %in% rownames(x),]
#   x <- x[rnsgs$ensembl_transcript_id_version,]
#   ret <- aggregate(as.integer(x$NumReads), by=list(rnsgs$mgi_symbol),sum)
#   rownames(ret) <- ret[,1]
#   ret[,-1,drop=F]
# })

dsAgList <- lapply(colnames(dsList[[1]]),function(n){
  Reduce(cbind,lapply(dsList,function(x)x[,n,drop=F]))
  })
names(dsAgList) <- names(txList)[c(3,4,1,2)]
dsAgList <- lapply(dsAgList,function(x){
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  x
  })
boneMat <- summarizeToGene(dsAgList,rnSymbolGenes)$counts[-1,]
```

    ## removing duplicated transcript rows from tx2gene

    ## transcripts missing from tx2gene: 13832

    ## summarizing abundance

    ## summarizing counts

    ## summarizing length

``` r
sampleNames <- c("1807_BM_fl_A1","1810_BM_KO_E1","1811_BM_KO_G1","1815_BM_fl_B1","1818_BM_fl_C1","1825_BM_KO_F1","1984_BM_fl_D1","1985_BM_KO_H1")
SampleNameMat <- sapply(strsplit(sampleNames,"_"),function(i)i)
#boneMat <-  as.matrix(Reduce(rn.merge,dsAgList))
colnames(boneMat) <- paste(SampleNameMat[3,],gsub("[[:digit:]]","",SampleNameMat[4,]),sep = "_")
boneMatNorm <-  median.normalize(boneMat)
boneMatNorm <- boneMatNorm[,order(colnames(boneMatNorm))]
heatmap.2(cor(boneMatNorm,method = "spe"),col=cols,trace="none")
```

![](BoneNotebook_files/figure-markdown_github/BiomartSeparate-1.png)

``` r
std.heatmap(cor(rn.merge(boneMatNorm,ambrosiMatNorm),method = "spe"),main="Spearman Correlation\n Ambrosi vs Candice")
```

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Rowv is FALSE, while dendrogram is `both'. Omitting row
    ## dendogram.

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Colv is FALSE, while dendrogram is `column'. Omitting column
    ## dendogram.

![](BoneNotebook_files/figure-markdown_github/BiomartSeparate-2.png)

### Differential expression of Candice's data

Use a DESeq2 False Discovery Rate of .1, breaking into up and down in KO groups.

``` r
cond <- as.factor(SampleNameMat[3,])
dds <- DESeqDataSetFromMatrix(round(boneMat),colData = DataFrame(cond),design = ~cond)
```

    ## converting counts to integer mode

``` r
DESeqOutput <-  DESeq(dds)
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

``` r
resCandice <-  results(DESeqOutput)
res <- resCandice[!is.na(resCandice$padj),]
write.table(res[order(res$log2FoldChange,decreasing = F),],paste0("~/Desktop/DEG","Ingraham",".txt"),sep = "\t",quote=F)
res <- res[res$log2FoldChange<0,]


std.heatmap(log(boneMatNorm[rownames(res[order(res$pvalue,decreasing = F),])[1:25],]+1,2),main="Most significant DE genes\ndown in KO\nlog2(NormalizedCounts+1)")
```

![](BoneNotebook_files/figure-markdown_github/differentialExpression-1.png)

``` r
std.heatmap(log(boneMatNorm[rownames(res[order(res$pvalue,decreasing = F),])[1:25],]+1,2)-rowMeans(log(boneMatNorm[rownames(res[order(res$pvalue,decreasing = F),])[1:25],]+1,2)),main="Most significant DE genes\ndown in KO\nlog2(FC)")
```

![](BoneNotebook_files/figure-markdown_github/differentialExpression-2.png)

``` r
res <-  results(DESeqOutput)
res <- res[!is.na(res$padj),]
boneUpDown <- list(rownames(res[res$padj<.1&res$log2FoldChange>0,]),rownames(res[res$padj<.1&res$log2FoldChange<0,]))
res <- res[res$log2FoldChange>0,]
std.heatmap(log(boneMatNorm[rownames(res[order(res$pvalue,decreasing = F),])[1:25],]+1,2)-rowMeans(log(boneMatNorm[rownames(res[order(res$pvalue,decreasing = F),])[1:25],]+1,2)),main="Most significant DE genes\nup in KO\nlog2(FC)")
```

![](BoneNotebook_files/figure-markdown_github/differentialExpression-3.png)

``` r
std.heatmap(log(boneMatNorm[rownames(res[order(res$pvalue,decreasing = F),])[1:25],]+1,2),main="Most significant DE genes\nup in KO\nlog2(normalized counts + 1)")
```

![](BoneNotebook_files/figure-markdown_github/differentialExpression-4.png)

``` r
res <-  results(DESeqOutput)
res <- res[!is.na(res$padj),]


hist(res$log2FoldChange,main = "Log2 Fold Changes Detected",breaks = 40)
```

![](BoneNotebook_files/figure-markdown_github/differentialExpression-5.png)

``` r
plot(res$log2FoldChange,-log(res$padj),ylab="-logPadj",xlab="logFC",main="Volcano Plot")
```

![](BoneNotebook_files/figure-markdown_github/differentialExpression-6.png)

``` r
#ESR1 not differentially expressed
barplot((boneMatNorm["Esr1",]),main="Esr1")
```

![](BoneNotebook_files/figure-markdown_github/differentialExpression-7.png)

``` r
barplot((boneMatNorm["Gper1",]),main="Gper1")
```

![](BoneNotebook_files/figure-markdown_github/differentialExpression-8.png)

``` r
#Just a sanity Check
barplot(boneMatNorm["Ncoa1",],las=2,main="Ncoa1")
```

![](BoneNotebook_files/figure-markdown_github/differentialExpression-9.png)

``` r
barplot(boneMatNorm["Ncoa2",],las=2,main="Ncoa2")
```

![](BoneNotebook_files/figure-markdown_github/differentialExpression-10.png)

``` r
barplot(boneMatNorm["Ncoa3",],las=2,main="Ncoa3")
```

![](BoneNotebook_files/figure-markdown_github/differentialExpression-11.png)

``` r
barplot((boneMatNorm["Kiss1",]),main="Kiss1")
```

![](BoneNotebook_files/figure-markdown_github/differentialExpression-12.png)

``` r
#also a sanity check
std.heatmap(cor(ambrosiMatNorm,method = "spearman"))
```

![](BoneNotebook_files/figure-markdown_github/differentialExpression-13.png)

``` r
entrezmm <- getBM(attributes = c('mgi_symbol', "entrezgene"), filters = "mgi_symbol",
                   values = rownames(res[res$padj<.1,]), mart = mart)
eKG <- enrichKEGG(entrezmm$entrezgene[!is.na(entrezmm$entrezgene)],organism = "mmu",pvalueCutoff = .1)
eKG
```

    ## #
    ## # over-representation test
    ## #
    ## #...@organism     mmu 
    ## #...@ontology     KEGG 
    ## #...@keytype      kegg 
    ## #...@gene     chr [1:258] "67851" "668661" "224904" "232345" "11370" "11504" ...
    ## #...pvalues adjusted by 'BH' with cutoff <0.1 
    ## #...9 enriched terms found
    ## 'data.frame':    9 obs. of  9 variables:
    ##  $ ID         : chr  "mmu03010" "mmu00190" "mmu05012" "mmu04714" ...
    ##  $ Description: chr  "Ribosome" "Oxidative phosphorylation" "Parkinson's disease" "Thermogenesis" ...
    ##  $ GeneRatio  : chr  "14/134" "12/134" "12/134" "13/134" ...
    ##  $ BgRatio    : chr  "177/8204" "134/8204" "144/8204" "230/8204" ...
    ##  $ pvalue     : num  1.01e-06 1.70e-06 3.64e-06 9.23e-05 7.51e-04 ...
    ##  $ p.adjust   : num  0.00014 0.00014 0.0002 0.00381 0.02478 ...
    ##  $ qvalue     : num  0.000132 0.000132 0.000189 0.003595 0.023397 ...
    ##  $ geneID     : chr  "66845/68565/19896/19899/19921/19951/57808/67248/11837/68052/267019/57294/27050/20102" "12859/20463/17706/17711/17716/17717/17719/17720/17721/17722/54405/104130" "12859/20463/17706/17711/17716/17717/17719/17720/17721/17722/54405/104130" "12859/20463/12894/17706/17711/17716/17717/17719/17720/17721/17722/54405/104130" ...
    ##  $ Count      : int  14 12 12 13 9 3 7 8 6
    ## #...Citation
    ##   Guangchuang Yu, Li-Gen Wang, Yanyan Han and Qing-Yu He.
    ##   clusterProfiler: an R package for comparing biological themes among
    ##   gene clusters. OMICS: A Journal of Integrative Biology
    ##   2012, 16(5):284-287

``` r
eKG$Description
```

    ## [1] "Ribosome"                            
    ## [2] "Oxidative phosphorylation"           
    ## [3] "Parkinson's disease"                 
    ## [4] "Thermogenesis"                       
    ## [5] "Retrograde endocannabinoid signaling"
    ## [6] "Non-homologous end-joining"          
    ## [7] "Pyrimidine metabolism"               
    ## [8] "Measles"                             
    ## [9] "TGF-beta signaling pathway"

``` r
eKG$geneID
```

    ## [1] "66845/68565/19896/19899/19921/19951/57808/67248/11837/68052/267019/57294/27050/20102"
    ## [2] "12859/20463/17706/17711/17716/17717/17719/17720/17721/17722/54405/104130"            
    ## [3] "12859/20463/17706/17711/17716/17717/17719/17720/17721/17722/54405/104130"            
    ## [4] "12859/20463/12894/17706/17711/17716/17717/17719/17720/17721/17722/54405/104130"      
    ## [5] "14710/17716/17717/17719/17720/17721/17722/54405/104130"                              
    ## [6] "227525/21673/14375"                                                                  
    ## [7] "22169/66422/18102/54369/66420/20022/331487"                                          
    ## [8] "70192/19106/27103/71586/23960/246728/246727/20846"                                   
    ## [9] "12159/12166/13179/15903/16323/19877"

``` r
library(ReactomePA)
getMatrixWithSelectedIds <- function(df, type, keys,db){
require("AnnotationDbi")
require(db,character.only = TRUE)
#library(AnnotationDbi)
#library(db,character.only = TRUE) 
db <- get(db)
geneSymbols <- mapIds(db, keys=rownames(df), column=type, keytype=keys, multiVals="first")
# get the entrez ids with gene symbols i.e. remove those with NA's for gene symbols
inds <- which(!is.na(geneSymbols))

found_genes <- geneSymbols[inds]
 
# subset your data frame based on the found_genes
df2 <- df[names(found_genes), ]
rownames(df2) <- found_genes
return(df2)
}

entrezres <- getMatrixWithSelectedIds(res,type = "ENTREZID","SYMBOL","org.Mm.eg.db")
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
as.data.frame(enrichPathway(rownames(entrezres)[entrezres$padj<.1],organism="mouse"))
```

    ##                          ID
    ## R-MMU-72689     R-MMU-72689
    ## R-MMU-156827   R-MMU-156827
    ## R-MMU-72706     R-MMU-72706
    ## R-MMU-1799339 R-MMU-1799339
    ## R-MMU-975956   R-MMU-975956
    ## R-MMU-72613     R-MMU-72613
    ## R-MMU-72737     R-MMU-72737
    ## R-MMU-927802   R-MMU-927802
    ## R-MMU-975957   R-MMU-975957
    ## R-MMU-72766     R-MMU-72766
    ## R-MMU-72695     R-MMU-72695
    ## R-MMU-6791226 R-MMU-6791226
    ## R-MMU-72312     R-MMU-72312
    ## R-MMU-8868773 R-MMU-8868773
    ## R-MMU-72649     R-MMU-72649
    ## R-MMU-72702     R-MMU-72702
    ## R-MMU-72662     R-MMU-72662
    ## R-MMU-3000178 R-MMU-3000178
    ## R-MMU-2022090 R-MMU-2022090
    ## R-MMU-1474244 R-MMU-1474244
    ## R-MMU-216083   R-MMU-216083
    ## R-MMU-1474228 R-MMU-1474228
    ## R-MMU-8948216 R-MMU-8948216
    ##                                                                                                          Description
    ## R-MMU-72689                                                                 Formation of a pool of free 40S subunits
    ## R-MMU-156827                                       L13a-mediated translational silencing of Ceruloplasmin expression
    ## R-MMU-72706                                                  GTP hydrolysis and joining of the 60S ribosomal subunit
    ## R-MMU-1799339                                            SRP-dependent cotranslational protein targeting to membrane
    ## R-MMU-975956                            Nonsense Mediated Decay (NMD) independent of the Exon Junction Complex (EJC)
    ## R-MMU-72613                                                                        Eukaryotic Translation Initiation
    ## R-MMU-72737                                                                     Cap-dependent Translation Initiation
    ## R-MMU-927802                                                                           Nonsense-Mediated Decay (NMD)
    ## R-MMU-975957                               Nonsense Mediated Decay (NMD) enhanced by the Exon Junction Complex (EJC)
    ## R-MMU-72766                                                                                              Translation
    ## R-MMU-72695                                      Formation of the ternary complex, and subsequently, the 43S complex
    ## R-MMU-6791226                                          Major pathway of rRNA processing in the nucleolus and cytosol
    ## R-MMU-72312                                                                                          rRNA processing
    ## R-MMU-8868773                                                             rRNA processing in the nucleus and cytosol
    ## R-MMU-72649                                                                 Translation initiation complex formation
    ## R-MMU-72702                                                           Ribosomal scanning and start codon recognition
    ## R-MMU-72662   Activation of the mRNA upon binding of the cap-binding complex and eIFs, and subsequent binding to 43S
    ## R-MMU-3000178                                                                                      ECM proteoglycans
    ## R-MMU-2022090                                           Assembly of collagen fibrils and other multimeric structures
    ## R-MMU-1474244                                                                      Extracellular matrix organization
    ## R-MMU-216083                                                                      Integrin cell surface interactions
    ## R-MMU-1474228                                                                Degradation of the extracellular matrix
    ## R-MMU-8948216                                                                           Collagen chain trimerization
    ##               GeneRatio  BgRatio       pvalue     p.adjust       qvalue
    ## R-MMU-72689      12/121  88/8657 2.686592e-09 9.056897e-07 8.071199e-07
    ## R-MMU-156827     12/121  98/8657 9.449668e-09 9.056897e-07 8.071199e-07
    ## R-MMU-72706      12/121  99/8657 1.062758e-08 9.056897e-07 8.071199e-07
    ## R-MMU-1799339    11/121  81/8657 1.335826e-08 9.056897e-07 8.071199e-07
    ## R-MMU-975956     11/121  81/8657 1.335826e-08 9.056897e-07 8.071199e-07
    ## R-MMU-72613      12/121 106/8657 2.330642e-08 1.128697e-06 1.005856e-06
    ## R-MMU-72737      12/121 106/8657 2.330642e-08 1.128697e-06 1.005856e-06
    ## R-MMU-927802     11/121 101/8657 1.384501e-07 5.214955e-06 4.647390e-06
    ## R-MMU-975957     11/121 101/8657 1.384501e-07 5.214955e-06 4.647390e-06
    ## R-MMU-72766      14/121 212/8657 1.362494e-06 4.618856e-05 4.116167e-05
    ## R-MMU-72695       6/121  47/8657 4.426370e-05 1.111058e-03 9.901367e-04
    ## R-MMU-6791226     9/121 121/8657 4.588438e-05 1.111058e-03 9.901367e-04
    ## R-MMU-72312       9/121 121/8657 4.588438e-05 1.111058e-03 9.901367e-04
    ## R-MMU-8868773     9/121 121/8657 4.588438e-05 1.111058e-03 9.901367e-04
    ## R-MMU-72649       6/121  54/8657 9.831420e-05 2.083032e-03 1.856327e-03
    ## R-MMU-72702       6/121  54/8657 9.831420e-05 2.083032e-03 1.856327e-03
    ## R-MMU-72662       6/121  55/8657 1.091046e-04 2.175674e-03 1.938887e-03
    ## R-MMU-3000178     5/121  40/8657 2.184614e-04 4.114357e-03 3.666575e-03
    ## R-MMU-2022090     5/121  44/8657 3.448500e-04 6.152850e-03 5.483210e-03
    ## R-MMU-1474244    12/121 272/8657 4.026846e-04 6.825504e-03 6.082657e-03
    ## R-MMU-216083      6/121  76/8657 6.482689e-04 1.046491e-02 9.325973e-03
    ## R-MMU-1474228     7/121 129/8657 2.144219e-03 3.304047e-02 2.944454e-02
    ## R-MMU-8948216     4/121  41/8657 2.466764e-03 3.635796e-02 3.240098e-02
    ##                                                                                             geneID
    ## R-MMU-72689               53356/19896/19899/19921/19951/57808/11837/68052/267019/57294/27050/20102
    ## R-MMU-156827              53356/19896/19899/19921/19951/57808/11837/68052/267019/57294/27050/20102
    ## R-MMU-72706               53356/19896/19899/19921/19951/57808/11837/68052/267019/57294/27050/20102
    ## R-MMU-1799339                   19896/19899/19921/19951/57808/11837/68052/267019/57294/27050/20102
    ## R-MMU-975956                    19896/19899/19921/19951/57808/11837/68052/267019/57294/27050/20102
    ## R-MMU-72613               53356/19896/19899/19921/19951/57808/11837/68052/267019/57294/27050/20102
    ## R-MMU-72737               53356/19896/19899/19921/19951/57808/11837/68052/267019/57294/27050/20102
    ## R-MMU-927802                    19896/19899/19921/19951/57808/11837/68052/267019/57294/27050/20102
    ## R-MMU-975957                    19896/19899/19921/19951/57808/11837/68052/267019/57294/27050/20102
    ## R-MMU-72766   53356/66845/68565/19896/19899/19921/19951/57808/11837/68052/267019/57294/27050/20102
    ## R-MMU-72695                                                   53356/68052/267019/57294/27050/20102
    ## R-MMU-6791226                                14791/72544/19896/19899/19921/19951/57808/11837/20826
    ## R-MMU-72312                                  14791/72544/19896/19899/19921/19951/57808/11837/20826
    ## R-MMU-8868773                                14791/72544/19896/19899/19921/19951/57808/11837/20826
    ## R-MMU-72649                                                   53356/68052/267019/57294/27050/20102
    ## R-MMU-72702                                                   53356/68052/267019/57294/27050/20102
    ## R-MMU-72662                                                   53356/68052/267019/57294/27050/20102
    ## R-MMU-3000178                                                        12111/12842/12824/12832/13179
    ## R-MMU-2022090                                                        12842/12824/12832/12837/16948
    ## R-MMU-1474244            232345/12111/12159/12842/12824/12832/12837/13179/15891/319480/16613/16948
    ## R-MMU-216083                                                  12842/12824/12832/12837/15891/319480
    ## R-MMU-1474228                                           232345/12842/12824/12832/12837/13179/16613
    ## R-MMU-8948216                                                              12842/12824/12832/12837
    ##               Count
    ## R-MMU-72689      12
    ## R-MMU-156827     12
    ## R-MMU-72706      12
    ## R-MMU-1799339    11
    ## R-MMU-975956     11
    ## R-MMU-72613      12
    ## R-MMU-72737      12
    ## R-MMU-927802     11
    ## R-MMU-975957     11
    ## R-MMU-72766      14
    ## R-MMU-72695       6
    ## R-MMU-6791226     9
    ## R-MMU-72312       9
    ## R-MMU-8868773     9
    ## R-MMU-72649       6
    ## R-MMU-72702       6
    ## R-MMU-72662       6
    ## R-MMU-3000178     5
    ## R-MMU-2022090     5
    ## R-MMU-1474244    12
    ## R-MMU-216083      6
    ## R-MMU-1474228     7
    ## R-MMU-8948216     4

``` r
eG <- enrichGO(rownames(res[res$padj<.1,]),OrgDb ='org.Mm.eg.db',keyType = "SYMBOL",ont = "BP")
dfGO <- as.data.frame(eG)
print(dfGO[1:30,])
```

    ##                    ID
    ## GO:0051607 GO:0051607
    ## GO:0009615 GO:0009615
    ## GO:0071346 GO:0071346
    ## GO:0034341 GO:0034341
    ## GO:0006220 GO:0006220
    ## GO:0001649 GO:0001649
    ## GO:0009147 GO:0009147
    ## GO:0006303 GO:0006303
    ## GO:0000726 GO:0000726
    ## GO:0030199 GO:0030199
    ## GO:0006221 GO:0006221
    ## GO:0048333 GO:0048333
    ## GO:0035282 GO:0035282
    ## GO:0050688 GO:0050688
    ## GO:0002831 GO:0002831
    ## GO:0072527 GO:0072527
    ## GO:0043900 GO:0043900
    ## GO:0072528 GO:0072528
    ## GO:0001958 GO:0001958
    ## GO:0036075 GO:0036075
    ## GO:0044419 GO:0044419
    ## GO:0042832 GO:0042832
    ## GO:0050830 GO:0050830
    ## GO:0043901 GO:0043901
    ## GO:0042455 GO:0042455
    ## GO:0009148 GO:0009148
    ## GO:0009163 GO:0009163
    ## GO:0009220 GO:0009220
    ## GO:0046132 GO:0046132
    ## GO:0001562 GO:0001562
    ##                                                         Description
    ## GO:0051607                                defense response to virus
    ## GO:0009615                                        response to virus
    ## GO:0071346                    cellular response to interferon-gamma
    ## GO:0034341                             response to interferon-gamma
    ## GO:0006220                  pyrimidine nucleotide metabolic process
    ## GO:0001649                               osteoblast differentiation
    ## GO:0009147     pyrimidine nucleoside triphosphate metabolic process
    ## GO:0006303 double-strand break repair via nonhomologous end joining
    ## GO:0000726                               non-recombinational repair
    ## GO:0030199                             collagen fibril organization
    ## GO:0006221               pyrimidine nucleotide biosynthetic process
    ## GO:0048333                          mesodermal cell differentiation
    ## GO:0035282                                             segmentation
    ## GO:0050688                  regulation of defense response to virus
    ## GO:0002831                regulation of response to biotic stimulus
    ## GO:0072527         pyrimidine-containing compound metabolic process
    ## GO:0043900                     regulation of multi-organism process
    ## GO:0072528      pyrimidine-containing compound biosynthetic process
    ## GO:0001958                                endochondral ossification
    ## GO:0036075                                 replacement ossification
    ## GO:0044419               interspecies interaction between organisms
    ## GO:0042832                            defense response to protozoan
    ## GO:0050830              defense response to Gram-positive bacterium
    ## GO:0043901            negative regulation of multi-organism process
    ## GO:0042455                      ribonucleoside biosynthetic process
    ## GO:0009148  pyrimidine nucleoside triphosphate biosynthetic process
    ## GO:0009163                          nucleoside biosynthetic process
    ## GO:0009220           pyrimidine ribonucleotide biosynthetic process
    ## GO:0046132           pyrimidine ribonucleoside biosynthetic process
    ## GO:0001562                                    response to protozoan
    ##            GeneRatio   BgRatio       pvalue     p.adjust       qvalue
    ## GO:0051607    15/243 202/23577 3.157935e-09 9.152741e-06 8.009334e-06
    ## GO:0009615    16/243 246/23577 6.314413e-09 9.152741e-06 8.009334e-06
    ## GO:0071346     8/243  75/23577 1.055996e-06 1.020445e-03 8.929654e-04
    ## GO:0034341     8/243  94/23577 5.889357e-06 4.268312e-03 3.735092e-03
    ## GO:0006220     5/243  29/23577 1.082890e-05 6.278595e-03 5.494240e-03
    ## GO:0001649    11/243 224/23577 2.254231e-05 9.764248e-03 8.544448e-03
    ## GO:0009147     4/243  17/23577 2.357700e-05 9.764248e-03 8.544448e-03
    ## GO:0006303     5/243  38/23577 4.243819e-05 1.537854e-02 1.345737e-02
    ## GO:0000726     5/243  43/23577 7.804379e-05 1.935407e-02 1.693626e-02
    ## GO:0030199     5/243  43/23577 7.804379e-05 1.935407e-02 1.693626e-02
    ## GO:0006221     4/243  23/23577 8.355669e-05 1.935407e-02 1.693626e-02
    ## GO:0048333     4/243  23/23577 8.355669e-05 1.935407e-02 1.693626e-02
    ## GO:0035282     7/243 103/23577 9.690644e-05 1.935407e-02 1.693626e-02
    ## GO:0050688     6/243  72/23577 9.974272e-05 1.935407e-02 1.693626e-02
    ## GO:0002831     8/243 139/23577 1.001418e-04 1.935407e-02 1.693626e-02
    ## GO:0072527     5/243  48/23577 1.331328e-04 2.412199e-02 2.110855e-02
    ## GO:0043900    13/243 378/23577 1.634455e-04 2.787227e-02 2.439032e-02
    ## GO:0072528     4/243  28/23577 1.855446e-04 2.988299e-02 2.614985e-02
    ## GO:0001958     4/243  29/23577 2.134991e-04 3.094670e-02 2.708068e-02
    ## GO:0036075     4/243  29/23577 2.134991e-04 3.094670e-02 2.708068e-02
    ## GO:0044419    12/243 350/23577 2.989092e-04 3.921132e-02 3.431284e-02
    ## GO:0042832     4/243  32/23577 3.155078e-04 3.921132e-02 3.431284e-02
    ## GO:0050830     6/243  89/23577 3.207643e-04 3.921132e-02 3.431284e-02
    ## GO:0043901     8/243 165/23577 3.246194e-04 3.921132e-02 3.431284e-02
    ## GO:0042455     4/243  33/23577 3.561396e-04 4.035305e-02 3.531194e-02
    ## GO:0009148     3/243  14/23577 3.619108e-04 4.035305e-02 3.531194e-02
    ## GO:0009163     4/243  35/23577 4.484112e-04 4.487979e-02 3.927317e-02
    ## GO:0009220     3/243  15/23577 4.489527e-04 4.487979e-02 3.927317e-02
    ## GO:0046132     3/243  15/23577 4.489527e-04 4.487979e-02 3.927317e-02
    ## GO:0001562     4/243  36/23577 5.004110e-04 4.679650e-02 4.095044e-02
    ##                                                                                                      geneID
    ## GO:0051607        Ddx60/Eif2ak2/Eif2ak4/Gbp4/Ifih1/Ifit1/Oas2/Oas3/Oasl1/Oasl2/Parp9/Rtp4/Stat1/Tspan6/Zbp1
    ## GO:0009615 Ddx60/Eif2ak2/Eif2ak4/Gbp4/Ifih1/Ifit1/Oas2/Oas3/Oasl1/Oasl2/Parp9/Rps15a/Rtp4/Stat1/Tspan6/Zbp1
    ## GO:0071346                                                       Gbp10/Gbp4/Gbp6/Gbp7/Gbp8/Irf8/Parp9/Stat1
    ## GO:0034341                                                       Gbp10/Gbp4/Gbp6/Gbp7/Gbp8/Irf8/Parp9/Stat1
    ## GO:0006220                                                                      Cmpk2/Dctpp1/Nme1/Nme6/Uprt
    ## GO:0001649                                    Bmp3/Bmp4/Bmpr1a/Cat/Col1a1/Ibsp/Id3/Itga11/Runx2/Satb2/Sfrp2
    ## GO:0009147                                                                           Cmpk2/Dctpp1/Nme1/Nme6
    ## GO:0006303                                                                 Dclre1c/Ercc1/Parp9/Prpf19/Xrcc6
    ## GO:0000726                                                                 Dclre1c/Ercc1/Parp9/Prpf19/Xrcc6
    ## GO:0030199                                                                   Col1a1/Col2a1/Col5a2/Lox/Sfrp2
    ## GO:0006221                                                                             Cmpk2/Nme1/Nme6/Uprt
    ## GO:0048333                                                                          Bmp4/Bmpr1a/Inhba/Sfrp2
    ## GO:0035282                                                            Bmp4/Bmpr1a/Mafb/Nrp2/Pcsk6/Sfrp2/Ttn
    ## GO:0050688                                                            Ddx60/Eif2ak4/Gbp4/Parp9/Stat1/Tspan6
    ## GO:0002831                                                  Ddx60/Eif2ak4/Gbp4/Mif/Parp9/Stat1/Trib1/Tspan6
    ## GO:0072527                                                                      Cmpk2/Dctpp1/Nme1/Nme6/Uprt
    ## GO:0043900                    Ddx60/Eif2ak2/Eif2ak4/Gbp4/Inhba/Irf8/Mif/Oas3/Oasl1/Parp9/Stat1/Trib1/Tspan6
    ## GO:0072528                                                                             Cmpk2/Nme1/Nme6/Uprt
    ## GO:0001958                                                                         Bmp4/Col1a1/Col2a1/Runx2
    ## GO:0036075                                                                         Bmp4/Col1a1/Col2a1/Runx2
    ## GO:0044419                           Cul4a/Eif2ak2/Eif2ak4/Eif3g/Gbp6/Gbp7/Irf8/Oas3/Oasl1/Romo1/Stat1/Zbp1
    ## GO:0042832                                                                             Gbp10/Gbp6/Gbp7/Irf8
    ## GO:0050830                                                              Gbp10/Gbp6/Gbp7/Rarres2/Romo1/Rpl39
    ## GO:0043901                                               Eif2ak2/Eif2ak4/Irf8/Oas3/Oasl1/Stat1/Trib1/Tspan6
    ## GO:0042455                                                                              Aprt/Nme1/Nme6/Uprt
    ## GO:0009148                                                                                  Cmpk2/Nme1/Nme6
    ## GO:0009163                                                                              Aprt/Nme1/Nme6/Uprt
    ## GO:0009220                                                                                   Nme1/Nme6/Uprt
    ## GO:0046132                                                                                   Nme1/Nme6/Uprt
    ## GO:0001562                                                                             Gbp10/Gbp6/Gbp7/Irf8
    ##            Count
    ## GO:0051607    15
    ## GO:0009615    16
    ## GO:0071346     8
    ## GO:0034341     8
    ## GO:0006220     5
    ## GO:0001649    11
    ## GO:0009147     4
    ## GO:0006303     5
    ## GO:0000726     5
    ## GO:0030199     5
    ## GO:0006221     4
    ## GO:0048333     4
    ## GO:0035282     7
    ## GO:0050688     6
    ## GO:0002831     8
    ## GO:0072527     5
    ## GO:0043900    13
    ## GO:0072528     4
    ## GO:0001958     4
    ## GO:0036075     4
    ## GO:0044419    12
    ## GO:0042832     4
    ## GO:0050830     6
    ## GO:0043901     8
    ## GO:0042455     4
    ## GO:0009148     3
    ## GO:0009163     4
    ## GO:0009220     3
    ## GO:0046132     3
    ## GO:0001562     4

``` r
ifnGenes <- Reduce(union,strsplit(dfGO[which(grepl(pattern = "defense|interferon|immune",dfGO[,2])),"geneID"],"/"))
repairGenes <- Reduce(union,strsplit(dfGO[which(grepl(pattern = "pyrimidine|repair",dfGO[,2])),"geneID"],"/"))
bmpGenes <- Reduce(union,strsplit(dfGO[which(grepl(pattern = "ossi|osteoblast|collagen|muscle",dfGO[,2])),"geneID"],"/"))
std.heatmap(log(boneMatNorm[ifnGenes,]+1,2)-rowMeans(log(boneMatNorm[ifnGenes,]+1,2)),main="IFN response\nLog2(FC) from mean")
```

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Rowv is FALSE, while dendrogram is `both'. Omitting row
    ## dendogram.

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Colv is FALSE, while dendrogram is `column'. Omitting column
    ## dendogram.

![](BoneNotebook_files/figure-markdown_github/EnrichGO-1.png)

``` r
std.heatmap(log(boneMatNorm[repairGenes,]+1,2)-rowMeans(log(boneMatNorm[repairGenes,]+1,2)),main="DNA synth/repair\nLog2(FC) from mean")
```

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Rowv is FALSE, while dendrogram is `both'. Omitting row
    ## dendogram.

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Colv is FALSE, while dendrogram is `column'. Omitting column
    ## dendogram.

![](BoneNotebook_files/figure-markdown_github/EnrichGO-2.png)

``` r
std.heatmap(log(boneMatNorm[bmpGenes,]+1,2)-rowMeans(log(boneMatNorm[bmpGenes,]+1,2)),main="BMP Related\nLog2(FC) from mean")
```

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Rowv is FALSE, while dendrogram is `both'. Omitting row
    ## dendogram.

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Colv is FALSE, while dendrogram is `column'. Omitting column
    ## dendogram.

![](BoneNotebook_files/figure-markdown_github/EnrichGO-3.png)

``` r
std.heatmap(log(boneMatNorm[ifnGenes,]+1,2),main="IFN response\nLog2(normalized counts+1)")
```

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Rowv is FALSE, while dendrogram is `both'. Omitting row
    ## dendogram.

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Colv is FALSE, while dendrogram is `column'. Omitting column
    ## dendogram.

![](BoneNotebook_files/figure-markdown_github/EnrichGO-4.png)

``` r
std.heatmap(log(boneMatNorm[repairGenes,]+1,2),main="DNA synth/repair\nLog2(normalized counts+1)")
```

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Rowv is FALSE, while dendrogram is `both'. Omitting row
    ## dendogram.

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Colv is FALSE, while dendrogram is `column'. Omitting column
    ## dendogram.

![](BoneNotebook_files/figure-markdown_github/EnrichGO-5.png)

``` r
std.heatmap(log(boneMatNorm[bmpGenes,]+1,2),main="BMP Related\nLog2(normalized counts+1)")
```

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Rowv is FALSE, while dendrogram is `both'. Omitting row
    ## dendogram.

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Colv is FALSE, while dendrogram is `column'. Omitting column
    ## dendogram.

![](BoneNotebook_files/figure-markdown_github/EnrichGO-6.png)

#### Candice DE in the Ambrosi

``` r
std.heatmap(log(ambrosiMatNorm[ifnGenes[ifnGenes%in%rownames(ambrosiMatNorm)],]+1,2)-rowMeans(log(ambrosiMatNorm[ifnGenes[ifnGenes%in%rownames(ambrosiMatNorm)],]+1,2)),main="IFN response\nLog2(FC) from mean")
```

![](BoneNotebook_files/figure-markdown_github/checkInAmbrosi-1.png)

``` r
std.heatmap(log(ambrosiMatNorm[bmpGenes[bmpGenes%in%rownames(ambrosiMatNorm)],]+1,2)-rowMeans(log(ambrosiMatNorm[bmpGenes[bmpGenes%in%rownames(ambrosiMatNorm)],]+1,2)),main="BMP Related\nLog2(FC) from mean")
```

![](BoneNotebook_files/figure-markdown_github/checkInAmbrosi-2.png)

``` r
std.heatmap(log(ambrosiMatNorm[repairGenes[repairGenes%in%rownames(ambrosiMatNorm)],]+1,2)-rowMeans(log(ambrosiMatNorm[repairGenes[repairGenes%in%rownames(ambrosiMatNorm)],]+1,2)),main="DNA synth/repair\nLog2(FC) from mean")
```

![](BoneNotebook_files/figure-markdown_github/checkInAmbrosi-3.png)

``` r
std.heatmap(log(ambrosiMatNorm[ifnGenes[ifnGenes%in%rownames(ambrosiMatNorm)],]+1,2),main="IFN response\nLog2(normalized counts+1)")
```

![](BoneNotebook_files/figure-markdown_github/checkInAmbrosi-4.png)

``` r
std.heatmap(log(ambrosiMatNorm[repairGenes[repairGenes%in%rownames(ambrosiMatNorm)],]+1,2),main="DNA synth/repair\nLog2(normalized counts+1)")
```

![](BoneNotebook_files/figure-markdown_github/checkInAmbrosi-5.png)

``` r
std.heatmap(log(ambrosiMatNorm[bmpGenes[bmpGenes%in%rownames(ambrosiMatNorm)],]+1,2),main="BMP Related\nLog2(normalized counts+1)")
```

![](BoneNotebook_files/figure-markdown_github/checkInAmbrosi-6.png)

Ncoa1/2/3, the Steroid receptor coactivators are equally expressed across the conditions as well.

Overlap
-------

So what overlaps in the up/down for each sorted cell type and the KO vs fl data?

Note: Sca1minus is osteo, ZFP (mature) and CD24- (less mature) are adipocytes, cd24+ is multipotent progenitor

``` r
overlaps <- lapply(1:2,function(u){
  a <- sapply(ambrosiUpDown,function(x){
    c("up"=sum(boneUpDown[[u]]%in%x[[1]]), "down"=sum(boneUpDown[[u]]%in%x[[2]]))
  })
  colnames(a) <- sort(unique(sapply(strsplit(colnames(ambrosiMatNorm),"_"),function(x)x[2])))
  a
})
#Up in KO
print(overlaps[[1]])
```

    ##      CD24- CD24+ Sca1- ZFP+
    ## up       8    17    19    9
    ## down     4    18    11   17

``` r
#Down in KO
print(overlaps[[2]])
```

    ##      CD24- CD24+ Sca1- ZFP+
    ## up       1     7     2    6
    ## down     3     4     2    2

Not much. If you get loose, maybe one could say that that are up in the new data are more likely to be up in the osteocytes and down in the preadipocytes. Which genes are they?

``` r
overlaps <- lapply(1:2,function(u){
  a <- sapply(ambrosiUpDown,function(x){
    c("up"=boneUpDown[[u]][boneUpDown[[u]]%in%x[[1]]], "down"=boneUpDown[[u]][boneUpDown[[u]]%in%x[[2]]])
  })
  names(a) <- sort(unique(sapply(strsplit(colnames(ambrosiMatNorm),"_"),function(x)x[2])))
  a
})
#Up in KO
print(overlaps[[1]])
```

    ## $`CD24-`
    ##        up1        up2        up3        up4        up5        up6 
    ##    "Cadm3" "Calcoco1"     "Ccr9"   "Clec4e"      "Dcn"    "Evi2b" 
    ##        up7        up8      down1      down2      down3      down4 
    ##   "Jchain"    "Pqlc3"    "Cmpk2"     "Myl1"    "Oasl2"   "Tceanc" 
    ## 
    ## $`CD24+`
    ##       up1       up2       up3       up4       up5       up6       up7 
    ##     "A2m"  "Col2a1"   "Evi2b"  "Fkbp1b"   "Gvin1"    "Ibsp"    "Lifr" 
    ##       up8       up9      up10      up11      up12      up13      up14 
    ##   "Lpin1" "mt-Cytb" "mt-Nd4l"  "mt-Nd5"  "mt-Nd6"   "Oasl2"  "Papss2" 
    ##      up15      up16      up17     down1     down2     down3     down4 
    ##    "Peg3"   "Satb2"     "Ttn"   "Cadm3"   "Chit1"   "Cmpk2"  "Col1a1" 
    ##     down5     down6     down7     down8     down9    down10    down11 
    ##  "Col5a2"     "Dcn"  "Fam78b"  "Ifi207"   "Ifit1"   "Inhba"     "Lox" 
    ##    down12    down13    down14    down15    down16    down17    down18 
    ##    "Mlip"    "Myl1"    "Oas3"   "Postn"   "Pqlc3"      "Xk" "Zscan29" 
    ## 
    ## $`Sca1-`
    ##        up1        up2        up3        up4        up5        up6 
    ##     "Bmp3"    "Chit1"     "Cir1"    "Cmpk2"   "Col1a1"   "Col2a1" 
    ##        up7        up8        up9       up10       up11       up12 
    ##   "Gm2810"     "Ibsp"     "Mlip"  "mt-Cytb"  "mt-Nd4l"   "mt-Nd5" 
    ##       up13       up14       up15       up16       up17       up18 
    ##   "mt-Nd6"     "Myl1"     "Oas3"    "Satb2"    "Smpd3"   "Tceanc" 
    ##       up19      down1      down2      down3      down4      down5 
    ##       "Xk" "Calcoco1"     "Ccr9"    "Cerkl"    "Ddx60"    "Evi2b" 
    ##      down6      down7      down8      down9     down10     down11 
    ##   "Fkbp1b"    "Gvin1"   "Jchain"     "Pi15"     "Zbp1"   "Zfp125" 
    ## 
    ## $`ZFP+`
    ##        up1        up2        up3        up4        up5        up6 
    ##    "Cadm3"    "Cmpk2"      "Dcn"   "Fam78b"   "Fkbp1b"     "Myl1" 
    ##        up7        up8        up9      down1      down2      down3 
    ##     "Oas2"    "Pcsk6"    "Pqlc3" "Cacna2d4"    "Cadm1"     "Ccr9" 
    ##      down4      down5      down6      down7      down8      down9 
    ##   "Col2a1"   "Dixdc1"     "Ibsp"     "Lifr"  "mt-Cytb"   "mt-Nd4" 
    ##     down10     down11     down12     down13     down14     down15 
    ##  "mt-Nd4l"   "mt-Nd5"   "mt-Nd6"   "Papss2"     "Peg3"    "Satb2" 
    ##     down16     down17 
    ##    "Smpd3"       "Xk"

``` r
#Down in KO
print(overlaps[[2]])
```

    ## $`CD24-`
    ##       up    down1    down2    down3 
    ##   "Ly6d" "Rpl35a"  "Rpl39" "Rps15a" 
    ## 
    ## $`CD24+`
    ##       up1       up2       up3       up4       up5       up6       up7 
    ## "Adamts1"    "Cst3" "Eif2ak4"   "Nop10"    "Scd2"  "Tspan6"  "Zfp108" 
    ##     down1     down2     down3     down4 
    ##   "Romo1" "S100a10"  "S100a4"   "Sap30" 
    ## 
    ## $`Sca1-`
    ##       up1       up2     down1     down2 
    ##    "Emg1" "Tmem147"    "Ly6d" "Rarres2" 
    ## 
    ## $`ZFP+`
    ##       up1       up2       up3       up4       up5       up6     down1 
    ##   "Atox1"  "Mrpl33" "Rarres2"   "Romo1" "S100a10"  "S100a4"    "Cst3" 
    ##     down2 
    ##  "Gemin6"

Nothing jumps out at me...

Receptor search
---------------

Let's check the expression of a list of hormone receptors I compiled:

``` r
save.image("~/code/IngrahamLab/BoneNotebook_cache/markdown_github/everything.RData")
#load("~/code/IngrahamLab/BoneNotebook_cache/markdown_github/everything.RData")
#I looked through the literature and found what may be all the hormone receptors
receptors <- c("Esr1","Esr2","Gper1","Esrra","Esrrb","Pgr","Gnrhr","Trhr","Trhr2","Lhcgr","Ghrhr","Ghr","Ghsr","Nr4a1","Fshr","Prlhr","Pth1r","Pth2r","Prlr","Thra","Thrb","Trhr","Tshr","Crhr1","Crhr2","Mc2r",    "Mchr1","Trhr2","Mc1r","Znhit3","Kiss1r","Ar")
print(receptors)
```

    ##  [1] "Esr1"   "Esr2"   "Gper1"  "Esrra"  "Esrrb"  "Pgr"    "Gnrhr" 
    ##  [8] "Trhr"   "Trhr2"  "Lhcgr"  "Ghrhr"  "Ghr"    "Ghsr"   "Nr4a1" 
    ## [15] "Fshr"   "Prlhr"  "Pth1r"  "Pth2r"  "Prlr"   "Thra"   "Thrb"  
    ## [22] "Trhr"   "Tshr"   "Crhr1"  "Crhr2"  "Mc2r"   "Mchr1"  "Trhr2" 
    ## [29] "Mc1r"   "Znhit3" "Kiss1r" "Ar"

``` r
std.heatmap(log(ambrosiMatNorm[receptors[receptors%in%rownames(ambrosiMatNorm)],]+1,2))
```

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Rowv is FALSE, while dendrogram is `both'. Omitting row
    ## dendogram.

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Colv is FALSE, while dendrogram is `column'. Omitting column
    ## dendogram.

![](BoneNotebook_files/figure-markdown_github/hormonereceptors-1.png)

Now let's broaden the search to all the paracrine, autocrine etc receptors annotated!

``` r
descriptions <- getBM(c("mgi_symbol","mgi_description"),filters =c("mgi_symbol"),values=rownames(boneMat) ,mart = mart)
rownames(descriptions) <- descriptions$mgi_symbol
descriptions[descriptions$mgi_symbol %in% receptors,]
```

    ##        mgi_symbol                                 mgi_description
    ## Ar             Ar                               androgen receptor
    ## Crhr1       Crhr1      corticotropin releasing hormone receptor 1
    ## Crhr2       Crhr2      corticotropin releasing hormone receptor 2
    ## Esr1         Esr1                     estrogen receptor 1 (alpha)
    ## Esr2         Esr2                      estrogen receptor 2 (beta)
    ## Esrra       Esrra                estrogen related receptor, alpha
    ## Esrrb       Esrrb                 estrogen related receptor, beta
    ## Fshr         Fshr           follicle stimulating hormone receptor
    ## Ghr           Ghr                         growth hormone receptor
    ## Ghrhr       Ghrhr       growth hormone releasing hormone receptor
    ## Ghsr         Ghsr            growth hormone secretagogue receptor
    ## Gnrhr       Gnrhr         gonadotropin releasing hormone receptor
    ## Gper1       Gper1           G protein-coupled estrogen receptor 1
    ## Kiss1r     Kiss1r                                  KISS1 receptor
    ## Lhcgr       Lhcgr luteinizing hormone/choriogonadotropin receptor
    ## Mc1r         Mc1r                         melanocortin 1 receptor
    ## Mc2r         Mc2r                         melanocortin 2 receptor
    ## Mchr1       Mchr1        melanin-concentrating hormone receptor 1
    ## Nr4a1       Nr4a1 nuclear receptor subfamily 4, group A, member 1
    ## Pgr           Pgr                           progesterone receptor
    ## Prlhr       Prlhr            prolactin releasing hormone receptor
    ## Prlr         Prlr                              prolactin receptor
    ## Pth1r       Pth1r                  parathyroid hormone 1 receptor
    ## Pth2r       Pth2r                  parathyroid hormone 2 receptor
    ## Thra         Thra                  thyroid hormone receptor alpha
    ## Trhr         Trhr          thyrotropin releasing hormone receptor
    ## Trhr2       Trhr2        thyrotropin releasing hormone receptor 2
    ## Tshr         Tshr            thyroid stimulating hormone receptor
    ## Znhit3     Znhit3                         zinc finger, HIT type 3

``` r
recdesc <- descriptions[grepl("receptor",descriptions$mgi_description),]
recdesc <- recdesc[!grepl("interactor|non-receptor|interacting|ligand|associated",recdesc$mgi_description),]

#There are lots of receptors expressed in the bone stromal cell populations
heatmap.2(log(ambrosiMatNorm[recdesc$mgi_symbol[recdesc$mgi_symbol%in%rownames(ambrosiMatNorm)],]+1,2),Rowv=T,Colv = F,trace = "none",col=cols)
```

    ## Warning in heatmap.2(log(ambrosiMatNorm[recdesc$mgi_symbol[recdesc
    ## $mgi_symbol %in% : Discrepancy: Colv is FALSE, while dendrogram is `both'.
    ## Omitting column dendogram.

![](BoneNotebook_files/figure-markdown_github/moreReceptors-1.png)

``` r
hordesc <- descriptions[grepl("hormone",descriptions$mgi_description),]

resCandiceSub <- resCandice[!is.na(resCandice$padj),] 
resCandiceSub <- resCandiceSub[resCandiceSub$padj<.15,] 
resCandiceSub <-  resCandiceSub[order(resCandiceSub$log2FoldChange,decreasing = T),]

#DE receptors hard coded above
print(rownames(resCandiceSub)[rownames(resCandiceSub)%in%receptors])
```

    ## [1] "Pgr" "Ghr"

``` r
#From the list of all receptors
print(rownames(resCandiceSub)[rownames(resCandiceSub)%in%recdesc$mgi_symbol])
```

    ##  [1] "Pgr"      "Ryr3"     "Ghr"      "Csf2ra"   "Vldlr"    "Epha7"   
    ##  [7] "Bmpr1a"   "Olfr419"  "Ccr9"     "Acvr1"    "Lilr4b"   "Lilrb4a" 
    ## [13] "Rtp4"     "Rara"     "Klri2"    "Ptpre"    "Klrk1"    "Adgrg7"  
    ## [19] "Lifr"     "Ptger4"   "Tlr7"     "Il18rap"  "Tnfrsf22" "Rack1"   
    ## [25] "Rarres2"

``` r
dereceptors <-  c(rownames(resCandiceSub)[rownames(resCandiceSub)%in%recdesc$mgi_symbol],rownames(resCandiceSub)[rownames(resCandiceSub)%in%receptors])
#reverse order
dereceptors <- rownames(resCandiceSub)[rownames(resCandiceSub)%in%dereceptors]

std.heatmap(log(boneMatNorm[dereceptors[dereceptors%in%rownames(boneMatNorm)],]+1,2)-rowMeans(log(boneMatNorm[dereceptors[dereceptors%in%rownames(boneMatNorm)],]+1,2)),main="Differentially expressed receptors\n Bone marrow (FDR .15)\nlog2FC from mean")
```

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Rowv is FALSE, while dendrogram is `both'. Omitting row
    ## dendogram.

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Colv is FALSE, while dendrogram is `column'. Omitting column
    ## dendogram.

![](BoneNotebook_files/figure-markdown_github/moreReceptors-2.png)

``` r
std.heatmap(log(ambrosiMatNorm[dereceptors[dereceptors%in%rownames(ambrosiMatNorm)],]+1,2)-rowMeans(log(ambrosiMatNorm[dereceptors[dereceptors%in%rownames(ambrosiMatNorm)],]+1,2)),main="Differentially expressed receptors\n Ambrosi  (FDR .15)\nlog2FC from mean")
```

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Rowv is FALSE, while dendrogram is `both'. Omitting row
    ## dendogram.

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Colv is FALSE, while dendrogram is `column'. Omitting column
    ## dendogram.

![](BoneNotebook_files/figure-markdown_github/moreReceptors-3.png)

``` r
print(descriptions[dereceptors,])
```

    ##          mgi_symbol
    ## Pgr             Pgr
    ## Ryr3           Ryr3
    ## Ghr             Ghr
    ## Csf2ra       Csf2ra
    ## Vldlr         Vldlr
    ## Epha7         Epha7
    ## Bmpr1a       Bmpr1a
    ## Olfr419     Olfr419
    ## Ccr9           Ccr9
    ## Acvr1         Acvr1
    ## Lilr4b       Lilr4b
    ## Lilrb4a     Lilrb4a
    ## Rtp4           Rtp4
    ## Rara           Rara
    ## Klri2         Klri2
    ## Ptpre         Ptpre
    ## Klrk1         Klrk1
    ## Adgrg7       Adgrg7
    ## Lifr           Lifr
    ## Ptger4       Ptger4
    ## Tlr7           Tlr7
    ## Il18rap     Il18rap
    ## Tnfrsf22   Tnfrsf22
    ## Rack1         Rack1
    ## Rarres2     Rarres2
    ##                                                                             mgi_description
    ## Pgr                                                                   progesterone receptor
    ## Ryr3                                                                   ryanodine receptor 3
    ## Ghr                                                                 growth hormone receptor
    ## Csf2ra   colony stimulating factor 2 receptor, alpha, low-affinity (granulocyte-macrophage)
    ## Vldlr                                                 very low density lipoprotein receptor
    ## Epha7                                                                       Eph receptor A7
    ## Bmpr1a                                         bone morphogenetic protein receptor, type 1A
    ## Olfr419                                                              olfactory receptor 419
    ## Ccr9                                                       chemokine (C-C motif) receptor 9
    ## Acvr1                                                            activin A receptor, type 1
    ## Lilr4b                       leukocyte immunoglobulin-like receptor, subfamily B, member 4B
    ## Lilrb4a                      leukocyte immunoglobulin-like receptor, subfamily B, member 4A
    ## Rtp4                                                         receptor transporter protein 4
    ## Rara                                                          retinoic acid receptor, alpha
    ## Klri2                                    killer cell lectin-like receptor family I member 2
    ## Ptpre                                        protein tyrosine phosphatase, receptor type, E
    ## Klrk1                                killer cell lectin-like receptor subfamily K, member 1
    ## Adgrg7                                               adhesion G protein-coupled receptor G7
    ## Lifr                                                    leukemia inhibitory factor receptor
    ## Ptger4                                             prostaglandin E receptor 4 (subtype EP4)
    ## Tlr7                                                                   toll-like receptor 7
    ## Il18rap                                           interleukin 18 receptor accessory protein
    ## Tnfrsf22                              tumor necrosis factor receptor superfamily, member 22
    ## Rack1                                                     receptor for activated C kinase 1
    ## Rarres2                             retinoic acid receptor responder (tazarotene induced) 2

``` r
eacivector <- resCandiceSub$log2FoldChange
names(eacivector) <- rownames(resCandiceSub)


boneEACI <- eacitest(eacivector,"org.Mm.eg","SYMBOL",sets = "GO")
```

    ## Loading necessary libraries...

    ## Loaded Package org.Mm.eg.db

    ## Converting annotations to data.frames ...

    ## iteration 1 done; time  1.35 sec 
    ## iteration 2 done; time  0.14 sec 
    ## iteration 3 done; time  0.12 sec 
    ## iteration 4 done; time  0.12 sec 
    ## iteration 5 done; time  0.13 sec 
    ## iteration 6 done; time  0.13 sec 
    ## iteration 7 done; time  0.13 sec 
    ## iteration 8 done; time  0.13 sec 
    ## iteration 9 done; time  0.28 sec 
    ## iteration 10 done; time  0.13 sec

    ## Labeling output ...

    ## Loaded Package GO.db

``` r
print(boneEACI$setscores[1:30,])
```

    ##                                                                                                             Term
    ## GO:0098589                                                                                       membrane region
    ## GO:0005509                                                                                   calcium ion binding
    ## GO:0003725                                                                           double-stranded RNA binding
    ## GO:0005581                                                                                       collagen trimer
    ## GO:0030016                                                                                             myofibril
    ## GO:0015078                                                       hydrogen ion transmembrane transporter activity
    ## GO:0019843                                                                                          rRNA binding
    ## GO:0002449                                                                          lymphocyte mediated immunity
    ## GO:0048285                                                                                     organelle fission
    ## GO:0034341                                                                          response to interferon-gamma
    ## GO:0071346                                                                 cellular response to interferon-gamma
    ## GO:0007600                                                                                    sensory perception
    ## GO:0008289                                                                                         lipid binding
    ## GO:0005681                                                                                  spliceosomal complex
    ## GO:0000981                        RNA polymerase II transcription factor activity, sequence-specific DNA binding
    ## GO:0005911                                                                                    cell-cell junction
    ## GO:1990204                                                                                oxidoreductase complex
    ## GO:0110053                                                             regulation of actin filament organization
    ## GO:0010769                                          regulation of cell morphogenesis involved in differentiation
    ## GO:0000785                                                                                             chromatin
    ## GO:0070011                                                   peptidase activity, acting on L-amino acid peptides
    ## GO:0061695                                        transferase complex, transferring phosphorus-containing groups
    ## GO:0042113                                                                                     B cell activation
    ## GO:0016705 oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen
    ## GO:0001654                                                                                       eye development
    ## GO:0051301                                                                                         cell division
    ## GO:0043467                                          regulation of generation of precursor metabolites and energy
    ## GO:0042578                                                                   phosphoric ester hydrolase activity
    ## GO:0016055                                                                                 Wnt signaling pathway
    ## GO:0198738                                                                            cell-cell signaling by wnt
    ##            Ontology   set.mean    set.sd set.size         pval
    ## GO:0098589       CC  1.1974654 0.6277870        9 0.000000e+00
    ## GO:0005509       MF  1.1497516 1.1351717       12 0.000000e+00
    ## GO:0003725       MF  0.8704810 0.5197230       10 0.000000e+00
    ## GO:0005581       CC  0.8574391 0.5255333        6 2.220446e-16
    ## GO:0030016       CC  0.8458390 0.5763350        8 4.440892e-16
    ## GO:0015078       MF -0.8006743 0.2144011        9 1.449466e-14
    ## GO:0019843       MF -0.7558211 0.1986168       11 3.823360e-13
    ## GO:0002449       BP -0.5898845 1.2793780        8 1.428357e-08
    ## GO:0048285       BP  0.5833068 0.8816670       13 2.354918e-08
    ## GO:0034341       BP  0.5460414 0.4676675        9 1.730189e-07
    ## GO:0071346       BP  0.5460414 0.4676675        9 1.730189e-07
    ## GO:0007600       BP  0.5240520 0.3683686        9 5.295603e-07
    ## GO:0008289       MF  0.4912054 0.7669608       11 2.598986e-06
    ## GO:0005681       CC -0.4907059 0.3476755        7 2.374721e-06
    ## GO:0000981       MF  0.4729394 0.2945684       13 6.039757e-06
    ## GO:0005911       CC -0.4543680 1.6562069        9 1.242261e-05
    ## GO:1990204       CC -0.4054323 0.2838993        8 9.594303e-05
    ## GO:0110053       BP  0.3972824 0.2889531       10 1.451251e-04
    ## GO:0010769       BP  0.3918799 0.2843453       13 1.786456e-04
    ## GO:0000785       CC -0.3585742 0.3936405       13 5.581995e-04
    ## GO:0070011       MF  0.3450697 1.7231954        6 9.722858e-04
    ## GO:0061695       CC -0.3336838 0.2622644        9 1.316335e-03
    ## GO:0042113       BP -0.3320844 0.6587143       10 1.388396e-03
    ## GO:0016705       MF  0.3188707 0.3514336        5 2.310484e-03
    ## GO:0001654       BP  0.3042500 0.2518834       10 3.650647e-03
    ## GO:0051301       BP  0.3034110 0.1987903       11 3.745675e-03
    ## GO:0043467       BP -0.2898783 0.2357221        9 5.235120e-03
    ## GO:0042578       MF  0.2788441 0.6169008        8 7.740857e-03
    ## GO:0016055       BP  0.2777861 0.8015413        8 7.977514e-03
    ## GO:0198738       BP  0.2777861 0.8015413        8 7.977514e-03

``` r
print("done")
```

    ## [1] "done"

``` r
std.heatmap(cor.compare(boneMatNorm,ambrosiMatNorm,method="spearman"))
```

    ## [1] "Num Genes:"
    ## [1] 22492

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Rowv is FALSE, while dendrogram is `both'. Omitting row
    ## dendogram.

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Colv is FALSE, while dendrogram is `column'. Omitting column
    ## dendogram.

![](BoneNotebook_files/figure-markdown_github/moreReceptors-4.png)

``` r
IngrahamMatLog <- log(boneMatNorm+1,2)
AmbrosiMatLog <- log(ambrosiMatNorm+1,2)

informativeList <- unlist(union(unlist(ambrosiUpDown),unlist(boneUpDown)))
std.heatmap(cor.compare(AmbrosiMatLog-rowMeans(AmbrosiMatLog),IngrahamMatLog-rowMeans(IngrahamMatLog),interest.set = informativeList,method="spearman"))
```

    ## [1] "Num Genes:"
    ## [1] 3718

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Rowv is FALSE, while dendrogram is `both'. Omitting row
    ## dendogram.

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Colv is FALSE, while dendrogram is `column'. Omitting column
    ## dendogram.

![](BoneNotebook_files/figure-markdown_github/moreReceptors-5.png)

``` r
goresultsup <- getBM(attributes = c( 'mgi_symbol', 'go_id', 'name_1006'), filters = 'mgi_symbol',
                 values = rownames(resCandiceSub[resCandiceSub$log2FoldChange>0,]), mart = mart)

goresultsdown <- getBM(attributes = c( 'mgi_symbol', 'go_id', 'name_1006'), filters = 'mgi_symbol',
                 values = rownames(resCandiceSub[resCandiceSub$log2FoldChange<0,]), mart = mart)
```

``` r
sort(table(goresultsup$name_1006),decreasing = T)[15:50]
```

    ## 
    ##                                                          ATP binding 
    ##                                                                   28 
    ##                                                   biological_process 
    ##                                                                   28 
    ##                                                   cellular_component 
    ##                                                                   28 
    ##                                                        intracellular 
    ##                                                                   25 
    ##                                                        mitochondrion 
    ##                                                                   24 
    ##                                                 transferase activity 
    ##                                                                   24 
    ##                           regulation of transcription, DNA-templated 
    ##                                                                   23 
    ##                                            defense response to virus 
    ##                                                                   21 
    ##                                                          nucleoplasm 
    ##                                                                   21 
    ## positive regulation of transcription from RNA polymerase II promoter 
    ##                                                                   21 
    ##                                          oxidation-reduction process 
    ##                                                                   20 
    ##                                                            transport 
    ##                                                                   19 
    ##                                                endoplasmic reticulum 
    ##                                                                   18 
    ##                                               innate immune response 
    ##                                                                   18 
    ##                                              oxidoreductase activity 
    ##                                                                   18 
    ##                                    protein homodimerization activity 
    ##                                                                   18 
    ##                                                          RNA binding 
    ##                                                                   18 
    ##                                            identical protein binding 
    ##                                                                   17 
    ##                                                immune system process 
    ##                                                                   17 
    ## negative regulation of transcription from RNA polymerase II promoter 
    ##                                                                   17 
    ##                                                  signal transduction 
    ##                                                                   16 
    ##                                         transcription, DNA-templated 
    ##                                                                   16 
    ##                                                          DNA binding 
    ##                                                                   15 
    ##                                   multicellular organism development 
    ##                                                                   15 
    ##                                                 nucleic acid binding 
    ##                                                                   15 
    ##                                                   hydrolase activity 
    ##                                                                   14 
    ##         transcription factor activity, sequence-specific DNA binding 
    ##                                                                   14 
    ##                                                 cell differentiation 
    ##                                                                   13 
    ##                                                         cytoskeleton 
    ##                                                                   13 
    ##                                                          GTP binding 
    ##                                                                   13 
    ##                                                      immune response 
    ##                                                                   13 
    ##                  positive regulation of transcription, DNA-templated 
    ##                                                                   13 
    ##                                                 extracellular matrix 
    ##                                                                   11 
    ##                                                      Golgi apparatus 
    ##                                                                   11 
    ##                                integral component of plasma membrane 
    ##                                                                   11 
    ##                             intracellular membrane-bounded organelle 
    ##                                                                   11

``` r
sort(table(goresultsdown$name_1006),decreasing = T)[15:50]
```

    ## 
    ##                              intracellular 
    ##                                         37 
    ##                                nucleoplasm 
    ##                                         35 
    ##                          metal ion binding 
    ##                                         33 
    ##             integral component of membrane 
    ##                                         29 
    ##                         biological_process 
    ##                                         26 
    ##                            plasma membrane 
    ##                                         23 
    ##               mitochondrial inner membrane 
    ##                                         22 
    ##                                  nucleolus 
    ##                                         22 
    ##                      endoplasmic reticulum 
    ##                                         21 
    ##                             focal adhesion 
    ##                                         21 
    ## regulation of transcription, DNA-templated 
    ##                                         19 
    ##               transcription, DNA-templated 
    ##                                         18 
    ##          cytosolic large ribosomal subunit 
    ##                                         17 
    ##                                DNA binding 
    ##                                         17 
    ##                       extracellular matrix 
    ##                                         17 
    ##                       transferase activity 
    ##                                         17 
    ##                       extracellular region 
    ##                                         14 
    ##                         cellular_component 
    ##                                         13 
    ##          cytosolic small ribosomal subunit 
    ##                                         13 
    ##                         hydrolase activity 
    ##                                         13 
    ##                            rRNA processing 
    ##                                         13 
    ##                                  transport 
    ##                                         13 
    ##                  identical protein binding 
    ##                                         12 
    ##                oxidation-reduction process 
    ##                                         12 
    ##                         nucleotide binding 
    ##                                         11 
    ##          protein homodimerization activity 
    ##                                         11 
    ##                            mRNA processing 
    ##                                         10 
    ##                          apoptotic process 
    ##                                          9 
    ##             endoplasmic reticulum membrane 
    ##                                          9 
    ##                        extracellular space 
    ##                                          9 
    ##            perinuclear region of cytoplasm 
    ##                                          9 
    ##  positive regulation of cell proliferation 
    ##                                          9 
    ##                               RNA splicing 
    ##                                          9 
    ##                                ATP binding 
    ##                                          8 
    ##                               cytoskeleton 
    ##                                          8 
    ##                       nucleic acid binding 
    ##                                          8

``` r
library(GOsummaries)
gs = gosummaries(list(rownames(resCandiceSub[resCandiceSub$log2FoldChange>0,]),rownames(resCandiceSub[resCandiceSub$log2FoldChange<0,])),organism = "mmusculus")
plot(gs, fontsize = 8)
```

![](BoneNotebook_files/figure-markdown_github/goPie-1.png)

``` r
ind <- intersect(rownames(boneMatNorm),rownames(ambrosiMatNorm))

coefMat <- sapply(1:ncol(boneMatNorm),function(i)coef(lm(y~.+0,data=log(data.frame(y= boneMatNorm[ind,i],ambrosiMatNorm[ind,])+1,2))))

std.heatmap(coefMat)
```

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Rowv is FALSE, while dendrogram is `both'. Omitting row
    ## dendogram.

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Colv is FALSE, while dendrogram is `column'. Omitting column
    ## dendogram.

![](BoneNotebook_files/figure-markdown_github/whichCellType-1.png)

``` r
boneMatNormLog <- log(boneMatNorm+1,2)
ambrosiMatNormLog <- log(ambrosiMatNorm+1,2)


costF <- function(x,y, theta){
  (cor(t(y*mean(theta)),colSums(x*theta),method = "pearson")+1)/2
}

x=t(ambrosiMatNormLog[ind,])
y=t(boneMatNormLog[ind,1])
burnin <- 2000

if(!file.exists("~/Optimat.Rds")){
optiMat <- apply(boneMatNormLog[ind,],2,function(y){
paramtocost <- matrix(integer(0),nrow = 0,ncol=12)
for(i in 1:3000){
  if(i%%100 == 0)print(i)
  if(nrow(paramtocost)>burnin){
    currentParam <- sapply(1:nrow(x),function(j){
      p <- smooth.spline(paramtocost[,j] , paramtocost[,12])
      p <- predict(p,1:10000)
      py <- (p$y - min(p$y))/(max(p$y)-min(p$y))
      sample(1:10000,1,replace = T,prob = (py)/(sum(py)))
    })
  }else if(i==1)currentParam <- rep(1,nrow(x))
  else if(i==2)currentParam <- rep(10000,nrow(x))
  else{currentParam <- sample(1:10000,nrow(x),replace = T)}
  cost <- costF(t(ambrosiMatNormLog[ind,]),t(boneMatNormLog[ind,1]),currentParam)
  paramtocost <- rbind(paramtocost,c(currentParam,cost))
}
paramtocost <- paramtocost[order(paramtocost[,ncol(paramtocost)]),]
print(paramtocost[nrow(paramtocost),ncol(paramtocost)]*2-1)
print(cor(t(x),y))
paramtocost[nrow(paramtocost),]
})

saveRDS(optiMat,file = "~/Optimat.Rds")
}
optiMat <- readRDS("~/Optimat.Rds")
std.heatmap(cor.compare(ambrosiMatNormLog,boneMatNormLog,min = 1),method="pearson")
```

    ## [1] "Num Genes:"
    ## [1] 18791

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Rowv is FALSE, while dendrogram is `both'. Omitting row
    ## dendogram.

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Colv is FALSE, while dendrogram is `column'. Omitting column
    ## dendogram.

    ## Warning in plot.window(...): "method" is not a graphical parameter

    ## Warning in plot.xy(xy, type, ...): "method" is not a graphical parameter

    ## Warning in title(...): "method" is not a graphical parameter

![](BoneNotebook_files/figure-markdown_github/whichCellType-2.png)

``` r
normOptiMat <- sapply(1:ncol(optiMat),function(i) optiMat[-12,i]/sum(optiMat[-12,i]))
colnames(normOptiMat) <- colnames(optiMat)
rownames(normOptiMat) <- colnames(ambrosiMatNorm)
std.heatmap(normOptiMat)
```

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Rowv is FALSE, while dendrogram is `both'. Omitting row
    ## dendogram.

    ## Warning in heatmap.2(M, Rowv = F, Colv = F, trace = "none", col = cols, :
    ## Discrepancy: Colv is FALSE, while dendrogram is `column'. Omitting column
    ## dendogram.

![](BoneNotebook_files/figure-markdown_github/whichCellType-3.png)

``` r
grad.descent <- function(x, maxit){
    theta <- matrix(c(0, 0), nrow=1) # Initialize the parameters
 
    alpha = .05 # set learning rate
    for (i in 1:maxit) {
      theta <- theta - alpha  * grad(x, y, theta)   
    }
 return(theta)
}

grad <- function(x, y, theta) {
  gradient <- (1/m)* (t(x) %*% ((x %*% t(theta)) - y))
  return(t(gradient))
}
```
