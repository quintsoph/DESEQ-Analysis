#nohup R CMD BATCH ./myprog.R &
library("RColorBrewer")
 library("DESeq2")
 library("geneplotter")
 library("EDASeq")
 library("genefilter")
 library("pheatmap")
 library ("locfit")
 source("https://bioconductor.org/biocLite.R")
biocLite()
#Load file
	#setwd("/DataDrives/dd5/dbGaP/ncbi/dbGaP-6000/files/phe000020.v1.GTEx_RNAseq.expression-data-matrixfmt.c1")
	#batchnormal<- read.table("GTEx_Data_20160115_v7_RNAseq_RNASeQCv1.1.8_gene_reads.gct", sep="\t", fill=TRUE, skip=2, header=TRUE)
	#normal<- data.frame(batchnormal[,-1])
#remove genes
	#normal_matrix <- as.matrix(normal[,-1])
	#rownames(normal_matrix) <- normal[,1]

#setup covairate tissue type							  
	#setwd("/home/quints1/gtexstuff/deseqstuff")
	#attrib<-read.table("patienttissue.txt", fill=TRUE, sep="\t", row.names=1)
	#colnames(attrib)= c("type")
	#library("stringr")
	#rownames(attrib)=str_replace_all(rownames(attrib), "-", ".")
	#matchcut_new <- normal_matrix[,colnames(normal_matrix) %in% rownames(attrib)]
	#currentpat<- as.matrix(colnames(matchcut_new))
	#attrib2<- as.matrix(cbind(rownames(attrib), attrib))
	#name.use<- rownames(attrib)[rownames(attrib) %in% currentpat]
	#subsetattrib<- attrib2[name.use,]
	#colnames(subsetattrib)= c("patient","type")
	#thing<- subsetattrib[order(colnames(matchcut_new)),]
	#rownames(thing)=NULL

#save the files
saveRDS(matchcut_new, "patientgene.rds")
saveRDS(thing, "patientsattrib.rds")

#create the R class
matchcut_new<- readRDS("patientgene.rds")
thing<-readRDS("patientsattrib.rds")

dds <- DESeqDataSetFromMatrix(countData = matchcut_new,
							colData = as.data.frame(thing),
                              design = ~ type)

# step 1
  FDR <- 0.1
  floorPDEG <- 0.05
  filter <- apply(counts(dds), 1, function(x) { all(x > 0)})
  dds <- dds[filter,]
  #size factors
ddsnormal <- estimateSizeFactors(dds)
    sizeFactors(ddsnormal) <- sizeFactors(ddsnormal) / mean(sizeFactors(ddsnormal))
    sf <- sizeFactors(ddsnormal) 
    write.table(sf, file="sizefactors.1.txt", quote=FALSE, row.names=TRUE, sep="\t")

  library("locfit")
# step 2
  if (file.exists("ddsGenes.LRT.RData")) {
    load("ddsGenes.LRT.RData")
  } else {
    ddsnormal <- estimateDispersions(ddsnormal)
    ddsnormal2 <- nbinomLRT(ddsnormal, reduced = ~ 1, maxit=500)
    if (any(mcols(ddsnormal2)$betaConv)) {
      ddsnormal2 <- ddsnormal2[which(mcols(ddsnormal2)$betaConv),]
    }
    save(ddsnormal2, file="ddsGenes.LRT.RData")
  }
 
  if (file.exists("res.LRT.RData")) {
    load("res.LRT.RData")
  } else {
    res <- results(ddsnormal2)
    save(res, file="res.LRT.RData")
  }
 
  FDR = 1.0e-50
  floorPDEG = 1.0e-50
  pval <- res$pvalue
  pval[is.na(pval)] <- 1
  qval <- p.adjust(pval, method = "BH")
  if (sum(qval < FDR) > (floorPDEG * nrow(ddsnormal2))) {
      is.DEG <- as.logical(qval < FDR)
  } else {
    is.DEG <- as.logical(rank(pval, ties.method = "min") <= nrow(ddsnormal2) * floorPDEG)
  }

 
  # step 3
  nonDEGgenes <- rownames(counts(ddsnormal2, normalized=FALSE)[!is.DEG,])
  write.table(nonDEGgenes, file="nonDEGgenes.txt", quote=FALSE, row.names=TRUE, sep="\t")
ctrlgenesidx <-  match(nonDEGgenes, rownames(matchcut_new))
  ddsGenes <- DESeqDataSetFromMatrix(
			countData = matchcut_new,
				colData = as.data.frame(thing),
                 design = ~ type)
  
  ddsGenes <- estimateSizeFactors(ddsGenes, controlGenes = ctrlgenesidx)
 write.table(sizeFactors(ddsGenes), file="sizefactors.2.txt", quote=FALSE, row.names=TRUE, sep="\t")
  #save(ddsGenes, file="ddsGenes.DEGES.RData")
 
  normalcnts <- counts(ddsGenes, normalized=TRUE)
  idx.nz <- apply(normalcnts, 1, function(x) { all(x > 0)})
  
  sf <- read.table("sizefactors.2.txt", row.names=1, sep="\t", header=TRUE)
  sf <- c(t(sf))
  sizeFactors(ddsGenes) <- sf
  
  colnormal = colData(ddsGenes)

###Input TEtranscripts
totaltissue<- read.table("notissueL1HS.txt", header=TRUE, row.names=1, sep="\t")
ttissue<- t(totaltissue)
library("stringr")
colnames(ttissue)=str_replace_all(colnames(ttissue), "-", ".")
t_TEcntmatrix = t(ttissue)
 t_cntmatrix = t(matchcut_new)
m <- merge(t_cntmatrix, t_TEcntmatrix, by.x="row.names", by.y="row.names")
new_patients = m[,1]
new_genenames = colnames(m[,-1])
new_cntmatrix = as.matrix(t(m[,-1]))
colnames(new_cntmatrix) = new_patients
#subset out the correct patients
library("dplyr")
subsetedmatrix<- data.frame(matchcut_new) %>% select(one_of(dput(as.character(new_patients))))
  subsetedmatrix<- as.matrix(subsetedmatrix)
 rownames(subsetedmatrix) = rownames(matchcut_new)
 new_thing<- data.frame(thing) %>% filter(patient %in% new_patients)

 ddsGenes <- DESeqDataSetFromMatrix(
			countData = subsetedmatrix,
				colData = as.data.frame(new_thing),
                 design = ~ type)
  
  ddsGenes <- estimateSizeFactors(ddsGenes, controlGenes = ctrlgenesidx)
 write.table(sizeFactors(ddsGenes), file="sizefactors.subset.txt", quote=FALSE, row.names=TRUE, sep="\t")
 colnormal = colData(ddsGenes)
 
 sf <- read.table("sizefactors.subset.txt", row.names=1, sep="\t", header=TRUE)
  sf <- c(t(sf))
  sizeFactors(ddsGenes) <- sf
  
#write.table(new_cntmatrix, file="cntmatrix.geneTE.txt", quote=FALSE, row.names=TRUE, sep="\t")
 
  stopifnot( as.character(colnormal$patient) == new_patients ) # need patient ids in same order to use previously estimated sizeFactors
  colnames(new_cntmatrix) = NULL
  ddsnew <- DESeqDataSetFromMatrix(
    countData = new_cntmatrix,
    colData = colnormal,
    design = ~ type )
  ddsnew$type <- droplevels(ddsnew$type)
  # set sizefactors 
  sizeFactors(ddsnew) <- sf
  save(ddsnew, file="ddsNew.DEGES.RData")
  ddsnewcnts <- counts (ddsnew, normalized=TRUE) 
  write.table(ddsnewcnts, file="normalizedcnt.txt", quote=FALSE, row.names=TRUE, sep="\t")
 
 
  L1HSnormalized <- ddsnewcnts["L1HS",]
  write.table(L1HSnormalized, file="L1HS.normalized.txt", quote=FALSE, row.names=TRUE, sep="\t")
  log2colsums <- log2(colSums(ddsnewcnts))
  write.table(log2colsums, file="log2colsums.txt", quote=FALSE, row.names=TRUE, sep="\t")
  
#log transformation (vst)
 if (file.exists("VSTcnts.DEGES.RData")) {
  load(file = "VSTcnts.DEGES.RData")
 } else {
  vsd <- vst(ddsnew, blind = FALSE)
  VSTcnts <- assay(vsd)
  write.table(VSTcnts["L1HS",], "L1HS.VST.cnts.txt", quote=FALSE, row.names = TRUE)
  save(VSTcnts, file = "VSTcnts.DEGES.RData")
  write.table(VSTcnts, file="VSTcnt.txt", quote=FALSE, row.names=TRUE, sep="\t")
  
  genenames = rownames(ddsnew)
  genenames <- gsub("ALR/Alpha", "ALR_Alpha", genenames)
  genenames <- gsub("BSR/Beta", "BSR_Beta", genenames)
  rownames(ddsnew) <- genenames
  coldataNew = colData(ddsnew)
  for (i in 1:nrow(VSTcnts)) {
    fname = paste0("VSTcnts.normal/",genenames[i], ".txt")
    cnts <- matrix(VSTcnts[i,], ncol=1)
    rownames(cnts) <- coldataNew$patient
    colnames(cnts) <- c("patient\tgene")
    write.table(cnts, file=fname, quote=FALSE, sep="\t", row.names=TRUE)
  }
  L1HSnormalized = read.table("L1HS.normalized.txt", row.names=1, sep="\t", check.names = FALSE)
  L1HStable <- cbind.data.frame(colData(ddsnew), L1HSnormalized )
  L1HStable <- cbind.data.frame(L1HStable, VSTcnts["L1HS",] )
  colnames(L1HStable)[7:8] = c("normalizedcnt", "VSTcnt")
  write.table(L1HStable, file="L1HS.VST.txt", quote=FALSE, row.names=TRUE, sep="\t")
 }
 
  
###Dr. Han script						  
  genenames <- rownames(read.table(paste("RawCnts",as.character(patient_types[1,2]),sep="/"), header=TRUE, row.names=1))
  #samplenames <- paste(rep(patient_types[,1], each=2), rep(c("cancer", "normal"), nrow(patient_types)), sep="_")
  #coldata <- cbind.data.frame(rep(patient_types[,1], each=2), rep(patient_types[,3], each=2), rep(c("cancer", "normal"), nrow(patient_types)), as.vector(t(patient_types[,4:5])))
  coldata <- cbind.data.frame(rep(patient_types[,1], each=1), rep(patient_types[,3], each=1), rep(c("normal"), nrow(patient_types)), as.vector(t(patient_types[,4])))
  #rownames(coldata) <- samplenames 
  rownames(coldata) <- patient_types[,1] 
  colnames(coldata) <- c("patient", "type", "condition", "batch")
  coldata$batch <- factor(coldata$batch)
  
  #need to make batch nested within type 
  coldata <- cbind.data.frame(coldata, rep(0, nrow(coldata)))
  colnames(coldata)[5] <- "nested.batch"
  for (i in 1:length(cancertypes)) {
    coldataLUAD <- coldata[coldata$type==cancertypes[i],]
    coldataLUAD$batch <- droplevels(coldataLUAD$batch)
    coldata[coldata$type==cancertypes[i],"nested.batch"]=as.numeric(coldataLUAD$batch)
  }
  coldata$nested.batch = factor(coldata$nested.batch)
  dim(coldata)
 
  cntmatrix = c() 
  if (file.exists("cntmatrix.gene.txt")) {
    cntmatrix = read.table(file="cntmatrix.gene.txt", header=TRUE, row.names=1, sep="\t", check.names = FALSE)
    cntGenes <- cntmatrix
  } else { 
    cntmatrix = matrix(rep(0, nrow(patient_types)*length(genenames)),  nrow=length(genenames))
    for (i in 1:nrow(patient_types)) {
      cnts <- read.table(paste("RawCnts",as.character(patient_types[i,2]),sep="/"), header=TRUE, row.names=1)
      #cntmatrix[,(i*2-1):(i*2)]=as.matrix(cnts)
      cntmatrix[,i]=as.matrix(cnts)
    }
    dim(cntmatrix)
    rownames(cntmatrix) <- genenames
    cntGenes <- cntmatrix
    colnames(cntmatrix) <- coldata$patient
    write.table(cntmatrix, file="cntmatrix.gene.txt", quote=FALSE, row.names=TRUE, sep="\t")
  } 
  
  ddsnormal <- DESeqDataSetFromMatrix(
    countData = cntGenes,
    colData = coldata,
    design = ~ nested.batch + type )
  #ddsnormal <- ddsnormal[, ddsnormal$type!="CESC"]
  #ddsnormal <- ddsnormal[, ddsnormal$type!="PAAD"]
  #ddsnormal <- ddsnormal[, ddsnormal$type!="PCPG"]
  #ddsnormal <- ddsnormal[, ddsnormal$type!="SARC"]
  #ddsnormal <- ddsnormal[, ddsnormal$type!="THYM"]
  
  ddsnormal$type <- droplevels(ddsnormal$type)
  #keep <- rowSums(counts(ddsnormal)) >= 10
  #dds <- dds[keep,]
 
  # step 1
  FDR <- 0.1
  floorPDEG <- 0.05
  filter <- apply(counts(ddsnormal), 1, function(x) { all(x > 0)})
  ddsnormal <- ddsnormal[filter,]

  
  if (file.exists("sizefactors.1.txt")) {
    sf <- read.table("sizefactors.1.txt", row.names=1, sep="\t", header=TRUE) 
    sizeFactors(ddsnormal) <- c(t(sf))
  } else {
    ddsnormal <- estimateSizeFactors(ddsnormal)
    sizeFactors(ddsnormal) <- sizeFactors(ddsnormal) / mean(sizeFactors(ddsnormal))
    sf <- sizeFactors(ddsnormal) 
    write.table(sf, file="sizefactors.1.txt", quote=FALSE, row.names=TRUE, sep="\t")
  }
 
# step 2
  if (file.exists("ddsGenes.LRT.RData")) {
    load("ddsGenes.LRT.RData")
  } else {
    ddsnormal <- estimateDispersions(ddsnormal)
    ddsnormal <- nbinomLRT(ddsnormal, reduced = ~ nested.batch, maxit=500)
    if (any(mcols(ddsnormal)$betaConv)) {
      ddsnormal <- ddsnormal[which(mcols(ddsnormal)$betaConv),]
    }
    save(ddsnormal, file="ddsGenes.LRT.RData")
  }
 
  if (file.exists("res.LRT.RData")) {
    load("res.LRT.RData")
  } else {
    res <- results(ddsnormal)
    save(res, file="res.LRT.RData")
  }
 
  FDR = 1.0e-50
  floorPDEG = 1.0e-50
  pval <- res$pvalue
  pval[is.na(pval)] <- 1
  qval <- p.adjust(pval, method = "BH")
  if (sum(qval < FDR) > (floorPDEG * nrow(ddsnormal))) {
      is.DEG <- as.logical(qval < FDR)
  } else {
    is.DEG <- as.logical(rank(pval, ties.method = "min") <= nrow(ddsnormal) * floorPDEG)
  }
 
  # step 3
  nonDEGgenes <- rownames(counts(ddsnormal2, normalized=FALSE)[!is.DEG,])
  write.table(nonDEGgenes, file="nonDEGgenes.txt", quote=FALSE, row.names=TRUE, sep="\t")
  ctrlgenesidx <-  match(nonDEGgenes, rownames(cntmatrix))
  ddsGenes <- DESeqDataSetFromMatrix(
    countData = cntmatrix,
    colData = coldata,
    design = ~ nested.batch + type )
  
  ddsGenes <- estimateSizeFactors(ddsGenes, controlGenes = ctrlgenesidx)
 
  #norm.factors <- sizeFactors(ddsGenes)/colSums(cntGenes)
  #norm.factors <- norm.factors/mean(norm.factors)
  #sizeFactors(ddsGenes) <- norm.factors
  write.table(sizeFactors(ddsGenes), file="sizefactors.2.txt", quote=FALSE, row.names=TRUE, sep="\t")
  save(ddsGenes, file="ddsGenes.DEGES.RData")
 
  normalcnts <- counts(ddsGenes, normalized=TRUE)
  idx.nz <- apply(normalcnts, 1, function(x) { all(x > 0)})
  sum(idx.nz)
  pdf("multidensity.pdf")
  multidensity( normalcnts[idx.nz ,],
                xlab="mean counts", xlim=c(0, 1000))
  dev.off()
  pdf("multiecdf.pdf")
  multiecdf( normalcnts[idx.nz ,],
             xlab="mean counts", xlim=c(0, 1000))
  dev.off()
  
 }
 
 
 
 if (file.exists("ddsNew.DEGES.RData")) {
  load(file = "ddsNew.DEGES.RData")
 } else {
 
  cntmatrix = c() 
  if (file.exists("cntmatrix.gene.txt")) {
    cntmatrix = read.table(file="cntmatrix.gene.txt", header=TRUE, row.names=1, sep="\t", check.names = FALSE)
    cntGenes <- cntmatrix
  } else { 
    cntmatrix = matrix(rep(0, nrow(patient_types)*length(genenames)),  nrow=length(genenames))
    for (i in 1:nrow(patient_types)) {
      cnts <- read.table(paste("RawCnts",as.character(patient_types[i,2]),sep="/"), header=TRUE, row.names=1)
      #cntmatrix[,(i*2-1):(i*2)]=as.matrix(cnts)
      cntmatrix[,i]=as.matrix(cnts)
    }
    dim(cntmatrix)
    rownames(cntmatrix) <- genenames
    cntGenes <- cntmatrix
    colnames(cntmatrix) <- coldata$patient
    write.table(cntmatrix, file="cntmatrix.gene.txt", quote=FALSE, row.names=TRUE, sep="\t")
  } 
  
  ddsGenes <- DESeqDataSetFromMatrix(
    countData = cntGenes,
    colData = coldata,
    design = ~ nested.batch + type )
  
  ddsGenes$type <- droplevels(ddsGenes$type)
 
  sf <- read.table("sizefactors.2.txt", row.names=1, sep="\t", header=TRUE)
  sf <- c(t(sf))
  sizeFactors(ddsGenes) <- sf
  
  colnormal = colData(ddsGenes)
  
 # construct counts from discounted elements 
 # create the discounted elements
  file.names <- dir("RawCnts", pattern =".discount.ele.cntTable")
  patientIDs <- substr(file.names, 9, 12)
  TEnames <- rownames(read.table(paste("RawCnts",as.character(file.names[1]),sep="/"), header=TRUE, row.names=1))
 
  TEcntmatrix = matrix(rep(0, length(file.names)*length(TEnames)),  nrow=length(TEnames))
  for (i in 1:length(file.names)) {
    cnts <- read.table(paste("RawCnts",as.character(file.names[i]),sep="/"), header=TRUE, row.names=1)
    TEcntmatrix[,i]=as.matrix(cnts)
  }
  dim(TEcntmatrix)
  rownames(TEcntmatrix) <- TEnames
  colnames(TEcntmatrix) <- patientIDs
  
 # merge with existing gene cntmatrix
  t_TEcntmatrix = t(TEcntmatrix)
  t_cntmatrix = t(cntmatrix)
  
  m <- merge(t_cntmatrix, t_TEcntmatrix, by.x="row.names", by.y="row.names")
  new_patients = m[,1]
  new_genenames = colnames(m[,-1])
  new_cntmatrix = t(m[,-1])
  colnames(new_cntmatrix) = new_patients
  write.table(new_cntmatrix, file="cntmatrix.geneTE.txt", quote=FALSE, row.names=TRUE, sep="\t")
 
  stopifnot( colnormal$patient == new_patients ) # need patient ids in same order to use previously estimated sizeFactors
  colnames(new_cntmatrix) = NULL
  ddsnew <- DESeqDataSetFromMatrix(
    countData = new_cntmatrix,
    colData = colnormal,
    design = ~ nested.batch + type )
  ddsnew$type <- droplevels(ddsnew$type)
  # set sizefactors 
  sizeFactors(ddsnew) <- sf
  save(ddsnew, file="ddsNew.DEGES.RData")
  ddsnewcnts <- counts (ddsnew, normalized=TRUE)
  write.table(ddsnewcnts, file="normalizedcnt.txt", quote=FALSE, row.names=TRUE, sep="\t")
 
  L1HSnormalized <- ddsnewcnts["L1HS:L1:LINE",]
  write.table(L1HSnormalized, file="L1HS.normalized.txt", quote=FALSE, row.names=TRUE, sep="\t")
  log2colsums <- log2(colSums(ddsnewcnts))
  write.table(log2colsums, file="log2colsums.txt", quote=FALSE, row.names=TRUE, sep="\t")
  
 }
 
 
 
 
 #log transformation (vst)
 if (file.exists("VSTcnts.DEGES.RData")) {
  load(file = "VSTcnts.DEGES.RData")
 } else {
  vsd <- vst(ddsnew, blind = FALSE)
  VSTcnts <- assay(vsd)
  write.table(VSTcnts["L1HS:L1:LINE",], "L1HS.VST.cnts.txt", quote=FALSE, row.names = TRUE)
  save(VSTcnts, file = "VSTcnts.DEGES.RData")
  write.table(VSTcnts, file="VSTcnt.txt", quote=FALSE, row.names=TRUE, sep="\t")
  
  genenames = rownames(ddsnew)
  genenames <- gsub("ALR/Alpha", "ALR_Alpha", genenames)
  genenames <- gsub("BSR/Beta", "BSR_Beta", genenames)
  rownames(ddsnew) <- genenames
  coldataNew = colData(ddsnew)
  for (i in 1:nrow(VSTcnts)) {
    fname = paste0("VSTcnts.normal/",genenames[i], ".txt")
    cnts <- matrix(VSTcnts[i,], ncol=1)
    rownames(cnts) <- coldataNew$patient
    colnames(cnts) <- c("patient\tgene")
    write.table(cnts, file=fname, quote=FALSE, sep="\t", row.names=TRUE)
  }
  L1HSnormalized = read.table("L1HS.normalized.txt", row.names=1, sep="\t", check.names = FALSE)
  L1HStable <- cbind.data.frame(colData(ddsnew), L1HSnormalized )
  L1HStable <- cbind.data.frame(L1HStable, VSTcnts["L1HS:L1:LINE",] )
  colnames(L1HStable)[7:8] = c("normalizedcnt", "VSTcnt")
  write.table(L1HStable, file="L1HS.VST.txt", quote=FALSE, row.names=TRUE, sep="\t")
 }
 