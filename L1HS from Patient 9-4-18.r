###Gather data for L1HS of each pat and put into one file
#Test with one file
setwd("/DataDrives/dd5/gtex/Adrenal_Gland")
pat<- read.table("GTEX-111CU_Adrenal_Gland/GTEX-111CU_Adrenal_Gland.discount.ele.cntTable")
rownames(pat)=pat[,1]
L1HSpat<-pat[c("L1HS:L1:LINE"),]
final<-cbind(L1HSpat, "GTEX-111CU_Adrenal_Gland")
colnames(final)=c("transposon", "expression", "patient")

#create for all tissue pats
#Adrenal_Gland
data.frame(patientlist<- read.table("/home/quints1/gtexstuff/patientnames/Adrenalpats.txt"))
setwd("/DataDrives/dd5/gtex/Adrenal_Gland")
allpats=NULL
for (i in patientlist[,1])
{
	final=NULL
	patfile=NULL
	patfile= paste(i,"/",i,".discount.ele.cntTable", sep="")
	pat<- read.table(patfile)
	rownames(pat)=pat[,1]
	L1HSpat<-pat[c("L1HS:L1:LINE"),]
	final<-cbind(L1HSpat, paste(i))
	colnames(final)=c("transposon", "expression", "patient")
	allpats<- rbind(allpats, final)
}

setwd("/home/quints1/gtexstuff/L1HSpatients")
write.table(allpats, "AdrenalGlandL1HSresults.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

#Breast
data.frame(patientlist<- read.table("/home/quints1/gtexstuff/patientnames/Breastpats.txt"))
setwd("/DataDrives/dd5/gtex/Breast")
allpats=NULL
for (i in patientlist[,1])
{
	final=NULL
	patfile=NULL
	patfile= paste(i,"/",i,".discount.ele.cntTable", sep="")
	pat<- read.table(patfile)
	rownames(pat)=pat[,1]
	L1HSpat<-pat[c("L1HS:L1:LINE"),]
	final<-cbind(L1HSpat, paste(i))
	colnames(final)=c("transposon", "expression", "patient")
	allpats<- rbind(allpats, final)
}

setwd("/home/quints1/gtexstuff/L1HSpatients")
write.table(allpats, "BreastL1HSresults.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

#Colon
data.frame(patientlist<- read.table("/home/quints1/gtexstuff/patientnames/Colonpats.txt"))
setwd("/DataDrives/dd5/gtex/Colon")
allpats=NULL
for (i in patientlist[,1])
{
	final=NULL
	patfile=NULL
	patfile= paste(i,"/",i,".discount.ele.cntTable", sep="")
	pat<- read.table(patfile)
	rownames(pat)=pat[,1]
	L1HSpat<-pat[c("L1HS:L1:LINE"),]
	final<-cbind(L1HSpat, paste(i))
	colnames(final)=c("transposon", "expression", "patient")
	allpats<- rbind(allpats, final)
}

setwd("/home/quints1/gtexstuff/L1HSpatients")
write.table(allpats, "ColonL1HSresults.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

#Kidney
data.frame(patientlist<- read.table("/home/quints1/gtexstuff/patientnames/Kidneypats.txt"))
setwd("/DataDrives/dd5/gtex/Kidney")
allpats=NULL
for (i in patientlist[,1])
{
	final=NULL
	patfile=NULL
	patfile= paste(i,"/",i,".discount.ele.cntTable", sep="")
	pat<- read.table(patfile)
	rownames(pat)=pat[,1]
	L1HSpat<-pat[c("L1HS:L1:LINE"),]
	final<-cbind(L1HSpat, paste(i))
	colnames(final)=c("transposon", "expression", "patient")
	allpats<- rbind(allpats, final)
}

setwd("/home/quints1/gtexstuff/L1HSpatients")
write.table(allpats, "KidneyL1HSresults.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

#Pancreas
data.frame(patientlist<- read.table("/home/quints1/gtexstuff/patientnames/Pancreaspats.txt"))
setwd("/DataDrives/dd5/gtex/Pancreas")
allpats=NULL
for (i in patientlist[,1])
{
	final=NULL
	patfile=NULL
	patfile= paste(i,"/",i,".discount.ele.cntTable", sep="")
	pat<- read.table(patfile)
	rownames(pat)=pat[,1]
	L1HSpat<-pat[c("L1HS:L1:LINE"),]
	final<-cbind(L1HSpat, paste(i))
	colnames(final)=c("transposon", "expression", "patient")
	allpats<- rbind(allpats, final)
}

setwd("/home/quints1/gtexstuff/L1HSpatients")
write.table(allpats, "PancreasL1HSresults.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

#Stomach
data.frame(patientlist<- read.table("/home/quints1/gtexstuff/patientnames/Stomachpats.txt"))
setwd("/DataDrives/dd5/gtex/Stomach")
allpats=NULL
for (i in patientlist[,1])
{
	final=NULL
	patfile=NULL
	patfile= paste(i,"/",i,".discount.ele.cntTable", sep="")
	pat<- read.table(patfile)
	rownames(pat)=pat[,1]
	L1HSpat<-pat[c("L1HS:L1:LINE"),]
	final<-cbind(L1HSpat, paste(i))
	colnames(final)=c("transposon", "expression", "patient")
	allpats<- rbind(allpats, final)
}

setwd("/home/quints1/gtexstuff/L1HSpatients")
write.table(allpats, "StomachL1HSresults.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

#Bladder
data.frame(patientlist<- read.table("/home/quints1/gtexstuff/patientnames/Bladderpats.txt"))
setwd("/DataDrives/dd5/gtex/Bladder")
allpats=NULL
for (i in patientlist[,1])
{
	final=NULL
	patfile=NULL
	patfile= paste(i,"/",i,".discount.ele.cntTable", sep="")
	pat<- read.table(patfile)
	rownames(pat)=pat[,1]
	L1HSpat<-pat[c("L1HS:L1:LINE"),]
	final<-cbind(L1HSpat, paste(i))
	colnames(final)=c("transposon", "expression", "patient")
	allpats<- rbind(allpats, final)
}

setwd("/home/quints1/gtexstuff/L1HSpatients")
write.table(allpats, "BladderL1HSresults.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

#Cervix
data.frame(patientlist<- read.table("/home/quints1/gtexstuff/patientnames/Cervixpats.txt"))
setwd("/DataDrives/dd5/gtex/Cervix_Uteri")
allpats=NULL
for (i in patientlist[,1])
{
	final=NULL
	patfile=NULL
	patfile= paste(i,"/",i,".discount.ele.cntTable", sep="")
	pat<- read.table(patfile)
	rownames(pat)=pat[,1]
	L1HSpat<-pat[c("L1HS:L1:LINE"),]
	final<-cbind(L1HSpat, paste(i))
	colnames(final)=c("transposon", "expression", "patient")
	allpats<- rbind(allpats, final)
}

setwd("/home/quints1/gtexstuff/L1HSpatients")
write.table(allpats, "CervixL1HSresults.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

#Esophagus
data.frame(patientlist<- read.table("/home/quints1/gtexstuff/patientnames/Esophaguspats.txt"))
setwd("/DataDrives/dd5/gtex/Esophagus")
allpats=NULL
for (i in patientlist[,1])
{
	final=NULL
	patfile=NULL
	patfile= paste(i,"/",i,".discount.ele.cntTable", sep="")
	pat<- read.table(patfile)
	rownames(pat)=pat[,1]
	L1HSpat<-pat[c("L1HS:L1:LINE"),]
	final<-cbind(L1HSpat, paste(i))
	colnames(final)=c("transposon", "expression", "patient")
	allpats<- rbind(allpats, final)
}

setwd("/home/quints1/gtexstuff/L1HSpatients")
write.table(allpats, "EsophagusL1HSresults.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

#Liver
data.frame(patientlist<- read.table("/home/quints1/gtexstuff/patientnames/Liverpats.txt"))
setwd("/DataDrives/dd5/gtex/Liver")
allpats=NULL
for (i in patientlist[,1])
{
	final=NULL
	patfile=NULL
	patfile= paste(i,"/",i,".discount.ele.cntTable", sep="")
	pat<- read.table(patfile)
	rownames(pat)=pat[,1]
	L1HSpat<-pat[c("L1HS:L1:LINE"),]
	final<-cbind(L1HSpat, paste(i))
	colnames(final)=c("transposon", "expression", "patient")
	allpats<- rbind(allpats, final)
}

setwd("/home/quints1/gtexstuff/L1HSpatients")
write.table(allpats, "LiverL1HSresults.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

#Prostate
data.frame(patientlist<- read.table("/home/quints1/gtexstuff/patientnames/Prostatepats.txt"))
setwd("/DataDrives/dd5/gtex/Prostate")
allpats=NULL
for (i in patientlist[,1])
{
	final=NULL
	patfile=NULL
	patfile= paste(i,"/",i,".discount.ele.cntTable", sep="")
	pat<- read.table(patfile)
	rownames(pat)=pat[,1]
	L1HSpat<-pat[c("L1HS:L1:LINE"),]
	final<-cbind(L1HSpat, paste(i))
	colnames(final)=c("transposon", "expression", "patient")
	allpats<- rbind(allpats, final)
}

setwd("/home/quints1/gtexstuff/L1HSpatients")
write.table(allpats, "ProstateL1HSresults.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

#Thyroid
data.frame(patientlist<- read.table("/home/quints1/gtexstuff/patientnames/Thyroidpats.txt"))
setwd("/DataDrives/dd5/gtex/Thyroid")
allpats=NULL
for (i in patientlist[,1])
{
	final=NULL
	patfile=NULL
	patfile= paste(i,"/",i,".discount.ele.cntTable", sep="")
	pat<- read.table(patfile)
	rownames(pat)=pat[,1]
	L1HSpat<-pat[c("L1HS:L1:LINE"),]
	final<-cbind(L1HSpat, paste(i))
	colnames(final)=c("transposon", "expression", "patient")
	allpats<- rbind(allpats, final)
}

setwd("/home/quints1/gtexstuff/L1HSpatients")
write.table(allpats, "ThyroidL1HSresults.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

