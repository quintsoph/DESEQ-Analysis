###Get patient names lined up
setwd("/home/quints1/gtexstuff/deseqstuff")
allp<- read.table("patienttissue.txt", header=FALSE, sep="\t")
library("stringi")

#Adrenal
setwd("/home/quints1/gtexstuff/L1HSpatients")
adrenal<- read.table("AdrenalGlandL1HSresults.txt", header=TRUE, row.names=1, sep="\t")
adrenalnames<- subset(allp, allp[,2] == "Adrenal Gland", select=c(1,2))
cutL1HS<-data.frame(stri_sub(adrenal[,3],0,-15))
cutL1HS<- data.frame(cutL1HS[-121,])
adrenal<- adrenal[-121,]
adrenal2<- cbind(adrenal, cutL1HS)
colnames(adrenal2)= c("transposon", "expression", "patient", "cutname")
#patient names
cutnames<- data.frame(stri_sub(adrenalnames[,1], 0, -15))
adrenalnames2<- cbind(adrenalnames, cutnames)
colnames(adrenalnames2)= c("patient", "type", "cutname")

patientstogether<- merge(adrenal2, adrenalnames2, by.y="cutname", by.x="cutname", all=TRUE)
patientsall<- patientstogether[,c("cutname", "transposon", "expression", "patient.y", "type")]
Adrenalfinal<- na.omit(patientsall)
Adrenalfinal<- Adrenalfinal[,c(4,3,5)]

setwd("/home/quints1/gtexstuff/L1HSpatients/L1HSpatientedits")
write.table(Adrenalfinal, "Adrenaledits.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

#Bladder
setwd("/home/quints1/gtexstuff/L1HSpatients")
bladder<- read.table("BladderL1HSresults.txt", header=TRUE, row.names=1, sep="\t")
bladdernames<- subset(allp, allp[,2] == "Bladder", select=c(1,2))
cutL1HS<-data.frame(stri_sub(bladder[,3],0,-9))
bladder2<- cbind(bladder, cutL1HS)
colnames(bladder2)= c("transposon", "expression", "patient", "cutname")
#patient names
cutnames<- data.frame(stri_sub(bladdernames[,1], 0, -15))
bladdernames2<- cbind(bladdernames, cutnames)
colnames(bladdernames2)= c("patient", "type", "cutname")

patientstogether<- merge(bladder2, bladdernames2, by.y="cutname", by.x="cutname", all=TRUE)
patientsall<- patientstogether[,c("cutname", "transposon", "expression", "patient.y", "type")]
Bladderfinal<- na.omit(patientsall)
Bladderfinal<- Bladderfinal[,c(4,3,5)]

setwd("/home/quints1/gtexstuff/L1HSpatients/L1HSpatientedits")
write.table(Bladderfinal, "Bladderedits.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

#Breast
setwd("/home/quints1/gtexstuff/L1HSpatients")
breast<- read.table("BreastL1HSresults.txt", header=TRUE, row.names=1, sep="\t")
breastnames<- subset(allp, allp[,2] == "Breast - Mammary Tissue", select=c(1,2))
cutL1HS<-data.frame(stri_sub(breast[,3],0,-23))
breast2<- cbind(breast, cutL1HS)
colnames(breast2)= c("transposon", "expression", "patient", "cutname")
#patient names
cutnames<- data.frame(stri_sub(breastnames[,1], 0, -15))
breastnames2<- cbind(breastnames, cutnames)
colnames(breastnames2)= c("patient", "type", "cutname")

patientstogether<- merge(breast2, breastnames2, by.y="cutname", by.x="cutname", all=TRUE)
patientsall<- patientstogether[,c("cutname", "transposon", "expression", "patient.y", "type")]
Breastfinal<- na.omit(patientsall)
Breastfinal<- Breastfinal[,c(4,3,5)]

setwd("/home/quints1/gtexstuff/L1HSpatients/L1HSpatientedits")
write.table(Breastfinal, "Breastedits.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

#Cervix
	#Cervix - Ectocervix
	setwd("/home/quints1/gtexstuff/L1HSpatients")
	cervix<- read.table("CervixL1HSresults.txt", header=TRUE, row.names=1, sep="\t")
	cervixnames<- subset(allp, allp[,2] == "Cervix - Ectocervix", select=c(1,2))
	cutnames<- data.frame(stri_sub(cervixnames[,1], 0, -15))
	ecto<- subset(cervix, grepl("*((Ectocervix))", cervix[,3]), select=c(1,2,3))
	ectocut<- data.frame(stri_sub(ecto[,3],0,-19))
	allecto<- cbind(ecto, ectocut)
	colnames(allecto)= c("transposon", "expression", "patient", "cutname")
	cervixnames2<- cbind(cervixnames, cutnames)
	colnames(cervixnames2)= c("patient", "type", "cutname")
	cervixecto<- merge(allecto, cervixnames2, by.y="cutname", by.x="cutname", all=TRUE) #merge file based on the beginning of pat name
	cervixectos<- cervixecto[,c("cutname", "transposon", "expression", "patient.y", "type")]

	#Cervix - Endocervix
	setwd("/home/quints1/gtexstuff/L1HSpatients")
	cervix<- read.table("CervixL1HSresults.txt", header=TRUE, row.names=1, sep="\t")
	cervixnames<- subset(allp, allp[,2] == "Cervix - Endocervix", select=c(1,2))
	cutnames<- data.frame(stri_sub(cervixnames[,1], 0, -15))
	endo<- subset(cervix, grepl("*((Endocervix))", cervix[,3]), select=c(1,2,3))
	endocut<- data.frame(stri_sub(endo[,3],0,-19))
	allendo<- cbind(endo, endocut)
	colnames(allendo)= c("transposon", "expression", "patient", "cutname")
	cervixnames2<- cbind(cervixnames, cutnames)
	colnames(cervixnames2)= c("patient", "type", "cutname")
	cervixendo<- merge(allendo, cervixnames2, by.y="cutname", by.x="cutname", all=TRUE) #merge file based on the beginning of pat name
	cervixendos<- cervixendo[,c("cutname", "transposon", "expression", "patient.y", "type")]

finalcervix<- rbind(cervixectos, cervixendos)
Cervixfinal<- na.omit(finalcervix)
Cervixfinal<- Cervixfinal[,c(4,3,5)]

setwd("/home/quints1/gtexstuff/L1HSpatients/L1HSpatientedits")
write.table(Cervixfinal, "Cervixedits.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

#Colon
	#Colon - Sigmoid
	setwd("/home/quints1/gtexstuff/L1HSpatients")
	colon<- read.table("ColonL1HSresults.txt", header=TRUE, row.names=1, sep="\t")
	colonnames<- subset(allp, allp[,2] == "Colon - Sigmoid", select=c(1,2))
	cutnames<- data.frame(stri_sub(colonnames[,1], 0, -15))
	sigs<- subset(colon, grepl("*((Colon_Sigmoid))", colon[,3]), select=c(1,2,3))
	sigcut<- data.frame(stri_sub(sigs[,3],0,-15))
	allsig<- cbind(sigs, sigcut)
	colnames(allsig)= c("transposon", "expression", "patient", "cutname")
	colonnames2<- cbind(colonnames, cutnames)
	colnames(colonnames2)= c("patient", "type", "cutname")
	colonsig<- merge(allsig, colonnames2, by.y="cutname", by.x="cutname", all=TRUE) #merge file based on the beginning of pat name
	colonsigs<- colonsig[,c("cutname", "transposon", "expression", "patient.y", "type")]

	#Colon-Transverse
	setwd("/home/quints1/gtexstuff/L1HSpatients")
	colon<- read.table("ColonL1HSresults.txt", header=TRUE, row.names=1, sep="\t")
	colonnames<- subset(allp, allp[,2] == "Colon - Transverse", select=c(1,2))
	cutnames<- data.frame(stri_sub(colonnames[,1], 0, -15))
	trans<- subset(colon, grepl("*((Colon_Transverse))", colon[,3]), select=c(1,2,3))
	transcut<-data.frame(stri_sub(trans[,3],0,-18))
	alltrans<- cbind(trans, transcut)
	colnames(alltrans)=c("transposon", "expression", "patient", "cutname")
	cutnames<- data.frame(stri_sub(colonnames[,1], 0, -15))
	colonnames2<- cbind(colonnames, cutnames)
	colnames(colonnames2)= c("patient", "type", "cutname")
	colontran<- merge(alltrans, colonnames2, by.y="cutname", by.x="cutname", all=TRUE) #merge file based on the beginning of pat name
	colontrans<- colontran[,c("cutname", "transposon", "expression", "patient.y", "type")]

finalcolon<- rbind(colonsigs, colontrans)
Colonfinal<- na.omit(finalcolon)
Colonfinal<- Colonfinal[,c(4,3,5)]

setwd("/home/quints1/gtexstuff/L1HSpatients/L1HSpatientedits")
write.table(Colonfinal, "Colonedits.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

#Esophagus
	#Esophagus - Gastroesophageal
	setwd("/home/quints1/gtexstuff/L1HSpatients")
	esoph<- read.table("EsophagusL1HSresults.txt", header=TRUE, row.names=1, sep="\t")
	esophnames<- subset(allp, allp[,2] == "Esophagus - Gastroesophageal Junction", select=c(1,2)) #keep gastro tissue type only
	cutnames<- data.frame(stri_sub(esophnames[,1], 0, -15))
	gas<- subset(esoph, grepl("*((Gastroesophageal_Junction))", esoph[,3]), select=c(1,2,3)) #select gastro from L1HS file
	gascut<- data.frame(stri_sub(gas[,3],0,-37))
	allgass<- cbind(gas, gascut)
	colnames(allgass)=c("transposon", "expression", "patient", "cutname")
	esophnames2<- cbind(esophnames, cutnames)
	colnames(esophnames2)= c("patient", "type", "cutname")
	esophogas<- merge(allgass, esophnames2, by.y="cutname", by.x="cutname", all=TRUE) #merge file based on the beginning of pat name
	esophogass<- esophogas[,c("cutname", "transposon", "expression", "patient.y", "type")]

	#Esophagus - mucosa
	setwd("/home/quints1/gtexstuff/L1HSpatients")
	esoph<- read.table("EsophagusL1HSresults.txt", header=TRUE, row.names=1, sep="\t")
	esophnames<- subset(allp, allp[,2] == "Esophagus - Mucosa", select=c(1,2)) #keep mucosa tissue type only
	cutnames<- data.frame(stri_sub(esophnames[,1], 0, -15))
	mucosa<- subset(esoph, grepl("*((Mucosa))", esoph[,3]), select=c(1,2,3)) #select mucosa from L1HS file
	mucosacut<-data.frame(stri_sub(mucosa[,3],0,-18))
	allmucosa<- cbind(mucosa, mucosacut)
	colnames(allmucosa)= c("transposon", "expression", "patient", "cutname")
	esophnames2<- cbind(esophnames, cutnames)
	colnames(esophnames2)= c("patient", "type", "cutname")
	esophomucus<- merge(allmucosa, esophnames2, by.y="cutname", by.x="cutname", all=TRUE) #merge file based on the beginning of pat name
	esophomucuss<- esophomucus[,c("cutname", "transposon", "expression", "patient.y", "type")]

	#Esophagus - muscularis
	setwd("/home/quints1/gtexstuff/L1HSpatients")
	esoph<- read.table("EsophagusL1HSresults.txt", header=TRUE, row.names=1, sep="\t")
	esophnames<- subset(allp, allp[,2] == "Esophagus - Muscularis", select=c(1,2)) #keep muscularis tissue type only
	cutnames<- data.frame(stri_sub(esophnames[,1], 0, -15))
	muscular<- subset(esoph, grepl("*((Muscularis))", esoph[,3]), select=c(1,2,3))
	muscularcut<-data.frame(stri_sub(muscular[,3],0,-22))
	allmuscular<- cbind(muscular, muscularcut)
	colnames(allmuscular)= c("transposon", "expression", "patient", "cutname")
	esophnames2<- cbind(esophnames, cutnames)
	colnames(esophnames2)= c("patient", "type", "cutname")
	esophomuscular<- merge(allmuscular, esophnames2, by.y="cutname", by.x="cutname", all=TRUE) #merge file based on the beginning of pat name
	esophomuscular<- esophomuscular[,c("cutname", "transposon", "expression", "patient.y", "type")]

finalesoph<- rbind(esophogass, esophomucuss, esophomuscular)
Esophfinal<- na.omit(finalesoph)
Esophfinal<- Esophfinal[,c(4,3,5)]

setwd("/home/quints1/gtexstuff/L1HSpatients/L1HSpatientedits")
write.table(Esophfinal, "Esophagusedits.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

#Kidney
setwd("/home/quints1/gtexstuff/L1HSpatients")
kidney<- read.table("KidneyL1HSresults.txt", header=TRUE, row.names=1, sep="\t")
kidneynames<- subset(allp, allp[,2] == "Kidney - Cortex", select=c(1,2))
cutL1HS<-data.frame(stri_sub(kidney[,3],0,-15))
kidney2<- cbind(kidney, cutL1HS)
colnames(kidney2)= c("transposon", "expression", "patient", "cutname")
#patient names
cutnames<- data.frame(stri_sub(kidneynames[,1], 0, -15))
kidneynames2<- cbind(kidneynames, cutnames)
colnames(kidneynames2)= c("patient", "type", "cutname")

patientstogether<- merge(kidney2, kidneynames2, by.y="cutname", by.x="cutname", all=TRUE)
patientsall<- patientstogether[,c("cutname", "transposon", "expression", "patient.y", "type")]
Kidneyfinal<- na.omit(patientsall)
Kidneyfinal<- Kidneyfinal[,c(4,3,5)]

setwd("/home/quints1/gtexstuff/L1HSpatients/L1HSpatientedits")
write.table(Kidneyfinal, "Kidneyedits.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

#Liver	
setwd("/home/quints1/gtexstuff/L1HSpatients")
liver<- read.table("LiverL1HSresults.txt", header=TRUE, row.names=1, sep="\t")
livernames<- subset(allp, allp[,2] == "Liver", select=c(1,2))
cutL1HS<-data.frame(stri_sub(liver[,3],0,-7))
liver2<- cbind(liver, cutL1HS)
colnames(liver2)= c("transposon", "expression", "patient", "cutname")
#patient names
cutnames<- data.frame(stri_sub(livernames[,1], 0, -15))
livernames2<- cbind(livernames, cutnames)
colnames(livernames2)= c("patient", "type", "cutname")

patientstogether<- merge(liver2, livernames2, by.y="cutname", by.x="cutname", all=TRUE)
patientsall<- patientstogether[,c("cutname", "transposon", "expression", "patient.y", "type")]
Liverfinal<- na.omit(patientsall)
Liverfinal<- Liverfinal[,c(4,3,5)]

setwd("/home/quints1/gtexstuff/L1HSpatients/L1HSpatientedits")
write.table(Liverfinal, "Liveredits.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

#Pancreas
setwd("/home/quints1/gtexstuff/L1HSpatients")
pancreas<- read.table("PancreasL1HSresults.txt", header=TRUE, row.names=1, sep="\t")
pancreasnames<- subset(allp, allp[,2] == "Pancreas", select=c(1,2))
cutL1HS<-data.frame(stri_sub(pancreas[,3],0,-10))
pancreas2<- cbind(pancreas, cutL1HS)
colnames(pancreas2)= c("transposon", "expression", "patient", "cutname")
#patient names
cutnames<- data.frame(stri_sub(pancreasnames[,1], 0, -15))
pancreasnames2<- cbind(pancreasnames, cutnames)
colnames(pancreasnames2)= c("patient", "type", "cutname")

patientstogether<- merge(pancreas2, pancreasnames2, by.y="cutname", by.x="cutname", all=TRUE)
patientsall<- patientstogether[,c("cutname", "transposon", "expression", "patient.y", "type")]
Pancreasfinal<- na.omit(patientsall)
Pancreasfinal<- Pancreasfinal[,c(4,3,5)]
#WFON

setwd("/home/quints1/gtexstuff/L1HSpatients/L1HSpatientedits")
write.table(Pancreasfinal, "Pancreasedits.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

#Prostate doubles the patients for some reason #solved - randomly 
setwd("/home/quints1/gtexstuff/L1HSpatients")
prostate<- read.table("ProstateL1HSresults.txt", header=TRUE, row.names=1, sep="\t")
prostatenames<- subset(allp, allp[,2] == "Prostate", select=c(1,2))
cutL1HS<-data.frame(stri_sub(prostate[,3],0,-10))
prostate2<- cbind(prostate, cutL1HS)
colnames(prostate2)= c("transposon", "expression", "patient", "cutname")
#patient names
cutnames<- data.frame(stri_sub(prostatenames[,1], 0, -15))
prostatenames2<- cbind(prostatenames, cutnames)
colnames(prostatenames2)= c("patient", "type", "cutname")

patientstogether<- merge(prostate2, prostatenames2, by.y="cutname", by.x="cutname", all=TRUE)
patientsall<- patientstogether[,c("cutname", "transposon", "expression", "patient.y", "type")]
even <- seq_len(nrow(patientsall)) %% 2   # index
x.loadings <- data.frame(x=patientsall[!even, ])
colnames(x.loadings)= c("cutname", "transposon", "expression", "patient.y", "type")
Prostatefinal<- na.omit(x.loadings)
Prostatefinal<- Prostatefinal[,c(4,3,5)]

setwd("/home/quints1/gtexstuff/L1HSpatients/L1HSpatientedits")
write.table(Prostatefinal, "Prostateedits.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

#Stomach
setwd("/home/quints1/gtexstuff/L1HSpatients")
stomach<- read.table("StomachL1HSresults.txt", header=TRUE, row.names=1, sep="\t")
stomachnames<- subset(allp, allp[,2] == "Stomach", select=c(1,2))
cutL1HS<-data.frame(stri_sub(stomach[,3],0,-9))
stomach2<- cbind(stomach, cutL1HS)
colnames(stomach2)= c("transposon", "expression", "patient", "cutname")
#patient names
cutnames<- data.frame(stri_sub(stomachnames[,1], 0, -15))
stomachnames2<- cbind(stomachnames, cutnames)
colnames(stomachnames2)= c("patient", "type", "cutname")

patientstogether<- merge(stomachnames2, stomach2, by.y="cutname", by.x="cutname", all=TRUE)
patientsall<- patientstogether[,c("cutname", "transposon", "expression", "patient.x", "type")]
Stomachfinal<- na.omit(patientsall)
Stomachfinal<- Stomachfinal[,c(4,3,5)]

setwd("/home/quints1/gtexstuff/L1HSpatients/L1HSpatientedits")
write.table(Stomachfinal, "Stomachedits.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

#Thyroid
setwd("/home/quints1/gtexstuff/L1HSpatients")
thyroid<- read.table("ThyroidL1HSresults.txt", header=TRUE, row.names=1, sep="\t")
thyroidnames<- subset(allp, allp[,2] == "Thyroid", select=c(1,2))
cutL1HS<-data.frame(stri_sub(thyroid[,3],0,-9))
thyroid2<- cbind(thyroid, cutL1HS)
colnames(thyroid2)= c("transposon", "expression", "patient", "cutname")
#patient names
cutnames<- data.frame(stri_sub(thyroidnames[,1], 0, -15))
thyroidnames2<- cbind(thyroidnames, cutnames)
colnames(thyroidnames2)= c("patient", "type", "cutname")

patientstogether<- merge(thyroidnames2, thyroid2, by.y="cutname", by.x="cutname", all=TRUE)
patientsall<- patientstogether[,c("cutname", "transposon", "expression", "patient.x", "type")]
Thyroidfinal<- na.omit(patientsall)
Thyroidfinal<- Thyroidfinal[,c(4,3,5)]

setwd("/home/quints1/gtexstuff/L1HSpatients/L1HSpatientedits")
write.table(Thyroidfinal, "Thyroidedits.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
#need to manually edit out

##Create final matrix with all tissue types
setwd("/home/quints1/gtexstuff/L1HSpatients/L1HSpatientedits")
adrenal<- read.table("Adrenaledits.txt", header=TRUE, sep="\t")
colnames(adrenal)= c("patient", "expression", "type")
bladder<- read.table("Bladderedits.txt", header=TRUE, sep="\t")
colnames(bladder)= c("patient", "expression", "type")
breast<- read.table("Breastedits.txt", header=TRUE, sep="\t")
colnames(breast)= c("patient", "expression", "type")
cervix<- read.table("Cervixedits.txt", header=TRUE, sep="\t")
colnames(cervix)= c("patient", "expression", "type")
colon<- read.table("Colonedits.txt", header=TRUE, sep="\t")
colnames(colon)= c("patient", "expression", "type")
esoph<- read.table("Esophagusedits.txt", header=TRUE, sep="\t")
colnames(esoph)= c("patient", "expression", "type")
kidney<- read.table("Kidneyedits.txt", header=TRUE, sep="\t")
colnames(kidney)= c("patient", "expression", "type")
liver<- read.table("Liveredits.txt", header=TRUE, sep="\t")
colnames(liver)= c("patient", "expression", "type")
pancreas<- read.table("Pancreasedits.txt", header=TRUE, sep="\t")
colnames(pancreas)= c("patient", "expression", "type")
prostate<- read.table("Prostateedits.txt", header=TRUE, sep="\t")
colnames(prostate)= c("patient", "expression", "type")
stomach<- read.table("Stomachedits.txt", header=TRUE, sep="\t")
colnames(stomach)= c("patient", "expression", "type")
thyroid<- read.table("Thyroidedits.txt", header=TRUE, sep="\t")
colnames(thyroid)= c("patient", "expression", "type")

alltissue<- rbind(adrenal,bladder,breast,cervix,colon,esoph,kidney,
				liver,pancreas,prostate,stomach,thyroid)
write.table(alltissue, "allL1HSmatrix.txt", col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
tissueall<- alltissue[,c(1,2)]
colnames(tissueall)=c("patients", "L1HS")
write.table(tissueall, "notissueL1HS.txt", col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
totaltissue<- read.table("notissueL1HS.txt", header=TRUE, row.names=1, sep="\t")

#make table for proposal
L1HSpatients<- read.table("allL1HSmatrix.txt", header=TRUE, sep="\t")
adrenal<- subset(L1HSpatients, L1HSpatients[,3] == "Adrenal Gland", select=c(1,2,3))
adrenalmean<- mean(adrenal[,2])
adrenalmedian<- median(adrenal[,2])
bladder<- subset(L1HSpatients, L1HSpatients[,3] == "Bladder", select=c(1,2,3))
bladdermean<- mean(bladder[,2])
bladdermedian<- median(bladder[,2])
breast<- subset(L1HSpatients, L1HSpatients[,3] == "Breast - Mammary Tissue", select=c(1,2,3))
breastmean<- mean(breast[,2])
breastmedian<- median(breast[,2])
cervixecto<- subset(L1HSpatients, L1HSpatients[,3] == "Cervix - Ectocervix", select=c(1,2,3))
cervixectomean<- mean(cervixecto[,2])
cervixectomedian<- median(cervixecto[,2])
cervixendo<- subset(L1HSpatients, L1HSpatients[,3] == "Cervix - Endocervix", select=c(1,2,3))
cervixendomean<- mean(cervixendo[,2])
cervixendomedian<- median(cervixendo[,2])
colonsig<- subset(L1HSpatients, L1HSpatients[,3] == "Colon - Sigmoid", select=c(1,2,3))
colonsigmean<- mean(colonsig[,2])
colonsigmedian<- median(colonsig[,2])
colontrans<- subset(L1HSpatients, L1HSpatients[,3] == "Colon - Transverse", select=c(1,2,3))
colontransmean<- mean(colontrans[,2])
colontransmedian<- median(cerivixendo[,2])
esophgast<- subset(L1HSpatients, L1HSpatients[,3] == "Esophagus - Gastroesophageal Junction", select=c(1,2,3))
esophgastmean<- mean(esophgast[,2])
esophgastmedian<- median(esophgast[,2])
esophmucosa<-subset(L1HSpatients, L1HSpatients[,3] == "Esophagus - Mucosa", select=c(1,2,3))
esophmucosamean<- mean(esophmucosa[,2])
esophmucosamedian<- median(esophmucosa[,2])
esophmuscular<- subset(L1HSpatients, L1HSpatients[,3] == "Esophagus - Muscularis", select=c(1,2,3))
esophmuscularmean<- mean(esophmuscular[,2])
esophmuscularmedian<- median(esophmuscular[,2])
kidney<- subset(L1HSpatients, L1HSpatients[,3] == "Kidney - Cortex", select=c(1,2,3))
kidneymean<- mean(kidney[,2])
kidneymedian<- median(kidney[,2])
liver<- subset(L1HSpatients, L1HSpatients[,3] == "Liver", select=c(1,2,3))
livermean<- mean(liver[,2])
livermedian<- median(liver[,2])
pancreas<- subset(L1HSpatients, L1HSpatients[,3] == "Pancreas", select=c(1,2,3))
pancreasmean<- mean(pancreas[,2])
pancreasmedian<- median(pancreas[,2])
prostate<- subset(L1HSpatients, L1HSpatients[,3] == "Prostate", select=c(1,2,3))
prostatemean<- mean(prostate[,2])
prostatemedian<- median(prostate[,2])
stomach<- subset(L1HSpatients, L1HSpatients[,3] == "Stomach", select=c(1,2,3))
stomachmean<- mean(stomach[,2])
stomachmedian<- median(stomach[,2])
thyroid<- subset(L1HSpatients, L1HSpatients[,3] == "Thyroid", select=c(1,2,3))
thyroidmean<- mean(thyroid[,2])
thyroidmedian<- median(thyroid[,2])
meanL1HS<- data.frame(adrenalmean, bladdermean, breastmean, cervixectomean, cervixendomean,
colonsigmean, colontransmean, esophgastmean, esophmucosamean, esophmuscularmean,
kidneymean, livermean, pancreasmean, prostatemean, stomachmean, thyroidmean)
colnames(meanL1HS)= c("Adrenal Gland", "Bladder", "Breast Mammary Tissue", "Cervix - Ectocervix",
"Cervix - Endocervix", "Colon - Sigmoid", "Colon - Transverse", "Esophagus - Gastroesophageal",
"Esophagus - Mucosa", "Esophagus - Muscularis", "Kidney", "Liver", "Pancreas", "Prostate", "Stomach", "Thyroid")
flipmean<- t(meanL1HS)
flipmean<- data.frame(rownames(flipmean), flipmean)

medianL1HS<- data.frame(adrenalmedian, bladdermedian, breastmedian, cervixectomedian, cervixendomedian,
colonsigmedian, colontransmedian, esophgastmedian, esophmucosamedian, esophmuscularmedian,
kidneymedian, livermedian, pancreasmedian, prostatemedian, stomachmedian, thyroidmedian)
colnames(medianL1HS)= c("Adrenal Gland", "Bladder", "Breast Mammary Tissue", "Cervix - Ectocervix",
"Cervix - Endocervix", "Colon - Sigmoid", "Colon - Transverse", "Esophagus - Gastroesophageal",
"Esophagus - Mucosa", "Esophagus - Muscularis", "Kidney", "Liver", "Pancreas", "Prostate", "Stomach", "Thyroid")
flipmedian<- t(medianL1HS)
flipmedian<- data.frame(rownames(flipmedian), flipmedian)

library("ggplot2")
L1HSmeangraph <- data.frame(Tissue= flipmean[,1], 
            Mean = flipmean[,2])
ggplot(L1HSmeangraph, aes(x=Tissue,y = Mean, fill=Tissue)) +geom_bar(stat = "identity") +
 
scale_fill_hue(c=45, l=80)+ theme_classic() + ggtitle("Mean L1HS Expression levels Across Tissue Types")+ 
theme(plot.title = element_text(face="bold"), axis.text.x = element_text(angle = 90, hjust = 1))

L1HSmediangraph <- data.frame(Tissue= flipmedian[,1], 
            Median = flipmedian[,2])
ggplot(L1HSmediangraph, aes(x=Tissue,y = Median, fill=Tissue)) +geom_bar(stat = "identity") +
 
scale_fill_hue(c=45, l=80)+ theme_classic() + ggtitle("Median L1HS Expression levels Across Tissue Types")+ 
theme(plot.title = element_text(face="bold"), axis.text.x = element_text(angle = 90, hjust = 1))