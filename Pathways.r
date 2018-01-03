library(jsonlite)
library(org.Mm.eg.db)

### Mouse genes from Monarch
mpList <- c(
  "MP_0005426",
  "MP_0001861",
  "MP_0013806",
  "MP_0001847",
  "MP_0001429",
  "MP_0001856",
  "MP_0003606",
  "MP_0005044",
  "MP_0001953",
  "MP_0001954")
MGI <- NULL
for (MP in mpList){
  address <- paste("http://monarchinitiative.org/phenotype/", MP, "/gene_associations.json", sep="")
  test <- fromJSON(address)
  ID <- NULL
  ID <- c(ID, test$gene_associations$gene[,1])
  mouse <- which(grepl("MGI:", ID))
  MGI <- c(MGI, ID[mouse])
}

MGI <- unique(MGI)

### Convert from MGI to Entrez Gene
EG<-mget(MGI, envir=org.Mm.egMGI2EG, ifnotfound=NA)
EG<-EG[which(!is.na(EG))]



###############################
### Human influenza entrez ID

library(xlsx)
library(GEOquery)
library(org.Hs.eg.db)
library(GOstats)
library(illuminaHumanv3.db)

setwd("..")
setwd("desktop/ohsu/thesis")

### load data

GeneTable <- read.xlsx("GeneTable.xlsx", 1)

dest <- paste(getwd(), "/data", sep="")

annot1 <- getGEO("GPL6947", destdir = dest)
annot2 <- getGEO("GPL6244", destdir = dest)

a1 <- Table(annot1)
a2 <- Table(annot2)

### Select only rows with GB accesion ID

GBu1 <- a1[,1]
GBu1 <- as.character(GBu1[which(GBu1!="")])

GBu2 <- a2[,2]
GBu2 <- as.character(GBu2[which(GBu2!="")])

GBu2 <- strsplit(GBu2, ",")
GBu2l <- NULL
for (x in 1:length(GBu2)){GBu2l <- c(GBu2l, unlist(GBu2[x])[1])}
GBu2 <- GBu2l


### Select GeneTable data for each dataset

set1 <- which(GeneTable[,7]=="x")
set2 <- which(GeneTable[,8]=="x")

GTl <- strsplit(as.character(GeneTable[,4]), ",")
GT<- NULL
GT1<-as.character(GeneTable[,3])
for (x in 1:length(GTl)){GT <- c(GT, unlist(GTl[x])[1])}



GS1 <- as.character(GT1[set1])
GS2 <- as.character(GT[set2])

GS1 <- GS1[which(GS1!=" ")]
GS2 <- GS2[which(!is.na(GS2))]


#convert lists of GenBank names to Entrez IDs
entGBu1<-mget(GBu1, envir=illuminaHumanv3ENTREZREANNOTATED, ifnotfound=NA)
entGBu2<-mget(GBu2, envir=org.Hs.egACCNUM2EG, ifnotfound=NA)

entGS1<-mget(GS1, envir=illuminaHumanv3ENTREZREANNOTATED, ifnotfound=NA)
entGS2<-mget(GS2, envir=org.Hs.egACCNUM2EG, ifnotfound=NA)


entGBu1<-entGBu1[which(!is.na(entGBu1))]
entGBu2<-entGBu2[which(!is.na(entGBu2))]

entGS1<-entGS1[which(!is.na(entGS1))]
entGS2<-entGS2[which(!is.na(entGS2))]

#combine human sets
entGBu <- c(entGBu1, entGBu2)
entGSflu <- c(entGS1, entGS2)



###############################
### Human influenza entrez ID

library(xlsx)
library(GEOquery)
library(org.Hs.eg.db)
library(GOstats)
library(illuminaHumanv3.db)


### load data

GeneTable <- read.xlsx("SarsGeneTable.xlsx", 1)

dest <- paste(getwd(), "/data", sep="")

annot1 <- getGEO("GPL201", destdir = dest)
annot2 <- getGEO("GPL4387", destdir = dest)

a1 <- Table(annot1)
a2 <- Table(annot2)

### Select only rows with accesion ID

GBu1 <- a1[,2]
GBu1 <- as.character(GBu1[which(GBu1!="")])

GBu2 <- a2[,17]
GBu2 <- as.character(GBu2[which(GBu2!="")])


### Select GeneTable data for each dataset

set1 <- which(GeneTable[,6]=="x")
set2 <- which(GeneTable[,7]=="x")
GTl <- strsplit(as.character(GeneTable[,4]), ",")


GS1 <- as.character(GTl[set1])
GS2 <- as.character(GTl[set2])

GS1 <- GS1[which(GS1!=" ")]
GS1 <- unlist(strsplit(GS1, "[.]"))
GS2 <- GS2[which(!is.na(GS2))]


#convert lists of GenBank names to Entrez IDs
entGBu1<-mget(GBu1, envir=org.Hs.egACCNUM2EG, ifnotfound=NA)
entGBu2<-mget(GBu2, envir=org.Hs.egACCNUM2EG, ifnotfound=NA)

entGS1<-mget(GS1, envir=org.Hs.egACCNUM2EG, ifnotfound=NA)
entGS2<-mget(GS2, envir=org.Hs.egUNIGENE2EG, ifnotfound=NA)


entGBu1<-entGBu1[which(!is.na(entGBu1))]
entGBu2<-entGBu2[which(!is.na(entGBu2))]

entGS1<-entGS1[which(!is.na(entGS1))]
entGS2<-entGS2[which(!is.na(entGS2))]

entGBu1 <- unlist(entGBu1)
entGBu2 <- unlist(entGBu2)
entGS1 <- unlist(entGS1)
entGS2 <- unlist(entGS2)

#combine human sets
entGBusars <- c(entGBu1, entGBu2)
entGSsars <- c(entGS1, entGS2)

### variables are EG, entGSflu, entGSsars
library(KEGG.db)

mouse <- mget(as.character(EG), envir=KEGGEXTID2PATHID, ifnotfound=NA)
flu <- mget(as.character(entGSflu), envir=KEGGEXTID2PATHID, ifnotfound=NA)
sars <- mget(as.character(entGSsars), envir=KEGGEXTID2PATHID, ifnotfound=NA)

mouse<-mouse[which(!is.na(mouse))]
flu<-flu[which(!is.na(flu))]
sars<-sars[which(!is.na(sars))]

subMouse<-substr(unlist(mouse), 4, 8)
subFlu<-substr(unlist(flu), 4, 8)
subSars<-substr(unlist(sars), 4, 8)

MF <- intersect(subMouse, subFlu)
MS <- intersect(subMouse, subSars)
FS <- intersect(subFlu, subSars)

MFn <- mget(MF, envir=KEGGPATHID2NAME, ifnotfound=NA)
MSn <- mget(MS, envir=KEGGPATHID2NAME, ifnotfound=NA)
FSn <- mget(FS, envir=KEGGPATHID2NAME, ifnotfound=NA)

tableFolder<-"C:/Users/Hunter/Desktop/ohsu/thesis/tables/Pathways/"
write.csv(cbind(MF, MFn), file=paste(tableFolder, "Mouse_Flu_pathwayIntersect.csv", sep=""), row.names=F)
write.csv(cbind(MS, MSn), file=paste(tableFolder, "Mouse_SARS_pathwayIntersect.csv", sep=""), row.names=F)
write.csv(cbind(FS, FSn), file=paste(tableFolder, "Flu_SARS_pathwayIntersect.csv", sep=""), row.names=F)


### Mouse genes from Monarch
library(biomaRt)
tableFolder<-"C:/Users/Hunter/Desktop/ohsu/thesis/tables/Pathways/ByMousePheno/"

entENflu<-mget(as.character(entGSflu), envir=org.Hs.egENSEMBLPROT, ifnotfound=NA)
entENflu<-entENflu[which(!is.na(entENflu))]

entENsars<-mget(as.character(entGSsars), envir=org.Hs.egENSEMBLPROT, ifnotfound=NA)
entENsars<-entENsars[which(!is.na(entENsars))]

  humanGenes = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouseGenes = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  flulist <- getLDS(attributes = "ensembl_gene_id", filters = "ensembl_peptide_id",
    values = unlist(entENflu), mart = humanGenes,
    attributesL = c("ensembl_gene_id", "mgi_symbol", "description"), martL = mouseGenes)
  sarslist <- getLDS(attributes = "ensembl_gene_id", filters = "ensembl_peptide_id",
    values = unlist(entENsars), mart = humanGenes,
    attributesL = c("ensembl_gene_id", "mgi_symbol", "description"), martL = mouseGenes)

respIntersect <- flulist[which(flulist[,2] %in% sarslist), 2:4]

print(c('phenotype', 'flu int', 'sars int'))
mpList <- c(
  "MP_0005426",
  "MP_0001861",
  "MP_0013806",
  "MP_0001847",
  "MP_0001429",
  "MP_0001856",
  "MP_0003606",
  "MP_0005044",
  "MP_0001953",
  "MP_0001954")
for (MP in mpList){
  address <- paste("http://monarchinitiative.org/phenotype/", MP, "/gene_associations.json", sep="")
  Laddress <- paste("http://monarchinitiative.org/phenotype/", MP, "/labels.json", sep="")
  test <- fromJSON(address)
  label <- fromJSON(Laddress)
  ID <- NULL
  ID <- c(ID, test$gene_associations$gene[,1])
  Mm <- which(grepl("MGI:", ID))
  mouseGene <- ID[Mm]
  lab <- label$labels
  EG<-mget(as.character(mouseGene), envir=org.Mm.egMGI2EG, ifnotfound=NA)
  EG<-EG[which(!is.na(EG))]
  
  EN<-mget(as.character(EG), envir=org.Mm.egENSEMBL, ifnotfound=NA)
  EN<-EN[which(!is.na(EN))]
  
  mouse <- mget(as.character(EG), envir=KEGGEXTID2PATHID, ifnotfound=NA)
  mouse<-mouse[which(!is.na(mouse))]
  subMouse<-substr(unlist(mouse), 4, 8)
  MF <- intersect(subMouse, subFlu)
  MS <- intersect(subMouse, subSars)
  MFn <- mget(MF, envir=KEGGPATHID2NAME, ifnotfound=NA)
  MSn <- mget(MS, envir=KEGGPATHID2NAME, ifnotfound=NA)
  
  fluIntersect <- flulist[which(flulist[,2] %in% EN), 2:4]
  sarsIntersect <- sarslist[which(sarslist[,2] %in% EN), 2:4]
  tableFolder<-"C:/Users/Hunter/Desktop/ohsu/thesis/tables/Monarch_Human/ByMousePheno/"
  fileName <- paste(tableFolder, lab, sep="")
  write.csv(fluIntersect, file=paste(fileName, "Monarch_Flu_Gene_Intersect.csv", sep=""), row.names=F)
  write.csv(sarsIntersect, file=paste(fileName, "Monarch_SARS_Gene_Intersect.csv", sep=""), row.names=F)
  print(c(lab, dim(fluIntersect), dim(sarsIntersect)))
  
  tableFolder<-"C:/Users/Hunter/Desktop/ohsu/thesis/tables/Pathways/ByMousePheno/"
  fileName <- paste(tableFolder, lab, sep="")
  write.csv(cbind(MF, MFn), file=paste(fileName, "Mouse_Flu_pathwayIntersect.csv", sep=""), row.names=F)
  write.csv(cbind(MS, MSn), file=paste(fileName, "Mouse_SARS_pathwayIntersect.csv", sep=""), row.names=F)
}


FluSarsP <- read.csv("C:/Users/Hunter/Desktop/ohsu/thesis/tables/Pathways/Flu_SARS_pathwayIntersect.csv")
MouseFluP <- read.csv("C:/Users/Hunter/Desktop/ohsu/thesis/tables/Pathways/Mouse_Flu_pathwayIntersect.csv")
MouseSarsP <- read.csv("C:/Users/Hunter/Desktop/ohsu/thesis/tables/Pathways/Mouse_SARS_pathwayIntersect.csv")

mfsInt<-which(MouseFluP[,1] %in% MouseSarsP[,1])
totInt<-intersect(FluSarsP[,1], MouseFluP[mfsInt,1])
totInt<-paste("0",as.character(totInt), sep='')
totIntn <- mget(as.character(totInt), envir=KEGGPATHID2NAME, ifnotfound=NA)

tableFolder<-"C:/Users/Hunter/Desktop/ohsu/thesis/tables/Pathways/"
write.csv(cbind(totInt, totIntn), file=paste(tableFolder, "Mouse_Respiration_pathwayIntersect.csv", sep=""), row.names=F)

#################################################
library(jsonlite)
library(org.Mm.eg.db)

setwd("..")
setwd("desktop/ohsu/thesis")

addpt1<-'http://sirius.monarchinitiative.org:8080/solr/golr/select?defType=edismax&qt=standard&indent=on&wt=csv&rows=100000&start=0&fl=subject&facet=true&facet.mincount=1&facet.sort=count&json.nl=arrarr&facet.limit=25&csv.encapsulator=&csv.separator=%09&csv.header=false&csv.mv.separator=|&fq=object_closure:%22'
addpt2<-'%22&fq=subject_category:%22gene%22&facet.field=relation_closure_label&facet.field=evidence_closure_label&facet.field=subject_closure_label&facet.field=object_closure_label&facet.field=subject_taxon_closure_label&facet.field=object_taxon_closure_label&q=*:*'


mpList <- c(
  "MP:0005426",
  "MP:0001861",
  ###"MP:0013806",
  "MP:0001847",
  "MP:0001429",
  "MP:0001856",
  "MP:0003606",
  "MP:0005044",
  "MP:0001953",
  "MP:0001954")
MGI <- NULL
for (MP in mpList){
  address <- paste(addpt1, MP, addpt2, sep="")
  geneTable<-read.table(address)
  MPr<-gsub(":", "_", MP)
  write.table(geneTable, file=paste("data/", MPr, ".csv", sep=""), col.names=FALSE, row.names=FALSE)
  mouse <- which(grepl("MGI:", geneTable[,1]))
  MGI <- c(MGI, as.character(geneTable[mouse,]))
}

MGI <- unique(MGI)

### Convert from MGI to Entrez Gene
EG<-mget(MGI, envir=org.Mm.egMGI2EG, ifnotfound=NA)
EG<-EG[which(!is.na(EG))]

entENflu<-mget(as.character(entGSflu), envir=org.Hs.egENSEMBLPROT, ifnotfound=NA)
entENflu<-entENflu[which(!is.na(entENflu))]
entENsars<-mget(as.character(entGSsars), envir=org.Hs.egENSEMBLPROT, ifnotfound=NA)
entENsars<-entENsars[which(!is.na(entENsars))]
entENmouse<-mget(as.character(EG), envir=org.Mm.egENSEMBL, ifnotfound=NA)
entENmouse<-entENmouse[which(!is.na(entENmouse))]



library(biomaRt)
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
Mmlist <- getLDS(attributes = "ensembl_gene_id", filters = "ensembl_gene_id",
  values = unlist(entENmouse), mart = mouse,
  attributesL = c("ensembl_gene_id", "ensembl_peptide_id", "hgnc_symbol", "description"), martL = human)

totalEnsembl <- c(unlist(entENflu), unlist(entENsars), Mmlist[,3])
ute<-unique(totalEnsembl)
epByCat <- c(NULL, NULL, NULL, NULL)
for (e in ute){
  fl<-e%in%unlist(entENflu)
  sl<-e%in%unlist(entENsars)
  ml<-e%in%unlist(Mmlist[,3])
  epByCat <- rbind(epByCat, c(e, fl, sl, ml))
}


### variables are EG, entGSflu, entGSsars
library(KEGG.db)

mouse <- mget(as.character(EG), envir=KEGGEXTID2PATHID, ifnotfound=NA)
flu <- mget(as.character(entGSflu), envir=KEGGEXTID2PATHID, ifnotfound=NA)
sars <- mget(as.character(entGSsars), envir=KEGGEXTID2PATHID, ifnotfound=NA)

MmPath <- cbind(substr(unlist(mouse), 4, 8), names(unlist(mouse)))
colnames(MmPath)<-c('Ontology','Genbank')
fluPath <- cbind(substr(unlist(flu), 4, 8), names(unlist(flu)))
colnames(fluPath)<-c('Ontology','Genbank')
sarsPath <- cbind(substr(unlist(sars), 4, 8), names(unlist(sars)))
colnames(sarsPath)<-c('Ontology','Genbank')

total<-rbind(MmPath, fluPath, sarsPath)
total<-cbind(total, FALSE)
total<-cbind(total, FALSE)
total<-cbind(total, FALSE)
colnames(total)<-c('Ontology', 'Genbank', 'Mouse', 'Flu', 'SARS')

Mmt<-which(total[,1] %in% MmPath[,1])
flut<-which(total[,1] %in% fluPath[,1])
sarst<-which(total[,1] %in% sarsPath[,1])
total[Mmt, 3] <- TRUE
total[flut, 4] <- TRUE
total[sarst, 5] <- TRUE

OntSym <- mget(as.character(total[,1]), envir=KEGGPATHID2NAME, ifnotfound=NA)
GeneSym <- mget(as.character(total[,2]), envir=org.Hs.egSYMBOL, ifnotfound=NA)
mGeneSym <- mget(as.character(total[,2]), envir=org.Mm.egSYMBOL, ifnotfound=NA)

FinalTotal<-cbind(total[,1], OntSym, total[,2], mGeneSym, total[,3], total[,4], total[,5])
colnames(FinalTotal)<-c('OntologyID','Ontology Description','GenBankID','Gene Symbol','Mouse','Flu','SARS')

tableFolder<-"C:/Users/Hunter/Desktop/ohsu/thesis/tables/Pathways/"
write.csv(FinalTotal, file=paste(tableFolder, "Gene_Pathway_Intersect_all.csv", sep=""), row.names=F)
