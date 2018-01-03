library(limma)
library(affy)
library(oligo)

### set directory and folders
setwd("C:/Users/Hunter/Desktop/ohsu/Thesis/data/Mouse")

FigureFolder<-"C:/Users/Hunter/Desktop/ohsu/Thesis/Figures/MouseExpression/"


### declare CEL files and read into R
target1<-c(
"16441x8024_f106_d7.CEL",
"16441x8024_f112_d7.CEL",
"16441x8024_f113_d7.CEL",
"16441x8024_f114_d7.CEL",
"16441x8024_f78_d4m.CEL",
"16441x8024_f79_d4m.CEL",
"16441x8024_f80_d4m.CEL")

target8<-c(
"8048x8026_23_f_Lu_H_IV_D7.CEL",
"8048x8026_28_f_Lu_H_M_D4.CEL",
"8048x8026_37_f_Lu_H_M_D4.CEL",
"8048x8026_38_f_Lu_H_M_D4.CEL",
"8048x8026_59_f_Lu_H_IV_D7.CEL",
"8048x8026_60_f_Lu_H_IV_D7.CEL")

rawData1<-read.celfiles(target1)
rawData8<-read.celfiles(target8)


### rma normalization

normData1<-rma(rawData1)
normData8<-rma(rawData8)


### unnormalized and normalized MA plots

par(mfrow=c(3,3))
MAplot(rawData1, main="")
dev.print(pdf, paste(FigureFolder,"16441x8024_MA_noNormalize.pdf", sep=""))

par(mfrow=c(3,3))
MAplot(rawData8, main="")
dev.print(pdf, paste(FigureFolder,"8048x8026_MA_noNormalize.pdf", sep=""))

par(mfrow=c(3,3))
MAplot(normData1, main="normalized")
dev.print(pdf, paste(FigureFolder,"16441x8024_MA_Normalize.pdf", sep=""))

par(mfrow=c(3,3))
MAplot(normData8, main="normalized")
dev.print(pdf, paste(FigureFolder,"8048x8026_MA_Normalize.pdf", sep=""))

dev.off()


### unnormalized and normalized box plots
color1<-c("blue", "blue", "blue", "blue", "blue", "blue", "blue")
color8<-c("red", "red", "red", "red", "red", "red")

color1[which(grepl("d7", sampleNames(rawData1)))]<-"red"
color8[which(grepl("M_D4", sampleNames(rawData8)))]<-"blue"

boxplot(rawData1, col=color1)
dev.print(pdf, paste(FigureFolder,"16441x8024_boxplot_noNormalize.pdf", sep=""))

boxplot(normData1, col=color1)
dev.print(pdf, paste(FigureFolder,"16441x8024_boxplot_Normalize.pdf", sep=""))

boxplot(rawData8, col=color8)
dev.print(pdf, paste(FigureFolder,"8048x8026_boxplot_noNormalize.pdf", sep=""))

boxplot(normData8, col=color8)
dev.print(pdf, paste(FigureFolder,"8048x8026_boxplot_Normalize.pdf", sep=""))


### create design
Cy5<-as.factor(color1)
levels(Cy5)<-c("control", "infected")
Cy5<-as.character(Cy5)
Cy3<-"ref"
design1<-modelMatrix(cbind(target1,Cy5,Cy3),ref="ref")
cont.matrix1<-makeContrasts(infected-control, levels=design1)

Cy5<-as.factor(color8)
levels(Cy5)<-c("control", "infected")
Cy5<-as.character(Cy5)
Cy3<-"ref"
design8<-modelMatrix(cbind(target8,Cy5,Cy3),ref="ref")
cont.matrix8<-makeContrasts(infected-control, levels=design8)

### Create linear model and bayes
fit<-lmFit(normData1, design1)
fit1<-contrasts.fit(fit, cont.matrix1)
fit1.bayes<-eBayes(fit1)

fit<-lmFit(normData8, design8)
fit8<-contrasts.fit(fit, cont.matrix8)
fit8.bayes<-eBayes(fit8)


#create table of differential expression
#adjusted p-value cutoff of 0.05
#minimum fold change of 1.5
top1<-topTable(fit1.bayes, p.value=0.05, lfc=log2(1.5), number=30000)
results1<-decideTests(fit1.bayes, p.value=0.05, lfc=log2(1.5))

top8<-topTable(fit8.bayes, p.value=0.05, lfc=log2(1.5), number=30000)
results8<-decideTests(fit8.bayes, p.value=0.05, lfc=log2(1.5))

### transcript IDs
trans1 <- rownames(top1)
trans8 <- rownames(top8)
transInt <- intersect(trans1, trans8)

annotTable <- read.csv("MoGene-2_1-st-v1.na35.mm10.transcript.csv", skip=22)

which1 <- which(testTable[,1] %in% trans1)
which8 <- which(testTable[,1] %in% trans8)
whichInt <- which(testTable[,1] %in% transInt)

geneAssign1 <- testTable[which1,8]
geneAssign8 <- testTable[which8,8]
geneAssignInt <- testTable[whichInt,8]


### save results
tableFolder <- "C:/Users/Hunter/Desktop/ohsu/Thesis/tables/MouseExpression/"

write.csv(geneAssign1, file=paste(tableFolder, "16441x8024_DElist.csv", sep=""), row.names=F)
write.csv(geneAssign8, file=paste(tableFolder, "8048x8026_DElist.csv", sep=""), row.names=F)
write.csv(geneAssignInt, file=paste(tableFolder, "intersect_DElist.csv", sep=""), row.names=F)

write.csv(top1, file=paste(tableFolder, "16441x8024_DEValueTable.csv", sep=""))
write.csv(top8, file=paste(tableFolder, "8048x8026_DEValueTable.csv", sep=""))
write.csv(geneAssignInt, file=paste(tableFolder, "intersect_DElist.csv", sep=""), row.names=F)


write.csv(trans1, file=paste(tableFolder, "16441x8024_DElist_TranscriptID.csv", sep=""), row.names=F)
write.csv(trans8, file=paste(tableFolder, "8048x8026_DElist_TranscriptID.csv", sep=""), row.names=F)
write.csv(transInt, file=paste(tableFolder, "intersect_DElist_TranscriptID.csv", sep=""), row.names=F)

library(mogene21sttranscriptcluster.db)

EG1<-mget(trans1, envir=mogene21sttranscriptclusterENSEMBL, ifnotfound=NA)
EG8<-mget(trans8, envir=mogene21sttranscriptclusterENSEMBL, ifnotfound=NA)
EGI<-mget(transInt, envir=mogene21sttranscriptclusterENSEMBL, ifnotfound=NA)

EG1<-EG1[which(!is.na(EG1))]
EG8<-EG8[which(!is.na(EG8))]
EGI<-EGI[which(!is.na(EGI))]
