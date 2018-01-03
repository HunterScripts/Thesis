library(plethy)

file.name <- "C:/Users/Hunter/Desktop/ohsu/Thesis/data/March_2014i buxco.txt"
chunk.size <- 500
db.name <- file.path("C:/Users/Hunter/Desktop/ohsu/Thesis/data", "bux_test.db")
parse.buxco(file.name=file.name, chunk.size=chunk.size, db.name=db.name, verbose=FALSE)
bux.db <- makeBuxcoDB(db.name=db.name)
addAnnotation(bux.db, query=day.infer.query, index=FALSE)
addAnnotation(bux.db, query=break.type.query, index=TRUE)


namesTable <- matrix(unlist(strsplit(samples(bux.db), " ")), ncol=3, byrow=TRUE)
samps <- unique(namesTable[,1])

transl.table <- data.frame ( sample_status = c( "mock" , "sars" , "flu" ),
  Inf_Status = c( "Mock" , "SARS" , "Flu" ), stringsAsFactors = F)

use.dta <- data.frame( samples = samples ( bux.db ), sample_status = 
  sapply ( strsplit ( samples ( bux.db ) , " "), "[", 3) , stringsAsFactors = F )


transl.table2 <- data.frame ( sample_name = samps,
  Sample = samps, stringsAsFactors = F)

use.dta2 <- data.frame( samples = samples ( bux.db ), sample_strain = 
  sapply ( strsplit ( samples ( bux.db ) , " "), "[", 1) , stringsAsFactors = F )


merge.dta <- merge ( use.dta , transl.table , by = "sample_status" )
merge.dta <- merge.dta [ , -which ( names ( merge.dta ) == "sample_status" ) ]



add.labels.by.sample ( bux.db , merge.dta )
add.labels.by.sample ( bux.db , use.dta2 )
annoLevels ( bux.db )


mvtsplot ( bux.db , plot.value = "f" , Break_type_label = "EXP" , outer.group.name = "Inf_Status" ,
  subset = bux.db$sample_strain == samps[1], samples = which(), outer.cols = c( Flu = "black" , SARS = "brown" , Mock = "blue" ) )

mvtsplot ( bux.db , plot.value = "Penh" , Break_type_label = "EXP" , outer.group.name = "sample_strain" ,
  outer.cols = c( Flu = "black" , SARS = "brown" , Mock = "blue" ) )


fileAdd <- "C:/Users/Hunter/Desktop/ohsu/Thesis/Figures/BuxcoBoxplots/"

for (i in 1:length(samps)){
  mvtsplot ( bux.db , plot.value = "f" , Break_type_label = "EXP" , outer.group.name = "Inf_Status" ,
    sample_strain = samps[i], outer.cols = c( Flu = "black" , SARS = "brown" , Mock = "blue" ))
  dev.print(pdf, paste(fileAdd, samps[i], "_FreqMVTSplot", ".pdf", sep=""))
  dev.off()

  mvtsplot ( bux.db , plot.value = "Penh" , Break_type_label = "EXP" , outer.group.name = "Inf_Status" ,
    sample_strain = samps[i], outer.cols = c( Flu = "black" , SARS = "brown" , Mock = "blue" ))
  dev.print(pdf, paste(fileAdd, samps[i], "_PenhMVTSplot", ".pdf", sep=""))
  dev.off()
}


exp.f <- retrieveData(bux.db, variables="f", Break_type_label ='EXP')
exp.penh <- retrieveData(bux.db, variables="Penh", Break_type_label ='EXP')

days <- unique(exp.f$Days)
status <- transl.table[,1]
table <- NULL
annotCol <- NULL

for (i in 1:length(samps)){
  for (j in 1:length(status)){
    MeanList <- NULL
    annotCol <- rbind(annotCol, c(samps[i], status[j], "mean"))
    for (k in 1:length(days)){
      subset<-which(grepl(samps[i],exp.f$Sample_Name) & grepl(status[j], exp.f$Sample_Name) & exp.f$Days==days[k])
      ave <- mean(exp.f$Value[subset], na.rm = TRUE)
      MeanList <- c(MeanList, ave)
    }
    table <- rbind(table, MeanList)
  }
  mock <- nrow(table) - 2
  sars <- nrow(table) - 1
  flu <- nrow(table)
  table <- rbind(table, diff(rbind(table[mock,],table[sars,])))
  table <- rbind(table, table[sars,]/table[mock,])
  table <- rbind(table, diff(rbind(table[mock,],table[flu,])))
  table <- rbind(table, table[flu,]/table[mock,])
  annotCol <- rbind(annotCol, c(samps[i], status[2], "Diff"))
  annotCol <- rbind(annotCol, c(samps[i], status[2], "Fold"))
  annotCol <- rbind(annotCol, c(samps[i], status[3], "Diff"))
  annotCol <- rbind(annotCol, c(samps[i], status[3], "Fold"))
  
}

CombTable <- cbind(annotCol, table)
colnames(CombTable) <- c("Sample", "Status", "Metric", paste("day", days))

file <- "C:/Users/Hunter/Desktop/ohsu/Thesis/tables/March_2014i_FreqSummaryTable.csv"
write.csv(CombTable, file, row.names = FALSE)

### Penh

days <- unique(exp.penh$Days)
status <- transl.table[,1]
table <- NULL
annotCol <- NULL

for (i in 1:length(samps)){
  for (j in 1:length(status)){
    MeanList <- NULL
    annotCol <- rbind(annotCol, c(samps[i], status[j], "mean"))
    for (k in 1:length(days)){
      subset<-which(grepl(samps[i],exp.penh$Sample_Name) & grepl(status[j], exp.penh$Sample_Name) & exp.penh$Days==days[k])
      ave <- mean(exp.penh$Value[subset], na.rm = TRUE)
      MeanList <- c(MeanList, ave)
    }
    table <- rbind(table, MeanList)
  }
  mock <- nrow(table) - 2
  sars <- nrow(table) - 1
  flu <- nrow(table)
  table <- rbind(table, diff(rbind(table[mock,],table[sars,])))
  table <- rbind(table, table[sars,]/table[mock,])
  table <- rbind(table, diff(rbind(table[mock,],table[flu,])))
  table <- rbind(table, table[flu,]/table[mock,])
  annotCol <- rbind(annotCol, c(samps[i], status[2], "Diff"))
  annotCol <- rbind(annotCol, c(samps[i], status[2], "Fold"))
  annotCol <- rbind(annotCol, c(samps[i], status[3], "Diff"))
  annotCol <- rbind(annotCol, c(samps[i], status[3], "Fold"))
  
}

CombTable <- cbind(annotCol, table)
colnames(CombTable) <- c("Sample", "Status", "Metric", paste("day", days))

file <- "C:/Users/Hunter/Desktop/ohsu/Thesis/tables/March_2014i_PenhSummaryTable.csv"
write.csv(CombTable, file, row.names = FALSE)


###########################################################


file.name <- "C:/Users/Hunter/Desktop/ohsu/Thesis/data/August_2014i buxco.txt"

chunk.size <- 500
db.name <- file.path("C:/Users/Hunter/Desktop/ohsu/Thesis/data", "bux_test2.db")
parse.buxco(file.name=file.name, chunk.size=chunk.size, db.name=db.name, verbose=FALSE)
bux.db <- makeBuxcoDB(db.name=db.name)
addAnnotation(bux.db, query=day.infer.query, index=FALSE)
addAnnotation(bux.db, query=break.type.query, index=TRUE)


namesTable <- matrix(unlist(strsplit(samples(bux.db), " ")), ncol=3, byrow=TRUE)
samps <- unique(namesTable[,1])


transl.table <- data.frame ( sample_status = c( "mock" , "sars", "ssars", "flu" ),
  Inf_Status = c( "Mock" , "SARS" , "SARS" , "Flu" ), stringsAsFactors = F)
  
use.dta <- data.frame( samples = samples ( bux.db ), sample_status = 
  sapply ( strsplit ( samples ( bux.db ) , " "), "[", 3) , stringsAsFactors = F )


transl.table2 <- data.frame ( sample_name = samps,
  Sample = samps, stringsAsFactors = F)

use.dta2 <- data.frame( samples = samples ( bux.db ), sample_strain = 
  sapply ( strsplit ( samples ( bux.db ) , " "), "[", 1) , stringsAsFactors = F )


merge.dta <- merge ( use.dta , transl.table , by = "sample_status" )
merge.dta <- merge.dta [ , -which ( names ( merge.dta ) == "sample_status" ) ]



add.labels.by.sample ( bux.db , merge.dta )
add.labels.by.sample ( bux.db , use.dta2 )
annoLevels ( bux.db )

transl.table <- data.frame ( sample_status = c( "mock" , "sars" , "flu" ),
  Inf_Status = c( "Mock" , "SARS" , "Flu" ), stringsAsFactors = F)



mvtsplot ( bux.db , plot.value = "f" , Break_type_label = "EXP" , outer.group.name = "sample_strain" ,
  outer.cols = c( "16012x15119" = "black" , "16072x3393" = "brown" , "16557x13067" = "blue" ))


fileAdd <- "C:/Users/Hunter/Desktop/ohsu/Thesis/Figures/BuxcoBoxplots/"

dev.print(pdf, paste(fileAdd, "August_2014i", "_FreqMVTSplot", ".pdf", sep=""))
dev.off()

mvtsplot ( bux.db , plot.value = "Penh" , Break_type_label = "EXP" , outer.group.name = "sample_strain" ,
  outer.cols = c( "16012x15119" = "black" , "16072x3393" = "brown" , "16557x13067" = "blue" ))

dev.print(pdf, paste(fileAdd, "August_2014i", "_PenhMVTSplot", ".pdf", sep=""))
dev.off()

exp.f <- retrieveData(bux.db, variables="f", Break_type_label ='EXP')
exp.penh <- retrieveData(bux.db, variables="Penh", Break_type_label ='EXP')

days <- unique(exp.f$Days)
status <- transl.table[,1]
table <- NULL
annotCol <- NULL

for (i in 1:length(samps)){
  for (j in 1:length(status)){
    MeanList <- NULL
    annotCol <- rbind(annotCol, c(samps[i], status[j], "mean"))
    for (k in 1:length(days)){
      subset<-which(grepl(samps[i],exp.f$Sample_Name) & grepl(status[j], exp.f$Sample_Name) & exp.f$Days==days[k])
      ave <- mean(exp.f$Value[subset], na.rm = TRUE)
      MeanList <- c(MeanList, ave)
    }
    table <- rbind(table, MeanList)
  }
  mock <- nrow(table) - 2
  sars <- nrow(table) - 1
  flu <- nrow(table)
  table <- rbind(table, diff(rbind(table[mock,],table[sars,])))
  table <- rbind(table, table[sars,]/table[mock,])
  table <- rbind(table, diff(rbind(table[mock,],table[flu,])))
  table <- rbind(table, table[flu,]/table[mock,])
  annotCol <- rbind(annotCol, c(samps[i], status[2], "Diff"))
  annotCol <- rbind(annotCol, c(samps[i], status[2], "Fold"))
  annotCol <- rbind(annotCol, c(samps[i], status[3], "Diff"))
  annotCol <- rbind(annotCol, c(samps[i], status[3], "Fold"))
  
}

CombTable <- cbind(annotCol, table)
colnames(CombTable) <- c("Sample", "Status", "Metric", paste("day", days))

file <- "C:/Users/Hunter/Desktop/ohsu/Thesis/tables/August_2014i_FreqSummaryTable.csv"
write.csv(CombTable, file, row.names = FALSE)

### Penh

days <- unique(exp.penh$Days)
status <- transl.table[,1]
table <- NULL
annotCol <- NULL

for (i in 1:length(samps)){
  for (j in 1:length(status)){
    MeanList <- NULL
    annotCol <- rbind(annotCol, c(samps[i], status[j], "mean"))
    for (k in 1:length(days)){
      subset<-which(grepl(samps[i],exp.penh$Sample_Name) & grepl(status[j], exp.penh$Sample_Name) & exp.penh$Days==days[k])
      ave <- mean(exp.penh$Value[subset], na.rm = TRUE)
      MeanList <- c(MeanList, ave)
    }
    table <- rbind(table, MeanList)
  }
  mock <- nrow(table) - 2
  sars <- nrow(table) - 1
  flu <- nrow(table)
  table <- rbind(table, diff(rbind(table[mock,],table[sars,])))
  table <- rbind(table, table[sars,]/table[mock,])
  table <- rbind(table, diff(rbind(table[mock,],table[flu,])))
  table <- rbind(table, table[flu,]/table[mock,])
  annotCol <- rbind(annotCol, c(samps[i], status[2], "Diff"))
  annotCol <- rbind(annotCol, c(samps[i], status[2], "Fold"))
  annotCol <- rbind(annotCol, c(samps[i], status[3], "Diff"))
  annotCol <- rbind(annotCol, c(samps[i], status[3], "Fold"))
  
}

CombTable <- cbind(annotCol, table)
colnames(CombTable) <- c("Sample", "Status", "Metric", paste("day", days))

file <- "C:/Users/Hunter/Desktop/ohsu/Thesis/tables/August_2014i_PenhSummaryTable.csv"
write.csv(CombTable, file, row.names = FALSE)




##################################


file.name <- "C:/Users/Hunter/Desktop/ohsu/Thesis/data/May_2014i buxco.txt"

chunk.size <- 500
db.name <- file.path("C:/Users/Hunter/Desktop/ohsu/Thesis/data", "bux_test3.db")
parse.buxco(file.name=file.name, chunk.size=chunk.size, db.name=db.name, verbose=FALSE)
bux.db <- makeBuxcoDB(db.name=db.name)
addAnnotation(bux.db, query=day.infer.query, index=FALSE)
addAnnotation(bux.db, query=break.type.query, index=TRUE)


namesTable <- matrix(unlist(strsplit(samples(bux.db), " ")), ncol=3, byrow=TRUE)
samps <- unique(namesTable[,1])

transl.table <- data.frame ( sample_status = c( "mock" , "sars" , "flu" ),
  Inf_Status = c( "Mock" , "SARS" , "Flu" ), stringsAsFactors = F)

use.dta <- data.frame( samples = samples ( bux.db ), sample_status = 
  sapply ( strsplit ( samples ( bux.db ) , " "), "[", 3) , stringsAsFactors = F )


transl.table2 <- data.frame ( sample_name = samps,
  Sample = samps, stringsAsFactors = F)

use.dta2 <- data.frame( samples = samples ( bux.db ), sample_strain = 
  sapply ( strsplit ( samples ( bux.db ) , " "), "[", 1) , stringsAsFactors = F )


merge.dta <- merge ( use.dta , transl.table , by = "sample_status" )
merge.dta <- merge.dta [ , -which ( names ( merge.dta ) == "sample_status" ) ]



add.labels.by.sample ( bux.db , merge.dta )
add.labels.by.sample ( bux.db , use.dta2 )
annoLevels ( bux.db )


  mvtsplot ( bux.db , plot.value = "f" , Break_type_label = "EXP" , outer.group.name = "Inf_Status" ,
    sample_strain = samps[1], outer.cols = c( Flu = "black" , SARS = "brown" , Mock = "blue" ))


fileAdd <- "C:/Users/Hunter/Desktop/ohsu/Thesis/Figures/BuxcoBoxplots/"

mvtsplot ( bux.db , plot.value = "f" , Break_type_label = "EXP" , outer.group.name = "sample_strain" ,
  outer.cols = c( "16441x8024" = "black" , "5489x16557" = "brown"))

dev.print(pdf, paste(fileAdd, "May_2014i", "_FreqMVTSplot", ".pdf", sep=""))
dev.off()


mvtsplot ( bux.db , plot.value = "Penh" , Break_type_label = "EXP" , outer.group.name = "sample_strain" ,
  outer.cols = c( "16441x8024" = "black" , "5489x16557" = "brown" ))

dev.print(pdf, paste(fileAdd, "May_2014i", "_PenhMVTSplot", ".pdf", sep=""))
dev.off()


exp.f <- retrieveData(bux.db, variables="f", Break_type_label ='EXP')
exp.penh <- retrieveData(bux.db, variables="Penh", Break_type_label ='EXP')

days <- unique(exp.f$Days)
status <- transl.table[,1]
table <- NULL
annotCol <- NULL

for (i in 1:length(samps)){
  for (j in 1:length(status)){
    MeanList <- NULL
    annotCol <- rbind(annotCol, c(samps[i], status[j], "mean"))
    for (k in 1:length(days)){
      subset<-which(grepl(samps[i],exp.f$Sample_Name) & grepl(status[j], exp.f$Sample_Name) & exp.f$Days==days[k])
      ave <- mean(exp.f$Value[subset], na.rm = TRUE)
      MeanList <- c(MeanList, ave)
    }
    table <- rbind(table, MeanList)
  }
  mock <- nrow(table) - 2
  sars <- nrow(table) - 1
  flu <- nrow(table)
  table <- rbind(table, diff(rbind(table[mock,],table[sars,])))
  table <- rbind(table, table[sars,]/table[mock,])
  table <- rbind(table, diff(rbind(table[mock,],table[flu,])))
  table <- rbind(table, table[flu,]/table[mock,])
  annotCol <- rbind(annotCol, c(samps[i], status[2], "Diff"))
  annotCol <- rbind(annotCol, c(samps[i], status[2], "Fold"))
  annotCol <- rbind(annotCol, c(samps[i], status[3], "Diff"))
  annotCol <- rbind(annotCol, c(samps[i], status[3], "Fold"))
  
}

CombTable <- cbind(annotCol, table)
colnames(CombTable) <- c("Sample", "Status", "Metric", paste("day", days))

file <- "C:/Users/Hunter/Desktop/ohsu/Thesis/tables/May_2014i_FreqSummaryTable.csv"
write.csv(CombTable, file, row.names = FALSE)
### Penh

days <- unique(exp.penh$Days)
status <- transl.table[,1]
table <- NULL
annotCol <- NULL

for (i in 1:length(samps)){
  for (j in 1:length(status)){
    MeanList <- NULL
    annotCol <- rbind(annotCol, c(samps[i], status[j], "mean"))
    for (k in 1:length(days)){
      subset<-which(grepl(samps[i],exp.penh$Sample_Name) & grepl(status[j], exp.penh$Sample_Name) & exp.penh$Days==days[k])
      ave <- mean(exp.penh$Value[subset], na.rm = TRUE)
      MeanList <- c(MeanList, ave)
    }
    table <- rbind(table, MeanList)
  }
  mock <- nrow(table) - 2
  sars <- nrow(table) - 1
  flu <- nrow(table)
  table <- rbind(table, diff(rbind(table[mock,],table[sars,])))
  table <- rbind(table, table[sars,]/table[mock,])
  table <- rbind(table, diff(rbind(table[mock,],table[flu,])))
  table <- rbind(table, table[flu,]/table[mock,])
  annotCol <- rbind(annotCol, c(samps[i], status[2], "Diff"))
  annotCol <- rbind(annotCol, c(samps[i], status[2], "Fold"))
  annotCol <- rbind(annotCol, c(samps[i], status[3], "Diff"))
  annotCol <- rbind(annotCol, c(samps[i], status[3], "Fold"))
  
}

CombTable <- cbind(annotCol, table)
colnames(CombTable) <- c("Sample", "Status", "Metric", paste("day", days))

file <- "C:/Users/Hunter/Desktop/ohsu/Thesis/tables/May_2014i_PenhSummaryTable.csv"
write.csv(CombTable, file, row.names = FALSE)



############################3



file.name <- "C:/Users/Hunter/Desktop/ohsu/Thesis/data/April_2014i buxco (1).txt"

chunk.size <- 500
db.name <- file.path("C:/Users/Hunter/Desktop/ohsu/Thesis/data", "bux_test4.db")
parse.buxco(file.name=file.name, chunk.size=chunk.size, db.name=db.name, verbose=FALSE)
bux.db <- makeBuxcoDB(db.name=db.name)
addAnnotation(bux.db, query=day.infer.query, index=FALSE)
addAnnotation(bux.db, query=break.type.query, index=TRUE)


namesTable <- matrix(unlist(strsplit(samples(bux.db), " ")), ncol=3, byrow=TRUE)
samps <- unique(namesTable[,1])

transl.table <- data.frame ( sample_status = c( "mock" , "sars" , "flu" ),
  Inf_Status = c( "Mock" , "SARS" , "Flu" ), stringsAsFactors = F)

use.dta <- data.frame( samples = samples ( bux.db ), sample_status = 
  sapply ( strsplit ( samples ( bux.db ) , " "), "[", 3) , stringsAsFactors = F )


transl.table2 <- data.frame ( sample_name = samps,
  Sample = samps, stringsAsFactors = F)

use.dta2 <- data.frame( samples = samples ( bux.db ), sample_strain = 
  sapply ( strsplit ( samples ( bux.db ) , " "), "[", 1) , stringsAsFactors = F )


merge.dta <- merge ( use.dta , transl.table , by = "sample_status" )
merge.dta <- merge.dta [ , -which ( names ( merge.dta ) == "sample_status" ) ]



add.labels.by.sample ( bux.db , merge.dta )
add.labels.by.sample ( bux.db , use.dta2 )
annoLevels ( bux.db )


  mvtsplot ( bux.db , plot.value = "f" , Break_type_label = "EXP" , outer.group.name = "Inf_Status" ,
    sample_strain = samps[1], outer.cols = c( Flu = "black" , SARS = "brown" , Mock = "blue" ))


fileAdd <- "C:/Users/Hunter/Desktop/ohsu/Thesis/Figures/BuxcoBoxplots/"

for (i in 1:length(samps)){
  mvtsplot ( bux.db , plot.value = "f" , Break_type_label = "EXP" , outer.group.name = "Inf_Status" ,
    sample_strain = samps[i], outer.cols = c( Flu = "black" , SARS = "brown" , Mock = "blue" ))
  dev.print(pdf, paste(fileAdd, samps[i], "_FreqMVTSplot", ".pdf", sep=""))
  dev.off()

  mvtsplot ( bux.db , plot.value = "Penh" , Break_type_label = "EXP" , outer.group.name = "Inf_Status" ,
    sample_strain = samps[i], outer.cols = c( Flu = "black" , SARS = "brown" , Mock = "blue" ))
  dev.print(pdf, paste(fileAdd, samps[i], "_PenhMVTSplot", ".pdf", sep=""))
  dev.off()
}


exp.f <- retrieveData(bux.db, variables="f", Break_type_label ='EXP')
exp.penh <- retrieveData(bux.db, variables="Penh", Break_type_label ='EXP')

days <- unique(exp.f$Days)
status <- transl.table[,1]
table <- NULL
annotCol <- NULL

for (i in 1:length(samps)){
  for (j in 1:length(status)){
    MeanList <- NULL
    annotCol <- rbind(annotCol, c(samps[i], status[j], "mean"))
    for (k in 1:length(days)){
      subset<-which(grepl(samps[i],exp.f$Sample_Name) & grepl(status[j], exp.f$Sample_Name) & exp.f$Days==days[k])
      ave <- mean(exp.f$Value[subset], na.rm = TRUE)
      MeanList <- c(MeanList, ave)
    }
    table <- rbind(table, MeanList)
  }
  mock <- nrow(table) - 2
  sars <- nrow(table) - 1
  flu <- nrow(table)
  table <- rbind(table, diff(rbind(table[mock,],table[sars,])))
  table <- rbind(table, table[sars,]/table[mock,])
  table <- rbind(table, diff(rbind(table[mock,],table[flu,])))
  table <- rbind(table, table[flu,]/table[mock,])
  annotCol <- rbind(annotCol, c(samps[i], status[2], "Diff"))
  annotCol <- rbind(annotCol, c(samps[i], status[2], "Fold"))
  annotCol <- rbind(annotCol, c(samps[i], status[3], "Diff"))
  annotCol <- rbind(annotCol, c(samps[i], status[3], "Fold"))
  
}

CombTable <- cbind(annotCol, table)
colnames(CombTable) <- c("Sample", "Status", "Metric", paste("day", days))

file <- "C:/Users/Hunter/Desktop/ohsu/Thesis/tables/April_2014i_FreqSummaryTable.csv"
write.csv(CombTable, file, row.names = FALSE)

### Penh

days <- unique(exp.penh$Days)
status <- transl.table[,1]
table <- NULL
annotCol <- NULL

for (i in 1:length(samps)){
  for (j in 1:length(status)){
    MeanList <- NULL
    annotCol <- rbind(annotCol, c(samps[i], status[j], "mean"))
    for (k in 1:length(days)){
      subset<-which(grepl(samps[i],exp.penh$Sample_Name) & grepl(status[j], exp.penh$Sample_Name) & exp.penh$Days==days[k])
      ave <- mean(exp.penh$Value[subset], na.rm = TRUE)
      MeanList <- c(MeanList, ave)
    }
    table <- rbind(table, MeanList)
  }
  mock <- nrow(table) - 2
  sars <- nrow(table) - 1
  flu <- nrow(table)
  table <- rbind(table, diff(rbind(table[mock,],table[sars,])))
  table <- rbind(table, table[sars,]/table[mock,])
  table <- rbind(table, diff(rbind(table[mock,],table[flu,])))
  table <- rbind(table, table[flu,]/table[mock,])
  annotCol <- rbind(annotCol, c(samps[i], status[2], "Diff"))
  annotCol <- rbind(annotCol, c(samps[i], status[2], "Fold"))
  annotCol <- rbind(annotCol, c(samps[i], status[3], "Diff"))
  annotCol <- rbind(annotCol, c(samps[i], status[3], "Fold"))
  
}

CombTable <- cbind(annotCol, table)
colnames(CombTable) <- c("Sample", "Status", "Metric", paste("day", days))

file <- "C:/Users/Hunter/Desktop/ohsu/Thesis/tables/April_2014i_PenhSummaryTable.csv"
write.csv(CombTable, file, row.names = FALSE)


#########################################################



file.name <- "C:/Users/Hunter/Desktop/ohsu/Thesis/data/May_2015i (1).txt"

chunk.size <- 500
db.name <- file.path("C:/Users/Hunter/Desktop/ohsu/Thesis/data", "bux_test5.db")
parse.buxco(file.name=file.name, chunk.size=chunk.size, db.name=db.name, verbose=FALSE)
bux.db <- makeBuxcoDB(db.name=db.name)
addAnnotation(bux.db, query=day.infer.query, index=FALSE)
addAnnotation(bux.db, query=break.type.query, index=TRUE)


namesTable <- matrix(unlist(strsplit(samples(bux.db), " ")), ncol=3, byrow=TRUE)
samps <- unique(namesTable[,1])

transl.table <- data.frame ( sample_status = c( "mock" , "sars" , "flu", "61" ),
  Inf_Statuss = c( "Mock" , "SARS" , "Flu", "Mock" ), stringsAsFactors = F)

use.dta <- data.frame( samples = samples ( bux.db ), sample_status = 
  sapply ( strsplit ( samples ( bux.db ) , " "), "[", 3) , stringsAsFactors = F )


transl.table2 <- data.frame ( sample_name = samps,
  Sample = samps, stringsAsFactors = F)

use.dta2 <- data.frame( samples = samples ( bux.db ), sample_strain = 
  sapply ( strsplit ( samples ( bux.db ) , " "), "[", 1) , stringsAsFactors = F )


merge.dta <- merge ( use.dta , transl.table , by = "sample_status" )
merge.dta <- merge.dta [ , -which ( names ( merge.dta ) == "sample_status" ) ]


add.labels.by.sample ( bux.db , merge.dta )
add.labels.by.sample ( bux.db , use.dta2 )
annoLevels ( bux.db )

transl.table <- data.frame ( sample_status = c( "mock" , "sars" , "flu"),
  Inf_Status = c( "Mock" , "SARS" , "Flu"), stringsAsFactors = F)

samps <- unique(use.dta2[,2])

  mvtsplot ( bux.db , plot.value = "f" , Break_type_label = "EXP" , outer.group.name = "Inf_Statuss" ,
    sample_strain = samps[1], outer.cols = c( Flu = "black" , SARS = "brown" , Mock = "blue" ))


fileAdd <- "C:/Users/Hunter/Desktop/ohsu/Thesis/Figures/BuxcoBoxplots/"

for (i in 1:length(samps)){
  mvtsplot ( bux.db , plot.value = "f" , Break_type_label = "EXP" , outer.group.name = "Inf_Statuss" ,
    sample_strain = samps[i], outer.cols = c( Flu = "black" , SARS = "brown" , Mock = "blue" ))
  dev.print(pdf, paste(fileAdd, samps[i], "_FreqMVTSplot", ".pdf", sep=""))
  dev.off()

  mvtsplot ( bux.db , plot.value = "Penh" , Break_type_label = "EXP" , outer.group.name = "Inf_Statuss" ,
    sample_strain = samps[i], outer.cols = c( Flu = "black" , SARS = "brown" , Mock = "blue" ))
  dev.print(pdf, paste(fileAdd, samps[i], "_PenhMVTSplot", ".pdf", sep=""))
  dev.off()
}


exp.f <- retrieveData(bux.db, variables="f", Break_type_label ='EXP')
exp.penh <- retrieveData(bux.db, variables="Penh", Break_type_label ='EXP')

days <- unique(exp.f$Days)
status <- transl.table[,1]
table <- NULL
annotCol <- NULL

for (i in 1:length(samps)){
  for (j in 1:length(status)){
    MeanList <- NULL
    annotCol <- rbind(annotCol, c(samps[i], status[j], "mean"))
    for (k in 1:length(days)){
      subset<-which(grepl(samps[i],exp.f$Sample_Name) & grepl(status[j], exp.f$Sample_Name) & exp.f$Days==days[k])
      ave <- mean(exp.f$Value[subset], na.rm = TRUE)
      MeanList <- c(MeanList, ave)
    }
    table <- rbind(table, MeanList)
  }
  mock <- nrow(table) - 2
  sars <- nrow(table) - 1
  flu <- nrow(table)
  table <- rbind(table, diff(rbind(table[mock,],table[sars,])))
  table <- rbind(table, table[sars,]/table[mock,])
  table <- rbind(table, diff(rbind(table[mock,],table[flu,])))
  table <- rbind(table, table[flu,]/table[mock,])
  annotCol <- rbind(annotCol, c(samps[i], status[2], "Diff"))
  annotCol <- rbind(annotCol, c(samps[i], status[2], "Fold"))
  annotCol <- rbind(annotCol, c(samps[i], status[3], "Diff"))
  annotCol <- rbind(annotCol, c(samps[i], status[3], "Fold"))
  
}

CombTable <- cbind(annotCol, table)
colnames(CombTable) <- c("Sample", "Status", "Metric", paste("day", days))

file <- "C:/Users/Hunter/Desktop/ohsu/Thesis/tables/May_2015i_FreqSummaryTable.csv"
write.csv(CombTable, file, row.names = FALSE)

### Penh

days <- unique(exp.penh$Days)
status <- transl.table[,1]
table <- NULL
annotCol <- NULL

for (i in 1:length(samps)){
  for (j in 1:length(status)){
    MeanList <- NULL
    annotCol <- rbind(annotCol, c(samps[i], status[j], "mean"))
    for (k in 1:length(days)){
      subset<-which(grepl(samps[i],exp.penh$Sample_Name) & grepl(status[j], exp.penh$Sample_Name) & exp.penh$Days==days[k])
      ave <- mean(exp.penh$Value[subset], na.rm = TRUE)
      MeanList <- c(MeanList, ave)
    }
    table <- rbind(table, MeanList)
  }
  mock <- nrow(table) - 2
  sars <- nrow(table) - 1
  flu <- nrow(table)
  table <- rbind(table, diff(rbind(table[mock,],table[sars,])))
  table <- rbind(table, table[sars,]/table[mock,])
  table <- rbind(table, diff(rbind(table[mock,],table[flu,])))
  table <- rbind(table, table[flu,]/table[mock,])
  annotCol <- rbind(annotCol, c(samps[i], status[2], "Diff"))
  annotCol <- rbind(annotCol, c(samps[i], status[2], "Fold"))
  annotCol <- rbind(annotCol, c(samps[i], status[3], "Diff"))
  annotCol <- rbind(annotCol, c(samps[i], status[3], "Fold"))
  
}

CombTable <- cbind(annotCol, table)
colnames(CombTable) <- c("Sample", "Status", "Metric", paste("day", days))

file <- "C:/Users/Hunter/Desktop/ohsu/Thesis/tables/May_2015i_PenhSummaryTable.csv"
write.csv(CombTable, file, row.names = FALSE)
