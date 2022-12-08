setwd("/Users/Shared/KKMDocuments/Documents/Karen/Structure/Pseudorca/Epigenetic.aging")

library(dplyr)
library(tidyr)
library(swfscMisc)

full.CRC.data <- read.csv("new.CRC.age.estimates.csv")
load("Rscripts/Pcra.epi.data.for.Eric.Rdata")
meth.samps <- meth.dat[[1]]$CpG.sum$coverage$id

names(full.CRC.data)[c(1,4,5,6,7,21,22,25,27,29,33)] <- 
  c("crc.id",
    "sex",
    "swfsc.labid",
    "biopsy.id",
    "date.biopsy",
    "age.best",
    "age.confidence",
    "age.min",
    "age.max.CRR",
    "ADULT.max.best",
    "pair.id")

# Eric wants:
#  [1] crc.id (first column in sheet)
#  [5] swfsc.labid (so it matches the id in the methylation data, eg., z00#####)
#  [40] age.class (Calf, Juvenile, Subadult, Adult) - this should only be if the age class (at time of biopsy) is certain
#  [4] sex (Male, Female)
#  [6] biopsy.id (to match with the field ID in the minimum tooth age spreadsheet from Kelly)
#  [7] date.biopsy (YY-MM-DD)
#  [21] age.best (CRC best estimate of age at biopsy)
#  [22] age.confidence (confidence rating 1:5)
#  [25] age.min (minimum age at biopsy)
#  [39] age.max (best maximum age at biopsy)
#  [33] pair.id (number grouping same individual biopsied multiple times)
             

full.CRC.data$age.max <- sapply(1:97, function(i){
  min(full.CRC.data$ADULT.max.best[i], full.CRC.data$age.max.CRR[i], na.rm = TRUE)
})

full.CRC.data$swfsc.labid <- as.numeric(full.CRC.data$swfsc.labid)
full.CRC.data$swfsc.labid[-c(31,35)] <- paste0("z0",zero.pad(full.CRC.data$swfsc.labid[-c(31,35)]))

age.class <- full.CRC.data[,c(5,4,21),]
age.class$best <- "adult"
age.class$best[which(age.class$sex=="Female" & age.class$age.best < 10)] <- "sub-adult"
age.class$best[which(age.class$sex=="Female" & age.class$age.best < 6)] <- "juvenile"
age.class$best[which(age.class$sex=="Female" & age.class$age.best < 3)] <- "calf"

age.class$best[which(age.class$sex=="Male")] <- "adult-PM"
age.class$best[which(age.class$sex=="Male" & age.class$age.best < 25)] <- "adult-SM"
age.class$best[which(age.class$sex=="Male" & age.class$age.best < 15)] <- "sub-adult"
age.class$best[which(age.class$sex=="Male" & age.class$age.best < 9)] <- "juvenile"
age.class$best[which(age.class$sex=="Male" & age.class$age.best < 3)] <- "calf"
full.CRC.data$age.class <- age.class$best

eia.datasheet.1 <- full.CRC.data[which(full.CRC.data$swfsc.labid %in% meth.samps),c(1,5,40,4,6,7,21,22,25,39,33)]

write.csv(eia.datasheet.1,"eia.datasheet.1.csv")