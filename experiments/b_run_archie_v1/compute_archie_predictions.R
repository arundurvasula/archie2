library(data.table)

load("full_model.Rdata")

args <- commandArgs(trailingOnly = T)
inputfile <- args[1]
outputfile <-args[2]

test <- read.table(inputfile)

HAP <- test$V209
CHR <- test$V210
WSTART <- test$V211
WEND <- test$V212
PREDICTION <- plogis(predict.glm(model, test[,1:209]))

keep <- data.frame(HAP, CHR, WSTART, WEND, PREDICTION)

write.table(keep, outputfile, quote=T, row.names = F, col.names = T)

