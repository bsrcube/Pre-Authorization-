#############################################
##          Pre - Authorization            ##
#############################################

## Clear the workspace
rm(list = ls(all.names = TRUE))

## Working Directory
setwd("F:/INSOFE/00Project/Data/")
# setwd(choose.dir())

## Some paths and file names
iPath = "F:/INSOFE/00Project/Data/"
iFile = "PriorAuth_Data.csv"

## Libraries
library(lubridate)
library(caret)
library(reshape2)
library(dplyr)
library(ggplot2)
library(forcats)
library(stringr)
library(scales)


## User defined function for getting the various metrics from the Confusion Matrix
getMetrics <- function (predVals, actuVals){
  
  confusionMatrix = table(actuVals, predVals)
  # print("Confusion Matrix - ", confusionMatrix)
  
  accuracy = sum(diag(confusionMatrix)) / sum(confusionMatrix)
  cat("Accuracy -", accuracy, "\n")
  
  recall = confusionMatrix[2,2] / sum(confusionMatrix[2,])
  cat("Recall -", recall, "\n")
  
  precision = confusionMatrix[2,2] / sum(confusionMatrix[,2])
  cat("Precision -", precision, "\n")
  
  specificity = confusionMatrix[1,1] / sum(confusionMatrix[1,])
  cat("Specificity -", specificity, "\n")
  
}


## Function to get the new data frame with reduced levels
dfWithReducedLevels <- function (df) {
  
  colnames(df)[grep(pattern = ".*(_CID)$", x = colnames(df))] = str_sub(string = colnames(df)[grep(pattern = ".*(_CID)$", x = colnames(df))], start = 0, end = nchar(x = colnames(df)[grep(pattern = ".*(_CID)$", x = colnames(df))]) - 4)
  
  df = data.frame(apply(X = df, MARGIN = 2, FUN = as.character), stringsAsFactors = FALSE)
  
  df[,"Drug"] = as.factor(str_c("Drug_",df$Drug))
  df[,"DrugSubClass"] = as.factor(str_c("DrugSubClass_",df$DrugSubClass))
  df[,"Drug_Chemical_Name"] = as.factor(str_c("Drug_Chemical_Name_",df$Drug_Chemical_Name))
  df[,"DrugGroup"] = as.factor(str_c("DrugGroup_",df$DrugGroup))
  df[,"DrugClass"] = as.factor(str_c("DrugClass_",df$DrugClass))
  # df[,"RxGroupId"] = as.factor(str_c("RxGroupID_",df$RxGroupId))
  # df[,"GPI"] = as.factor(str_c("GPI_",df$GPI))
  df[,"PCN"] = as.factor(str_c("PCN_",df$PCN))
  df[,"Bin"] = as.factor(str_c("Bin_",df$Bin))
  # df[,"NDC"] = as.factor(str_c("NDC_",df$NDC))
  
  return(df)
  
}

##########################################################################################################


## Read the input data
paData = read.csv(file = paste(iPath, iFile, sep = ""), header = TRUE, check.names = FALSE)
paData.org = paData
# > dim(unique(x = paData))
# [1] 6738   15


## Structure study and summary statistics
str(paData)
summary(paData)


## Setting the variables to their respective types
## Converting the Target variable to numeric and factor
paData[,"Target"] = as.factor(ifelse(paData[,"Target"] == TRUE, 1, 0))

## Plot the target distribution here 
# paData = unqData
dist.target = data.frame(Target = levels(paData$Target), Freq = c(table(paData$Target)), Percent = str_c(c(table(paData$Target)[1],table(paData$Target)[2]), "  ( ~ ", c(round(table(paData$Target)[1] / (table(paData$Target)[1] + table(paData$Target)[2]), 3) * 100, round(table(paData$Target)[2] / (table(paData$Target)[1] + table(paData$Target)[2]), 3) * 100), " % )"))

ggplot(dist.target, aes(x = Target, y = Freq, fill = Target)) + geom_bar(stat = "identity") +
    ggtitle("Data Distribution ~ Target Variable") + geom_text(aes(label = Percent), vjust = 1.5)
  
# ggsave(filename = "RawDataTargetDist.png")
  
rm(dist.target)

  
## Fix the Transdate column to have the same delimiters
paData[,"TransDate"] = mdy(paData[,"TransDate"])


## Handling duplicates and probable wrong entries
unqData = unique(paData)
paData = unqData

## Run chi-sq test to check the dependence of the predictors on the Target attribute
# chisq.test(paData$Drug, paData$Target) # 0.0223
# chisq.test(paData$DrugClass, paData$Target) # 9.999e-05
# chisq.test(paData$DrugSubClass, paData$Target) # 9.999e-05
# chisq.test(paData$DrugGroup, paData$Target) # 9.999e-05
# chisq.test(paData$Drug_Chemical_Name, paData$Target) # 9.999e-05
# chisq.test(paData$Bin, paData$Target) # 9.999e-05
# chisq.test(paData$PCN, paData$Target) # 9.999e-05
# chisq.test(paData$NDC, paData$Target) # 1
# chisq.test(paData$State, paData$Target) # 9.999e-05
# chisq.test(paData$TransDate, paData$Target) # 9.999e-05
# chisq.test(paData$RxGroupId, paData$Target) # 0.09969
# chisq.test(paData$GPI, paData$Target) # 0.8417
# chisq.test(paData$UserID, paData$Target) # 1
# chisq.test(paData$DoctorID, paData$Target) # 1



# chisq.test(paData$Drug, paData$Target, simulate.p.value = TRUE, B = 10000) # 0.0223
# chisq.test(paData$DrugClass, paData$Target, simulate.p.value = TRUE, B = 10000) # 9.999e-05
# chisq.test(paData$DrugSubClass, paData$Target, simulate.p.value = TRUE, B = 10000) # 9.999e-05
# chisq.test(paData$DrugGroup, paData$Target, simulate.p.value = TRUE, B = 10000) # 9.999e-05
# chisq.test(paData$Drug_Chemical_Name, paData$Target, simulate.p.value = TRUE, B = 10000) # 9.999e-05
# chisq.test(paData$Bin, paData$Target, simulate.p.value = TRUE, B = 10000) # 9.999e-05
# chisq.test(paData$PCN, paData$Target, simulate.p.value = TRUE, B = 10000) # 9.999e-05
# chisq.test(paData$NDC, paData$Target, simulate.p.value = TRUE, B = 10000) # 1
# chisq.test(paData$State, paData$Target, simulate.p.value = TRUE, B = 10000) # 9.999e-05
# chisq.test(paData$TransDate, paData$Target, simulate.p.value = TRUE, B = 10000) # 9.999e-05
# chisq.test(paData$RxGroupId, paData$Target, simulate.p.value = TRUE, B = 10000) # 0.09969
# chisq.test(paData$GPI, paData$Target, simulate.p.value = TRUE, B = 10000) # 0.8417
# chisq.test(paData$UserID, paData$Target, simulate.p.value = TRUE, B = 10000) # 1
# chisq.test(paData$DoctorID, paData$Target, simulate.p.value = TRUE, B = 10000) # 1


## Dropping columns that are not useful
paData = paData[,!colnames(paData) %in% c("UserID", "DoctorID", "NDC", "GPI", "RxGroupId")]
# > dim(unique(x = paData0))
# [1] 6649   13

str(paData)


dateBreaks = seq(from = ymd(levels(as.factor(paData$TransDate))[1]), 
                   to = ymd(levels(as.factor(paData$TransDate))[length(levels(as.factor(paData$TransDate)))]),
                   by = "1 day")
  
dist.TransDate = data.frame(Dates = levels(as.factor(paData$TransDate)), Freq = c(table(paData$TransDate)), row.names = NULL)
  
colors = c("#FF6666", "#009999")
ggplot(dist.TransDate, aes(x = as.Date(Dates), y = Freq, fill = Dates)) + 
    geom_bar(stat = "identity") + 
    ggtitle("Data Distribution ~ TransDate  <=> Train (Red) & Test (Green) data split") + 
    geom_text(aes(label = Freq), vjust = 1.2, size = 3.2) + 
    scale_fill_manual(values = c(rep(colors[1],15), rep(colors[2], 2))) +
    scale_x_date(breaks = dateBreaks, labels = date_format(format = "%b %d")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 8)) +
    xlab("Dates")
  
# ggsave("TrainTestSplit.png")

rm(dist.TransDate, dateBreaks)



## Time based sampling
set.seed(2357)
paData = paData[order(paData$TransDate), ]
trainData = paData[paData$TransDate < levels(as.factor(paData$TransDate))[length(levels(as.factor(paData$TransDate))) - 1], ]
testData = paData[paData$TransDate >= levels(as.factor(paData$TransDate))[length(levels(as.factor(paData$TransDate))) - 1], ]

## Set the Transdate attribute to factors
trainData$TransDate = as.factor(as.character(trainData$TransDate))
testData$TransDate = as.factor(as.character(testData$TransDate))

paData0 = trainData


## Check the data distribution wrt the Target
prop.table(table(testData$Target))
prop.table(table(trainData$Target))


## Plotting Elbow curve to get the k-values for the clusters

## K-means clustering on Drug (K = 4/5)
Temp = data.frame(paData0 %>% group_by(Drug) %>% select(Target, Drug))
Temp_gr = split(Temp, Temp$Target)
FDrug = data.frame(Temp_gr$`0`)
FDrug = FDrug %>% group_by(Drug) %>% summarize("FCount" = n())
TDrug = data.frame(Temp_gr$`1`)
TDrug = TDrug %>% group_by(Drug) %>% summarize("TCount" = n())
Drug_m = full_join(x = FDrug, y = TDrug, by = "Drug")
Drug_m[is.na(Drug_m)] = 0

dDrug = data.frame(Drug_m %>% select(everything()) %>% group_by(Drug) %>% summarize("FPercent" = FCount/sum(FCount + TCount), "TPercent" = TCount/sum(FCount + TCount))) 

wss_df <- data.frame(x = 1:15)
for (j in 1:10) {
  set.seed(2358 * j)
  tot.wss <- 0
  for (i in 1:15) {
    tot.wss[i] <-
      kmeans(x = dDrug$FPercent,
             centers = i,
             nstart = 10)$tot.withinss
  }
  wss_df[,paste0("Seed_",2358*j)] = tot.wss
}

wss_df = data.frame(wss_df, count = 1:nrow(wss_df))
wss_df = wss_df[,!colnames(wss_df) %in% c("x")]
wss_df.plot = melt(data = wss_df, id = "count")
ggplot(data = wss_df.plot, aes(x = count, y = value, colour = variable, group = variable)) + geom_line() + geom_point(size = 1.8, shape = 16) + xlab("k-values") + ylab("Within Sum of Square Errors") + ggtitle("Drug")
## K-means clustering on Drug (K = 4/5)

rm(i, tot.wss, wss_df, wss_df.plot, Temp, TDrug, FDrug, Drug_m)


## K-means clustering on PCN (k = 4)
Temp = data.frame(paData0 %>% group_by(PCN) %>% select(Target, PCN))
Temp_gr = split(Temp, Temp$Target)
FPCN = data.frame(Temp_gr$`0`)
FPCN = FPCN %>% group_by(PCN) %>% summarize("FCount" = n())
TPCN = data.frame(Temp_gr$`1`)
TPCN = TPCN %>% group_by(PCN) %>% summarize("TCount" = n())
PCN_m = full_join(x = FPCN, y = TPCN, by = "PCN")
PCN_m[is.na(PCN_m)] = 0

dPCN = data.frame(PCN_m %>% group_by(PCN) %>% summarize("FPercent" = FCount/sum(FCount + TCount), "TPercent" = TCount/sum(FCount + TCount))) 

wss_df <- data.frame(x = 1:15)
for (j in 1:10) {
  set.seed(2358 * j)
  tot.wss <- 0
  for (i in 1:15) {
    tot.wss[i] <-
      kmeans(x = dPCN$FPercent,
             centers = i,
             nstart = 10)$tot.withinss
  }
  wss_df[,paste0("Seed_",2358*j)] = tot.wss
}

wss_df = data.frame(wss_df, count = 1:nrow(wss_df))
wss_df = wss_df[,!colnames(wss_df) %in% c("x")]
wss_df.plot = melt(data = wss_df, id = "count")
ggplot(data = wss_df.plot, aes(x = count, y = value, colour = variable, group = variable)) + geom_line() + geom_point(size = 1.8, shape = 16) + xlab("k-values") + ylab("Within Sum of Square Errors") + ggtitle("PCN")
## K-means clustering on PCN (k = 4)

rm(i, tot.wss, wss_df, wss_df.plot, Temp, Temp_gr, TPCN, FPCN, PCN_m)



## K-means clustering on Bin (k = 4)
Temp = data.frame(paData0 %>% group_by(Bin) %>% select(Target, Bin))
Temp_gr = split(Temp, Temp$Target)
FBin = data.frame(Temp_gr$`0`)
FBin = FBin %>% group_by(Bin) %>% summarize("FCount" = n())
TBin = data.frame(Temp_gr$`1`)
TBin = TBin %>% group_by(Bin) %>% summarize("TCount" = n())
Bin_m = full_join(x = FBin, y = TBin, by = "Bin")
Bin_m[is.na(Bin_m)] = 0

dBin = data.frame(Bin_m %>% group_by(Bin) %>% summarize("FPercent" = FCount/sum(FCount + TCount), "TPercent" = TCount/sum(FCount + TCount))) 

wss_df <- data.frame(x = 1:15)
for (j in 1:10) {
  set.seed(2358 * j)
  tot.wss <- 0
  for (i in 1:15) {
    tot.wss[i] <-
      kmeans(x = dBin$FPercent,
             centers = i,
             nstart = 10)$tot.withinss
  }
  wss_df[,paste0("Seed_",2358*j)] = tot.wss
}

wss_df = data.frame(wss_df, count = 1:nrow(wss_df))
wss_df = wss_df[,!colnames(wss_df) %in% c("x")]
wss_df.plot = melt(data = wss_df, id = "count")
ggplot(data = wss_df.plot, aes(x = count, y = value, colour = variable, group = variable)) + geom_line() + geom_point(size = 1.8, shape = 16) + xlab("k-values") + ylab("Within Sum of Square Errors") + ggtitle("Bin")
## K-means clustering on Bin (k = 4)

rm(i, tot.wss, wss_df, wss_df.plot, Temp, Temp_gr, TBin, FBin, Bin_m)


## K-means clustering on DrugClass (k = 5)
Temp = data.frame(paData0 %>% group_by(DrugClass) %>% select(Target, DrugClass))
Temp_gr = split(Temp, Temp$Target)
FDrugClass = data.frame(Temp_gr$`0`)
FDrugClass = FDrugClass %>% group_by(DrugClass) %>% summarize("FCount" = n())
TDrugClass = data.frame(Temp_gr$`1`)
TDrugClass = TDrugClass %>% group_by(DrugClass) %>% summarize("TCount" = n())
DrugClass_m = full_join(x = FDrugClass, y = TDrugClass, by = "DrugClass")
DrugClass_m[is.na(DrugClass_m)] = 0

dDrugClass = data.frame(DrugClass_m %>% group_by(DrugClass) %>% summarize("FPercent" = FCount/sum(FCount + TCount), "TPercent" = TCount/sum(FCount + TCount))) 

wss_df <- data.frame(x = 1:15)
for (j in 1:10) {
  set.seed(2358 * j)
  tot.wss <- 0
  for (i in 1:15) {
    tot.wss[i] <-
      kmeans(x = dDrugClass$FPercent,
             centers = i,
             nstart = 10)$tot.withinss
  }
  wss_df[,paste0("Seed_",2358*j)] = tot.wss
}

wss_df = data.frame(wss_df, count = 1:nrow(wss_df))
wss_df = wss_df[,!colnames(wss_df) %in% c("x")]
wss_df.plot = melt(data = wss_df, id = "count")
ggplot(data = wss_df.plot, aes(x = count, y = value, colour = variable, group = variable)) + geom_line() + geom_point(size = 1.8, shape = 16) + xlab("k-values") + ylab("Within Sum of Square Errors") + ggtitle("DrugClass")
## K-means clustering on DrugClass (k = 5)

rm(i, tot.wss, wss_df, wss_df.plot, Temp, Temp_gr, TDrugClass, FDrugClass, DrugClass_m)


## K-means clustering on DrugSubClass (k = 4/5)
Temp = data.frame(paData0 %>% group_by(DrugSubClass) %>% select(Target, DrugSubClass))
Temp_gr = split(Temp, Temp$Target)
FDrugSubClass = data.frame(Temp_gr$`0`)
FDrugSubClass = FDrugSubClass %>% group_by(DrugSubClass) %>% summarize("FCount" = n())
TDrugSubClass = data.frame(Temp_gr$`1`)
TDrugSubClass = TDrugSubClass %>% group_by(DrugSubClass) %>% summarize("TCount" = n())
DrugSubClass_m = full_join(x = FDrugSubClass, y = TDrugSubClass, by = "DrugSubClass")
DrugSubClass_m[is.na(DrugSubClass_m)] = 0

dDrugSubClass = data.frame(DrugSubClass_m %>% group_by(DrugSubClass) %>% summarize("FPercent" = FCount/sum(FCount + TCount), "TPercent" = TCount/sum(FCount + TCount))) 

wss_df <- data.frame(x = 1:15)
for (j in 1:10) {
  set.seed(2358 * j)
  tot.wss <- 0
  for (i in 1:15) {
    tot.wss[i] <-
      kmeans(x = dDrugSubClass$FPercent,
             centers = i,
             nstart = 10)$tot.withinss
  }
  wss_df[,paste0("Seed_",2358*j)] = tot.wss
}

wss_df = data.frame(wss_df, count = 1:nrow(wss_df))
wss_df = wss_df[,!colnames(wss_df) %in% c("x")]
wss_df.plot = melt(data = wss_df, id = "count")
ggplot(data = wss_df.plot, aes(x = count, y = value, colour = variable, group = variable)) + geom_line() + geom_point(size = 1.8, shape = 16) + xlab("k-values") + ylab("Within Sum of Square Errors") + ggtitle("DrugSubClass")
## K-means clustering on DrugSubClass (k = 4/5)

rm(i, tot.wss, wss_df, wss_df.plot, Temp, Temp_gr, TDrugSubClass, FDrugSubClass, DrugSubClass_m)


## K-means clustering on Drug_Chemical_Name (k = 4)
Temp = data.frame(paData0 %>% group_by(Drug_Chemical_Name) %>% select(Target, Drug_Chemical_Name))
Temp_gr = split(Temp, Temp$Target)
FDrug_Chemical_Name = data.frame(Temp_gr$`0`)
FDrug_Chemical_Name = FDrug_Chemical_Name %>% group_by(Drug_Chemical_Name) %>% summarize("FCount" = n())
TDrug_Chemical_Name = data.frame(Temp_gr$`1`)
TDrug_Chemical_Name = TDrug_Chemical_Name %>% group_by(Drug_Chemical_Name) %>% summarize("TCount" = n())
Drug_Chemical_Name_m = full_join(x = FDrug_Chemical_Name, y = TDrug_Chemical_Name, by = "Drug_Chemical_Name")
Drug_Chemical_Name_m[is.na(Drug_Chemical_Name_m)] = 0

dDrug_Chemical_Name = data.frame(Drug_Chemical_Name_m %>% group_by(Drug_Chemical_Name) %>% summarize("FPercent" = FCount/sum(FCount + TCount), "TPercent" = TCount/sum(FCount + TCount))) 

wss_df <- data.frame(x = 1:15)
for (j in 1:10) {
  set.seed(2358 * j)
  tot.wss <- 0
  for (i in 1:15) {
    tot.wss[i] <-
      kmeans(x = dDrug_Chemical_Name$FPercent,
             centers = i,
             nstart = 10)$tot.withinss
  }
  wss_df[,paste0("Seed_",2358*j)] = tot.wss
}
wss_df = data.frame(wss_df, count = 1:nrow(wss_df))
wss_df = wss_df[,!colnames(wss_df) %in% c("x")]
wss_df.plot = melt(data = wss_df, id = "count")
ggplot(data = wss_df.plot, aes(x = count, y = value, colour = variable, group = variable)) + geom_line() + geom_point(size = 1.8, shape = 16) + xlab("k-values") + ylab("Within Sum of Square Errors") + ggtitle("Drug_Chemical_Name")
## K-means clustering on Drug_Chemical_Name (k = 4)

rm(i, tot.wss, wss_df, wss_df.plot, Temp, Temp_gr, TDrug_Chemical_Name, FDrug_Chemical_Name, Drug_Chemical_Name_m)


## K-means clustering on DrugGroup (k = 5/6)
Temp = data.frame(paData0 %>% group_by(DrugGroup) %>% select(Target, DrugGroup))
Temp_gr = split(Temp, Temp$Target)
FDrugGroup = data.frame(Temp_gr$`0`)
FDrugGroup = FDrugGroup %>% group_by(DrugGroup) %>% summarize("FCount" = n())
TDrugGroup = data.frame(Temp_gr$`1`)
TDrugGroup = TDrugGroup %>% group_by(DrugGroup) %>% summarize("TCount" = n())
DrugGroup_m = full_join(x = FDrugGroup, y = TDrugGroup, by = "DrugGroup")
DrugGroup_m[is.na(DrugGroup_m)] = 0

dDrugGroup = data.frame(DrugGroup_m %>% group_by(DrugGroup) %>% summarize("FPercent" = FCount/sum(FCount + TCount), "TPercent" = TCount/sum(FCount + TCount))) 

wss_df <- data.frame(x = 1:15)
for (j in 1:10) {
  set.seed(2358 * j)
  tot.wss <- 0
  for (i in 1:15) {
    tot.wss[i] <-
      kmeans(x = dDrugGroup$FPercent,
             centers = i,
             nstart = 10)$tot.withinss
  }
  wss_df[,paste0("Seed_",2358*j)] = tot.wss
}

wss_df = data.frame(wss_df, count = 1:nrow(wss_df))
wss_df = wss_df[,!colnames(wss_df) %in% c("x")]
wss_df.plot = melt(data = wss_df, id = "count")
ggplot(data = wss_df.plot, aes(x = count, y = value, colour = variable, group = variable)) + geom_line() + geom_point(size = 1.8, shape = 16) + xlab("k-values") + ylab("Within Sum of Square Errors") + ggtitle("DrugGroup")
## K-means clustering on DrugGroup (k = 5/6)

rm(i, tot.wss, wss_df, wss_df.plot, Temp, Temp_gr, TDrugGroup, FDrugGroup, DrugGroup_m)


## Running k-means on all the attributes after getting the k-values from elbow curve method
# set.seed(1985)
kmDrug = kmeans(x = dDrug[,!colnames(dDrug) %in% c("Drug")], centers = 5, nstart = 10)
dDrug[,"Drug_CID"] = as.factor(kmDrug$cluster)
dDrug = dDrug[,!colnames(dDrug) %in% c("FPercent", "TPercent")]
traindf = inner_join(x = paData0, y = dDrug, "Drug")
testdf = inner_join(x = testData, y = dDrug, "Drug") 
rm(kmDrug)

kmDrugSubClass = kmeans(x = dDrugSubClass[,!colnames(dDrugSubClass) %in% c("DrugSubClass")], centers = 5, nstart = 10)
dDrugSubClass[,"DrugSubClass_CID"] = as.factor(kmDrugSubClass$cluster)
dDrugSubClass = dDrugSubClass[,!colnames(dDrugSubClass) %in% c("FPercent", "TPercent")]
traindf = inner_join(x = traindf, y = dDrugSubClass, "DrugSubClass")
testdf = inner_join(x = testdf, y = dDrugSubClass, "DrugSubClass")
rm(kmDrugSubClass)

kmDrugClass = kmeans(x = dDrugClass[,!colnames(dDrugClass) %in% c("DrugClass")], centers = 5, nstart = 10)
dDrugClass[,"DrugClass_CID"] = as.factor(kmDrugClass$cluster)
dDrugClass = dDrugClass[,!colnames(dDrugClass) %in% c("FPercent", "TPercent")]
traindf = inner_join(x = traindf, y = dDrugClass, "DrugClass")
testdf = inner_join(x = testdf, y = dDrugClass, "DrugClass")
rm(kmDrugClass)

kmDrug_CN = kmeans(x = dDrug_Chemical_Name[,!colnames(dDrug_Chemical_Name) %in% c("Drug_Chemical_Name")], centers = 4, nstart = 10)
dDrug_Chemical_Name[,"Drug_Chemical_Name_CID"] = as.factor(kmDrug_CN$cluster)
dDrug_Chemical_Name = dDrug_Chemical_Name[,!colnames(dDrug_Chemical_Name) %in% c("FPercent", "TPercent")]
traindf = inner_join(x = traindf, y = dDrug_Chemical_Name, "Drug_Chemical_Name")
testdf = inner_join(x = testdf, y = dDrug_Chemical_Name, "Drug_Chemical_Name")
rm(kmDrug_CN)

kmDrugGroup = kmeans(x = dDrugGroup[,!colnames(dDrugGroup) %in% c("DrugGroup")], centers = 6, nstart = 10)
dDrugGroup[,"DrugGroup_CID"] = as.factor(kmDrugGroup$cluster)
dDrugGroup = dDrugGroup[,!colnames(dDrugGroup) %in% c("FPercent", "TPercent")]
traindf = inner_join(x = traindf, y = dDrugGroup, "DrugGroup")
testdf = inner_join(x = testdf, y = dDrugGroup, "DrugGroup")
rm(kmDrugGroup)

kmBin = kmeans(x = dBin[,!colnames(dBin) %in% c("Bin")], centers = 4, nstart = 10)
dBin[,"Bin_CID"] = as.factor(kmBin$cluster)
dBin = dBin[,!colnames(dBin) %in% c("FPercent", "TPercent")]
traindf = inner_join(x = traindf, y = dBin, "Bin")
testdf = inner_join(x = testdf, y = dBin, "Bin")
rm(kmBin)

kmPCN = kmeans(x = dPCN[,!colnames(dPCN) %in% c("PCN")], centers = 4, nstart = 10)
dPCN[,"PCN_CID"] = as.factor(kmPCN$cluster)
dPCN = dPCN[,!colnames(dPCN) %in% c("FPercent", "TPercent")]
traindf = inner_join(x = traindf, y = dPCN, "PCN")
testdf = inner_join(x = testdf, y = dPCN, "PCN")
rm(kmPCN)

rm(dDrug, dDrugGroup, dDrugClass, dDrug_Chemical_Name, dDrugSubClass, dPCN, dBin)



## Apply k-means and get the reduced levels for each of the attribute 
reqCols = grep(pattern = ".*(_CID)$", x = colnames(traindf), value = TRUE)

fTrainData = dfWithReducedLevels(traindf[,reqCols])
fTrainData = data.frame("Target" = traindf$Target, "State" = traindf$State, "TransDate" = traindf$TransDate, fTrainData)

fTestData = dfWithReducedLevels(testdf[,reqCols])
fTestData = data.frame("Target" = testdf$Target, "State" = testdf$State, "TransDate" = testdf$TransDate, fTestData)


paData0 = fTrainData
# str(paData0)

## For Benchmark models
# fTrainData = paData0
# fTestData = testData


## Setting the types & levels of the attributes right before feeding the input to the models

levels(fTrainData$Drug) = union(levels(fTrainData$Drug), levels(fTestData$Drug))
levels(fTrainData$DrugClass) = union(levels(fTrainData$DrugClass), levels(fTestData$DrugClass))
levels(fTrainData$DrugSubClass) = union(levels(fTrainData$DrugSubClass), levels(fTestData$DrugSubClass))
levels(fTrainData$DrugGroup) = union(levels(fTrainData$DrugGroup), levels(fTestData$DrugGroup))
levels(fTrainData$Drug_Chemical_Name) = union(levels(fTrainData$Drug_Chemical_Name), levels(fTestData$Drug_Chemical_Name))
levels(fTrainData$Bin) = union(levels(fTrainData$Bin), levels(fTestData$Bin))
levels(fTrainData$PCN) = union(levels(fTrainData$PCN), levels(fTestData$PCN))
levels(fTrainData$State) = union(levels(fTrainData$State), levels(fTestData$State))
levels(fTrainData$TransDate) = union(levels(fTrainData$TransDate), levels(fTestData$TransDate))


levels(fTestData$Drug) = union(levels(fTrainData$Drug), levels(fTestData$Drug))
levels(fTestData$DrugClass) = union(levels(fTrainData$DrugClass), levels(fTestData$DrugClass))
levels(fTestData$DrugSubClass) = union(levels(fTrainData$DrugSubClass), levels(fTestData$DrugSubClass))
levels(fTestData$DrugGroup) = union(levels(fTrainData$DrugGroup), levels(fTestData$DrugGroup))
levels(fTestData$Drug_Chemical_Name) = union(levels(fTrainData$Drug_Chemical_Name), levels(fTestData$Drug_Chemical_Name))
levels(fTestData$Bin) = union(levels(fTrainData$Bin), levels(fTestData$Bin))
levels(fTestData$PCN) = union(levels(fTrainData$PCN), levels(fTestData$PCN))
levels(fTestData$State) = union(levels(fTrainData$State), levels(fTestData$State))
levels(fTestData$TransDate) = union(levels(fTrainData$TransDate), levels(fTestData$TransDate))


trainData0 = fTrainData
testData0 = fTestData



## Naive Bayes classifier
library(e1071)
# set.seed(2358)
myModel = naiveBayes(Target ~ . , data = trainData0)
predTrain0 = predict(object = myModel, trainData0)
table(predTrain0, trainData0$Target)
predTest0 = predict(object = myModel, testData0[,!colnames(testData0) %in% c("Target")])
table(predTest0, testData0$Target)
resTrain0 = getMetrics(predTrain0, trainData0$Target)
resTest0 = getMetrics(predTest0, testData0$Target)


library(rpart)
# set.seed(2358)
myModel = rpart(formula = Target ~ . , data = trainData0, method = "class")
predTrain0 = predict(object = myModel, trainData0, type = "class")
table(predTrain0, trainData0$Target)
predTest0 = predict(object = myModel, testData0[,!colnames(testData0) %in% c("Target")], type = "class")
table(predTest0, testData0$Target)
resTrain0 = getMetrics(predTrain0, trainData0$Target)
resTest0 = getMetrics(predTest0, testData0$Target)


library(C50)
# set.seed(2358)
myModel = C5.0(formula = Target ~ . , data = trainData0, method = "class")
predTrain0 = predict(object = myModel, trainData0, type = "class")
table(predTrain0, trainData0$Target)
predTest0 = predict(object = myModel, testData0[,!colnames(testData0) %in% c("Target")], type = "class")
table(predTest0, testData0$Target)
resTrain0 = getMetrics(predTrain0, trainData0$Target)
resTest0 = getMetrics(predTest0, testData0$Target)


library(randomForest)
# set.seed(2358)
myModel = randomForest(Target ~ ., data = trainData0, ntree = 15, replace = TRUE)
predTrain0 = predict(object = myModel, trainData0, type = "class")
table(predTrain0, trainData0$Target)
predTest0 = predict(object = myModel, testData0[,!colnames(testData0) %in% c("Target")], type = "class")
table(predTest0, testData0$Target)
resTrain0 = getMetrics(predTrain0, trainData0$Target)
resTest0 = getMetrics(predTest0, testData0$Target)

