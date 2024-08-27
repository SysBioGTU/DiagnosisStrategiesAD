
# For ML based classification of AsymAD individuals
# Call related packages from library
library(tidyverse)
library(readr)
library(tibble)
library(doParallel)
library(caret)
library(randomForest)
library(openxlsx)
library(plyr)
library(stringr)
library(biomaRt)
setwd("C:\\path of file")
getwd()

set.seed(7/5/2023)
ROSMAP_ensmbl_unique <- openxlsx::read.xlsx("C:/file.xlsx")

gc()

#Filtering protein-coding genes in dataset.We used Biomart package for this purpose.
ensembl <- useMart('ensembl',dataset = 'hsapiens_gene_ensembl')
attributes.data <- c('hgnc_symbol','ensembl_gene_id','gene_biotype')


geneAnnot <- getBM(attributes = attributes.data,
                   filters = 'hgnc_symbol',
                   values = ROSMAP_ensmbl_unique$geneSymbol,
                   mart = ensembl)
proteinCoding <- geneAnnot[geneAnnot$gene_biotype=='protein_coding',]


ROSMAP_ensmbl_unique <- ROSMAP_ensmbl_unique %>% filter(.$geneSymbol %in% proteinCoding$hgnc_symbol)

ROSMAP_ensmbl_unique.1 <- ROSMAP_ensmbl_unique[,-c(2)]


# Read colnames of data which were separated into two files as train and prediction datasets.
trainData <- colnames(read_csv("C:/data.csv"))[-1]
predictionsData <- colnames(read_csv("C:/data.csv"))[-1]

# Data scaling before classification
normData <- function(df) {
  trainData <- matrix(data = sapply(df,as.numeric),ncol = ncol(df))
  colnames(trainData) <- colnames(df); rownames(trainData) <- row.names(df)
  trainData <- t(as.matrix(scale(t(trainData),scale = T, center = T)))
  trainData <- trainData[complete.cases(trainData),]
  return(trainData)
}

#The number of samples comes from consensus non-AD and consensus AD groups. 
#This is for transcriptome data. 
dataLabsRNA <- as.factor(c(rep('Control',178), rep('AD',204)))

# Train data for ROSMAP-transcriptome dataset
gc()
trainDataRNA <- ROSMAP_ensmbl_unique.1 |> 
  dplyr::select(geneSymbol,trainData) |>
  remove_rownames() |>
  column_to_rownames(var = 'geneSymbol')|>
  normData()
dim(trainDataRNA)

trainIndex <- createDataPartition(dataLabsRNA, p = .7, 
                                  list = FALSE, 
                                  times = 1)
sampleLabsTrainRNA <- dataLabsRNA[trainIndex]
sampleLabsTestRNA <- dataLabsRNA[-trainIndex]

trainRNA <- trainDataRNA[,trainIndex]
testRNA <- trainDataRNA[,-trainIndex]

dim(trainRNA)
dim(testRNA)

# Prediction data for ROSMAP-transcriptome dataset
precitonDataRNA <- ROSMAP_ensmbl_unique.1 %>% 
  dplyr::select(geneSymbol,predictionsData) %>% 
  remove_rownames()%>% 
  column_to_rownames(var = 'geneSymbol')|>
  normData()
dim(precitonDataRNA)

# Similar steps with ROSMAP-metabolome dataset
trainDataMet <- readr::read_delim('C:/data.csv') |>
  dplyr::rename(MetID = 1)|>
  dplyr::mutate(MetID = paste0('Met_',MetID)) |>
  remove_rownames() |>
  column_to_rownames(var = 'MetID')|>
  normData()

#The number of samples comes from consensus non-AD and consensus AD groups.
#This is for metabolome dataset.
dataLabsMet<- as.factor(c(rep('Control',53), rep('AD',85)))
trainIndex <- createDataPartition(dataLabsMet, p = .7, 
                                  list = FALSE, 
                                  times = 1)
sampleLabsTrainMet <- dataLabsMet[trainIndex]
sampleLabsTestMet <- dataLabsMet[-trainIndex]

trainMet <- trainDataMet[,trainIndex]
testMet <- trainDataMet[,-trainIndex]

dim(trainMet)
dim(testMet)

# Prediction with Metabolome Data
predictionDataMet <- readr::read_delim('C:/data.csv') %>% 
  dplyr::rename(MetID = 1)|>
  dplyr::mutate(MetID = paste0('Met_',MetID)) |>
  remove_rownames() |>
  column_to_rownames(var = 'MetID')|>
  normData()
dim(predictionDataMet)

# Save train,test and prediction (both data) in a list format
trainingData <- list(RNATrain = trainRNA,
                     MetTrain = trainMet)
testData <- list(RNATest = testRNA,
                 MetTest = testMet)
predictionData <- list(RNAPrediction = precitonDataRNA,
                       MetPrediction = predictionDataMet)

trainingLabs <- list(RNALabs = sampleLabsTrainRNA,
                     MetLabs = sampleLabsTrainMet)
testLabs <- list(RNALabs = sampleLabsTestRNA,
                 MetLabs = sampleLabsTestMet)

gc()

# Recursive Feature Elimination (RFE) method for feature selection
library(doParallel)
library(caret)

sizes <- c(50,100,300,500,1000,1500)
cores <- detectCores()
cl <- makePSOCKcluster(cores-3, outfile="")
registerDoParallel(cl)
resultsOverall <- lapply(seq_along(trainingData), function(y){
  lapply(sizes, function(x){
    seeds <- vector(mode = 'list', length = 51)
    for(i in 1:51) seeds[[i]]<- sample.int(n=1000, size = (x+1), replace = T)
    contr <- caret::rfeControl(functions=rfFuncs, 
                               method="cv", 
                               number=5, 
                               repeats=5, 
                               saveDetails = TRUE,
                               allowParallel = TRUE, 
                               seeds = seeds)
    # We use RF model for the RFE
    caret::rfe(x = t(trainingData[[y]]), 
               y = as.factor(trainingLabs[[y]]),
               sizes = 1:x, 
               method ='rf',
               ntree = round(sqrt(nrow(trainingData[[y]]))),
               # ntree = 250,
               metric= 'accuracy', 
               rfeControl=contr)
    
  })
})
stopCluster(cl)
registerDoSEQ()
names(resultsOverall) <- names(trainingData)
gc()

# Select the best feature set which balances number of features and accuracy
# Mean cross validated accuracy
modelPerf <- lapply(resultsOverall, function(x){
  lapply(x, function(y){mean(y$results$Accuracy)})
})
modelFeat <- lapply(resultsOverall, function(x){
  lapply(x, function(y){y$optVariables})
})

# Compute feature importance
modelFeatImport <- lapply(resultsOverall, function(x){
  lapply(x, function(y){
    df <- as.data.frame(varImp(y)) |>
      rownames_to_column(var = 'FeatureID')
    return(df)
  })
})

# Extract the model performance and number of selected features
modelFeat.df <- trainingData
for(i in 1:length(modelFeat)){
  df <- data.frame(matrix(data = NA, ncol = 5, nrow = 6)); colnames(df) <- c('Omic','FeatSpace','ModelAccuracy','NoFeatures','Features')
  df$Omic <- names(trainingData)[i]
  df$FeatSpace <- sizes
  df1 <- modelPerf[[i]]; df$ModelAccuracy <- do.call('c',df1)
  df2 <- modelFeat[[i]]; df$NoFeatures <- do.call('c', lapply(df2,length))
  for(j in 1:length(df2)){
    df$Features[j] <- toString(df2[[j]])
  }
  modelFeat.df[[i]] <- df
}
modelFeat.df <- do.call('rbind',modelFeat.df)|>
  as.data.frame()
openxlsx::write.xlsx(modelFeat.df, 'FeaturesSelection_results.xlsx', asTable = T)

# Selected features from transcriptome and metabolome data
selectFeats <- list(RNA = resultsOverall[[1]][[3]]$optVariables,
                    Met = resultsOverall[[2]][[2]]$optVariables)
openxlsx::write.xlsx(lapply(selectFeats,as.data.frame), 'FeaturesData.xlsx')


trainingDataFinal <- lapply(seq_along(sizes), function(x){
  feat1 <- resultsOverall[[1]][[x]]$optVariables; feat2 <- resultsOverall[[2]][[x]]$optVariables;
  df1 <- trainingData[[1]]; 
  df1 <- df1[which(rownames(df1) %in% as.character(feat1)),]
  df2 <- trainingData[[2]]; 
  df2 <- df2[which(rownames(df2) %in% as.character(feat2)),]
  return(list(RNA = df1, Met = df2))
})
names(trainingDataFinal) <- paste0(sizes)


testDataFinal <- lapply(seq_along(sizes), function(x){
  feat1 <- resultsOverall[[1]][[x]]$optVariables; feat2 <- resultsOverall[[2]][[x]]$optVariables;
  df1 <- testData[[1]]; 
  df1 <- df1[which(rownames(df1) %in% as.character(feat1)),]
  df2 <- testData[[2]]; 
  df2 <- df2[which(rownames(df2) %in% as.character(feat2)),]
  return(list(RNA = df1, Met = df2))
})
names(testDataFinal) <- paste0(sizes)

predDataFinal <- lapply(seq_along(sizes), function(x){
  feat1 <- resultsOverall[[1]][[x]]$optVariables; feat2 <- resultsOverall[[2]][[x]]$optVariables;
  df1 <- predictionData[[1]]; 
  df1 <- df1[which(rownames(df1) %in% as.character(feat1)),]
  df2 <- predictionData[[2]]; 
  df2 <- df2[which(rownames(df2) %in% as.character(feat2)),]
  return(list(RNA = df1, Met = df2))
})
names(predDataFinal) <- paste0(sizes)

# RF algorithm for classification
# Parameters and data were prepared above for RF classification.
set.seed(7/5/2023)
mlModels <- function(trainData, 
                     trainLabel, 
                     testData, 
                     testLabel=NULL,
                     resampling = 'cv',
                     predData) {
  optParam <- trainControl(method= paste0(resampling), 
                           number=5, 
                           repeats = 5, 
                           search = 'grid',
                           savePredictions=TRUE, 
                           classProbs=TRUE)
  
  # Train the models with RF
  rfFit <- caret::train(t(trainData), trainLabel, method = 'rf', trControl = optParam)
  # Training performance
  ModelTrainPerform <- mean(rfFit$results$Accuracy)
  # Test the models and measure performance
  Models <- list(rf = predict(rfFit, newdata = t(testData)))
  ModelsPredict <- list(rf = predict(rfFit, newdata = t(predData)))
  ModelsPredict <- lapply(ModelsPredict, function(x){
    table(x)
  })
  ModelsPredict <- do.call(rbind,ModelsPredict)|> as.data.frame()
  
  Models <- lapply(Models, function(x){
    caret::confusionMatrix(x,testLabel)$overall
  })
  Models2 <- do.call(rbind,Models) |> as.data.frame()
  Models2$TrainingPerformance <- ModelTrainPerform
  return(list(TrainedModels = list(svmLin = svmLin, knn = knnFit, RF = rfFit),
              ModelPredictions = ModelsPredict,
              ModelTrainPerf = ModelTrainPerform,
              ModelTestPerf = Models2))
}

trainTestModels <- lapply(seq_along(trainingDataFinal), function(x){
  trainData1 <- trainingDataFinal[[x]]; 
  testData1 <- testDataFinal[[x]];
  predData1 <- predDataFinal[[x]];
  return(lapply(seq_along(trainData1),function(y){
    train.Lab <- trainingLabs[[y]]; testLab = testLabs[[y]];
    df.train <- trainData1[[y]]; df.test <- testData1[[y]];
    df.pred<- predData1[[y]]
    mlModels(df.train, train.Lab, df.test, testLab, predData = df.pred)
  }))
})

names(trainTestModels) <- paste0('FeatSet_',names(trainingDataFinal))

lapply(seq_along(trainTestModels), function(x){
  fn <- names(trainTestModels)[x]
  df <- list()
  mlModel <- trainTestModels[[x]]
  mlModel.Perf <- do.call(rbind,lapply(mlModel,function(x){x$ModelTestPerf}))|> 
    as.data.frame()|>
    rownames_to_column(var = 'mlModel')|>
    dplyr::mutate(Omic = c(rep('RNA',3),rep('Metabolome',3)))|>
    dplyr::rename(TrainingAccuracy = 9,TestingAccuracy = 2)|>
    dplyr::mutate(Features = c(rep(length(resultsOverall[[1]][[x]]$optVariables),3),
                               rep(length(resultsOverall[[2]][[x]]$optVariables),3)))
  mlModel.Pred <- do.call(rbind,lapply(mlModel,function(x){x$ModelPredictions}))|> 
    as.data.frame()|>
    rownames_to_column(var = 'mlModel')|>
    dplyr::mutate(Omic = c(rep('RNA',3),rep('Metabolome',3)))
  
  df[['ModelPerformance']] <- cbind(mlModel.Perf,mlModel.Pred)
  openxlsx::write.xlsx(df, file = paste0('Results/FeatSet_',fn,'_Features.xlsx'))
})

save.image('Final.RData')
rm(list = ls())