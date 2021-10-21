#Installing necessary packages
install.packages("BiocManager")
BiocManager::install("limma")
BiocManager::install("edgeR")
library(limma)
library(edgeR)
library(ggpubr)
library(data.table)


#loading datasets
expresData <- read.delim("dgrp_expression_Female.txt", row.names = 1)
eyeSizeData <- read.delim("rpr.txt")


#data cleaning and categorizing in groups
colnames(expresData) <- gsub("line_", "", colnames(expresData))
group1names <- seq(1, length(expresData), 2)
splitData <- expresData[, group1names]
colnames(splitData) <- gsub(".1$", "", colnames(splitData))


#understanding target strains and creating a finalized data frame for both expression and eye sizes data
finalNames <- NULL
for (i in 1:nrow(eyeSizeData)) {
  target <- as.character(eyeSizeData[i,1])
  if(target %in% colnames(splitData)){
    finalNames <- append(finalNames, target)
  }
}
finalExpressData <- splitData[, finalNames]
droppedValues <- NULL
for (i in 1:nrow(eyeSizeData)) {
  target <- as.character(eyeSizeData[i, 1])
  if((target %in% finalNames)==FALSE){
    droppedValues <- append(droppedValues, i)
  }
}
finalSizeData <- eyeSizeData[-droppedValues,]


#Quantile calculation
q <- quantile(finalSizeData$Mean.Eye.Size, probs = c(0, 0.2, 0.8, 1))
groups <- c("small", "medium", "large")
eyeSizes <- cut(finalSizeData$Mean.Eye.Size, breaks = q, include.lowest = TRUE, labels = groups)


#DGElist creation and recalculation of normalization factors
expresList <- DGEList(as.matrix(finalExpressData))
expresList <- calcNormFactors(expresList)


#Eliminating low expression genes
#manually
droppedGenes <- which(apply(cpm(expresList), 1, max) < 35)
#automatically through edgeR
keptGenes <- filterByExpr(expresList, group = finalNames)
#manually calculated value is taken into consideration
expresListFilt <- expresList[-droppedGenes,]


#multidimensional scaling plot using column names as labels
plotMDS(expresList, col = as.numeric(eyeSizes))


#model matrix creation
modelMat <- model.matrix(~0 + eyeSizes)


#usage of voom with the model created to create table for more readability
vObj <- voom(expresListFilt, modelMat, plot = T)
vFit <- lmFit(vObj, design = modelMat)
vFitBayes <- eBayes(vFit)
fitTable <- topTable(vFitBayes)


#contrasting between different eye sizes and new table is created
contrasts <- makeContrasts(eyeSizeslarge - eyeSizesmedium, eyeSizesmedium - eyeSizessmall, eyeSizeslarge - eyeSizessmall, levels = modelMat)
contrFit <- contrasts.fit(vFitBayes, contrasts)
contrFitBayes <- eBayes(contrFit)
contrFitTable <- topTable(contrFitBayes, resort.by = "p")


#transposing data
transData <- rbind(finalExpressData, as.vector(finalSizeData[, 2]))
transCol <- colnames(transData)
transRow <- rownames(transData)
transData <- transpose(transData)
rownames(transData) <- transCol
colnames(transData) <- transRow
colnames(transData)[colnames(transData) == "18141"] <- "eyeSize"


#scatter plot for expressed genes in the table without contrasts
ggscatter(data = transData, x = "eyeSize", y = rownames(fitTable[1,]), add = "reg.line", xlab = "Mean Eye Size", ylab = "Expression Value", title = rownames(fitTable[1,]))
ggscatter(data = transData, x = "eyeSize", y = rownames(fitTable[2,]), add = "reg.line", xlab = "Mean Eye Size", ylab = "Expression Value", title = rownames(fitTable[2,]))
ggscatter(data = transData, x = "eyeSize", y = rownames(fitTable[3,]), add = "reg.line", xlab = "Mean Eye Size", ylab = "Expression Value", title = rownames(fitTable[3,]))


#scatter plot for expressed genes in the table with contrasts
ggscatter(data = transData, x = "eyeSize", y = rownames(contrFitTable[1,]), add = "reg.line", xlab = "Mean Eye Size", ylab = "Expression Value", title = rownames(contrFitTable[1,]))
ggscatter(data = transData, x = "eyeSize", y = rownames(contrFitTable[2,]), add = "reg.line", xlab = "Mean Eye Size", ylab = "Expression Value", title = rownames(contrFitTable[2,]))
ggscatter(data = transData, x = "eyeSize", y = rownames(contrFitTable[3,]), add = "reg.line", xlab = "Mean Eye Size", ylab = "Expression Value", title = rownames(contrFitTable[3,]))
ggscatter(data = transData, x = "eyeSize", y = rownames(contrFitTable[4,]), add = "reg.line", xlab = "Mean Eye Size", ylab = "Expression Value", title = rownames(contrFitTable[4,]))
ggscatter(data = transData, x = "eyeSize", y = rownames(contrFitTable[5,]), add = "reg.line", xlab = "Mean Eye Size", ylab = "Expression Value", title = rownames(contrFitTable[5,]))
ggscatter(data = transData, x = "eyeSize", y = rownames(contrFitTable[6,]), add = "reg.line", xlab = "Mean Eye Size", ylab = "Expression Value", title = rownames(contrFitTable[6,]))


#Table summary
print(fitTable)
summary(fitTable)
print(contrFitTable)
summary(contrFitTable)

#simple observations
newEyeData <- eyeSizeData[order(eyeSizeData$Mean.Eye.Size, decreasing = TRUE), ]
maxes <- NULL
means <- NULL
for (i in 1:nrow(newEyeData)) {
  target <- NULL
  target <- as.character(newEyeData[i,1])
  target <- paste("^", target, ".1$", sep = "")
  maxes <- append(maxes, which.max(expresData[, grep(target, colnames(expresData))]))
  means <- append(means, mean(expresData[, grep(target, colnames(expresData))]))
}
print(maxes)
print(means)
