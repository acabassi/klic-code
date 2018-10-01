##################################
### TCGA breast cancer dataset ###
# Variable selection with glmnet #
##################################

rm(list=ls())

library(glmnet)

# Load datasets
load("Breast.RData")

dim(X[[1]]) # 645 348
dim(X[[2]]) # 574 348
dim(X[[3]]) # 423 348
dim(X[[4]]) # 171 348

# Load PAM50 subtype
survivalBRCA <- read.csv("Table1Nature.csv") 
PAM50 <- survivalBRCA$PAM50.mRNA
names(PAM50) <- survivalBRCA$Complete.TCGA.ID
patientIDs <- colnames(X[[1]])
patientIDs <- substr(patientIDs,1,12)
patientIDs <- gsub(".", "-", patientIDs, fixed = TRUE)
PAM <- PAM50[patientIDs]
table(PAM)

pollo <- function(X, y, alpha = 0){
  cvfit=cv.glmnet(X, y, family = "multinomial", type.multinomial = "grouped", alpha = alpha) 
  plot(cvfit)
  coeff = coef(cvfit, s = "lambda.1se")
  # var = which(coeff$`1` != 0)
  var = which(coeff$`Basal-like`!= 0)
  var = var[2:length(var)] - 1 # remove the intercept
  Xred = X[,var]
  output <- list(var = var, Xred = Xred, nvar = length(var))
}

alpha <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
nVariables <- matrix(0, 4, 6)
output <- list()

# Gene expression data

output$expressionData <- list()

for(i in 2:5){
  output$expressionData[[i]] <-  pollo(t(X[[1]]), PAM, alpha = alpha[i])
  nVariables[1,i] <- output$expressionData[[i]]$nvar
}

# Methylation data

output$methylationData <- list()

for(i in 2:5){
  output$methylationData[[i]] <- pollo(t(X[[2]]), PAM, alpha = alpha[i])
  nVariables[2,i] <- output$methylationData[[i]]$nvar
}

# miRNA data

output$miRNA <- list()

for(i in 2:5){
  output$miRNA[[i]] <- pollo(t(X[[3]]), PAM, alpha = alpha[i])
  nVariables[3,i] <- output$miRNA[[i]]$nvar
}

# Protein expression data

output$protein <- list()

for(i in 2:5){
  output$protein[[i]] <- pollo(t(X[[4]]), PAM, alpha = alpha[i])
  nVariables[4,i] <- output$protein[[i]]$nvar
}

save(output, file = 'allSelectedGlmnet.RData')

### Forse bisognerebbe rifare tutto senza la classe piu' piccola
