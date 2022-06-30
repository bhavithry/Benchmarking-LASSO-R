# Clear environment
##############
rm(list=ls())
gc()
########

# Import libraries
##############
library(dplyr)
library(ROCR)
library(logr)
library(matrixStats)
library(Matrix)
library(glmnet)
library(droplasso)
library(SGL)
library(openxlsx)
library(dendextend)
########

# Import data
##############
setwd("C:/Users/bhavi/OneDrive/Desktop/Research/Data")

# Import data
cells <- read.csv("GSE60749_RnaSeq_single_cell_mESC_TPM.csv")
gene <- cells$Gene


# Transpose and change to data frame
cells <- data.frame(t(data.matrix(cells[,-1])))
colnames(cells) <- gene

class1 <- cells[1:183,]
class2 <- cells[184:267,]
########

# Data Pre-processing
##############
# Shuffle rows
set.seed(100)
class1<- class1[sample(nrow(class1)),]
set.seed(100)
class2<- class2[sample(nrow(class2)),]

class1['y'] <- 1
class2['y'] <- 0
l1 <- nrow(class1)
l2 <- nrow(class2)

class_all <- rbind(class1,class2)

# Remove genes not expressed or equally expressed in all cells
removeZeroVar <- function(df){
  df[, sapply(df, function(x) length(unique(x)) > 1)]
}
data <- removeZeroVar(class_all)

data['y'] <- class_all['y']

rm(class1)
rm(class2)
rm(class_all)
gc()
########

# Pre-allocate variables
##############
df_genes <- data.frame(beta = numeric())

k=1
j=l1+1
nfolds = 10
l1fold = floor(l1/nfolds)
l2fold = floor(l2/nfolds)
#########

# K-Fold Cross validation starts here
##############
for (i in 1:nfolds){
  
  if(i==nfolds){
    test_rows <- c(k:l1,j:(l1+l2))
  } else {
    test_rows <- c(k:(k+l1fold-1),j:(j+l2fold-1))
  }
  
  x_test <- as.matrix(data[test_rows,1:length(data)-1])
  y_test <- as.matrix(data[test_rows, length(data)])
  x_train <- as.matrix(data[-test_rows,1:length(data)-1])
  y_train <- as.matrix(data[-test_rows, length(data)])
  k=k+l1fold
  j=j+l2fold

# alpha = 0, Ridge
##############
set.seed(100)
ridge.cv <- cv.glmnet(x_train, y_train, type.measure="auc", alpha=0, nfolds = 5
                      , family="binomial")
ridge.fit <- glmnet(x_train, y_train, type.measure="auc"
                    , alpha=0, lambda = ridge.cv$lambda.1se
                    , family="binomial")
# Sort and Extract Coefficients and corresponding genes
ridge.coef <- ridge.fit$beta[, which(ridge.fit$lambda==ridge.cv$lambda.1se)]
ridge.coef <- abs(ridge.coef)
ridge.coef <- sort(ridge.coef, decreasing = TRUE)
ridge.genes <- names(ridge.coef)
ngenes <- length(ridge.coef[ridge.coef>=mean(ridge.coef)])
ridge.genes <- ridge.genes[1:ngenes]
########

# alpha = 1, Lasso
##############
set.seed(100)
lasso.cv <- cv.glmnet(x_train, y_train, type.measure="auc", alpha=1, nfolds = 5
                      , family="binomial")
lasso.fit <- glmnet(x_train, y_train, type.measure="auc"
                    , alpha=0, lambda = lasso.cv$lambda.1se
                    , family="binomial")
# Sort and Extract Coefficients and corresponding genes
lasso.coef <- lasso.fit$beta[, which(lasso.fit$lambda==lasso.cv$lambda.1se)]
lasso.coef <- abs(lasso.coef)
lasso.coef <- sort(lasso.coef, decreasing = TRUE)
lasso.genes <- names(lasso.coef)
ngenes <- length(lasso.coef[lasso.coef>=mean(lasso.coef)])
lasso.genes <- lasso.genes[1:ngenes]
########

# alpha = 0.5, Elasticnet
##############
set.seed(100)
elnet.cv <- cv.glmnet(x_train, y_train, type.measure="auc", alpha=0.5, nfolds = 5
                      , family="binomial")
elnet.fit <- glmnet(x_train, y_train, type.measure="auc"
                    , alpha=0, lambda = elnet.cv$lambda.1se
                    , family="binomial")
# Sort and Extract Coefficients and corresponding genes
elnet.coef <- elnet.fit$beta[, which(elnet.fit$lambda==elnet.cv$lambda.1se)]
elnet.coef <- abs(elnet.coef)
elnet.coef <- sort(elnet.coef, decreasing = TRUE)
elnet.genes <- names(elnet.coef)
ngenes <- length(elnet.coef[elnet.coef>=mean(elnet.coef)])
elnet.genes <- elnet.genes[1:ngenes]
########

# Drop Lasso
##############
set.seed(100)
droplasso.fit <- droplasso(x_train, y_train, family="binomial"
                           , keep_prob = 0.5
                           , lambda=0.001)
# Sort and Extract Coefficients and corresponding genes
droplasso.coef <- droplasso.fit$beta[,1]
names(droplasso.coef) <- colnames(x_train)
droplasso.coef <- abs(droplasso.coef)
droplasso.coef <- sort(droplasso.coef, decreasing = TRUE)
droplasso.genes <- names(droplasso.coef)
ngenes <- length(droplasso.coef[droplasso.coef>=mean(droplasso.coef)])
droplasso.genes <- droplasso.genes[1:ngenes]
########

# Create Gene pool 1 (Union of 4 methods)
##############
gene_pool1 <- unique(c(ridge.genes,lasso.genes,elnet.genes,droplasso.genes))
########

# Hierarchical Clustering with selected data
##############
Z.train <- x_train[,gene_pool1]
y.train <- y_train
Z.test <- x_test[,gene_pool1]
y.test <- y_test

cells_df_sc <- as.data.frame(scale(t(Z.train)))
dist_mat <- dist(cells_df_sc, method = 'euclidean')
hclust_avg <- hclust(dist_mat, method = 'average')
cut_avg <- cutree(hclust_avg, k = 200)

Z.train.sd <- colSds(Z.train)
Z.train    <- Z.train[, which(Z.train.sd > 1e-7)]
Z.test     <- Z.test[, which(Z.train.sd > 1e-7)]
grp        <- cut_avg[which(Z.train.sd > 1e-7)]
Z.train.sd <- Z.train.sd[which(Z.train.sd > 1e-7)]
########

# SGL
##############
Z.train <- scale(Z.train, center = TRUE, scale = Z.train.sd)
data.SGL <- list(x = Z.train, y = y.train)
set.seed(100)
sgl.fit <- SGL(data = data.SGL, index = grp, type = "logit")
m <- 6
sgl_lam <- numeric(m)
for(l in 1:m) {
  sgl.probabilities <- predictSGL(sgl.fit, lam=l, newX = Z.test) 
  sgl.y_pred <- ifelse(sgl.probabilities > 0.5, 1, 0)
  sgl.pred <- prediction(sgl.y_pred, y_test)
  sgl.auc.perf <- performance(sgl.pred, measure = "auc")
  sgl_lam[l] <- sgl.auc.perf@y.values
}
sgl.val <- max(unlist(sgl_lam))
print(sgl.val)
########


# Getting top genes
##############
beta <- sgl.fit$beta[,2]
sgl.coef <- data.frame(beta)
row.names(sgl.coef) <- colnames(Z.test)
sgl.coef <- abs(sgl.coef)
########

# Create Dataframe of top genes
###############
if (i==1){
df_genes<- sgl.coef
}
df_genes <- cbind(df_genes,data.frame(sgl.coef))
########

}
########

# Plot and print top genes
##############
top_genes <- data.frame(genes=row.names(df_genes),beta=rowSums(df_genes)/nfolds)
top_genes <- arrange(top_genes, -beta)
plot(top_genes$beta[1:10], xaxt="n"
     ,xlab="Genes"
     ,ylab="Sparse Group Lasso Coefficient")
xtick<-top_genes$genes[1:10]
axis(side=1, at=1:10, labels = xtick
     , las=2, cex.axis=0.7)

write.xlsx(as.data.frame(top_genes$genes[1:20])
           , file = "New_Pipeline_Genes_60749.xlsx", overwrite = TRUE)

#######################
# K-Means Clustering
#######################
library(readxl)
# Import data
##############
setwd("C:/Users/bhavi/OneDrive/Desktop/Research/Data")

# Import data
cells <- read.csv("GSE60749_RnaSeq_single_cell_mESC_TPM.csv")
gene <- cells$Gene

# Transpose and change to data frame
cells <- data.frame(t(data.matrix(cells[,-1])))
colnames(cells) <- gene

# Import top genes
genes <- data.frame(read_xlsx("New_Pipeline_Genes_60749.xlsx"))
genes <- unlist(list(genes[,1]))

# Subset data with top genes
cells <- cells[genes]

cl <- kmeans(cells, centers = 2, nstart = 25)

#par(mar=c(1,1,1,1))
plot(cells[,1], cells[,3]
     , col = cl$cluster
     , xlab = "piRNA - 44441 (TPM)"
     , ylab = "piRNA - 44260 (TPM)"
     , main = "GSE60749 - Pluripotent Stem Cells of Mus musculus"
      )
legend("topright",legend=c("mESCs Cell", "Dgcr8 -/- mESCs Cell")
       ,bty="n", pch = c(1,1)
       ,col=c("red", "black")
       ,cex = 0.9
       ,text.font=2
       ,title = "Cell Types"
      )
