#This is the code for our project


#Read in the data

input_data <- readRDS(file ="CLL-Bcells_list.RDS.gz")

annotation <- read.csv("sample_annotation.csv")

##Divide the data in 4 subgroups

genes <- input_data$genes

promoters <- input_data$promoters

cpgislands <- input_data$cpgislands

tiling <- input_data$tiling

#Rename our patients in "genes"

colnames(genes)[11:30] <- c("P1_healthy_beta","P2_healthy_beta","P3_healthy_beta","P4_healthy_beta","P5_healthy_beta","P6_CLL_beta","P7_CLL_beta","P8_CLL_beta","P9_CLL_beta","P10_CLL_beta","P1_healthy_coverage","P2_healthy_coverage","P3_healthy_coverage","P4_healthy_coverage","P5_healthy_coverage","P6_CLL_coverage","P7_CLL_coverage","P8_CLL_coverage","P9_CLL_coverage","P10_CLL_coverage")

#Create a data frame containing coverage values

cov_genes <- genes[ ,c(1,21:30)]

##Remove ChrX and ChrY

cov_genes_new <- cov_genes[-which(cov_genes =="chrX"),]
cov_genes <- cov_genes_new[-which(cov_genes_new =="chrY"),]
rm(cov_genes_new)
cov_genes <- cov_genes[ ,c(2:11)]

#Calculate coverage means for each gene

cov_genes_means <- rowMeans(data.matrix(cov_genes))


#logarithmic density plot of coverage values
plot(density(log10(cov_genes_means)), xlab = "coverage means", main = "coverage distribution")

#finding the upper threshold for deleting coverage values by quantiles we read about in literature
#we draw the lines of the quantiles in a logarithmic plot

quantile(cov_genes_means, probs = c(.95))
abline(v=log10(52529.75))

quantile(cov_genes_means, probs = c(.975))
abline(v=log10(82556.55))

quantile(cov_genes_means, probs = c(.999))
abline(v=log10(311230.4))

#we set the lower threshold to the coverage value 10 and 15
abline(v=log10(10))
abline(v=log10(15))


#Nested for loop

##Mean of every gene (check for NAÃÂ´s first and set them to zero if NAÃÂ´s available)

sum(is.na(cov_genes))

cov_genes[is.na(cov_genes)] <- 0

cov_genes_mean <- rowMeans(data.matrix(cov_genes))

##Set threshold

threshold1 <- 15

threshold2 <- quantile(cov_genes_mean, probs = seq(0.90,0.90,0.05))


##Nested for loop 

for(i in 1:nrow(cov_genes)){
  for(j in 1:ncol(cov_genes)){
    if(isTRUE(cov_genes[i,j] <= threshold1)){
      cov_genes[i,j] <- NA
    } 
    if(isTRUE(cov_genes[i,j] >= threshold2)){
      cov_genes[i,j] <- NA
    }
  }
}

rm(i,j,threshold1,threshold2)

##Dataframe for beta values(+Chromosomes)

beta_genes <- genes[,c(1,11:20)]

##Also remove chrX and chrY

beta_genes_new <- beta_genes[-which(beta_genes =="chrX"),]
beta_genes <- beta_genes_new[-which(beta_genes_new =="chrY"),]
rm(beta_genes_new)
beta_genes <- beta_genes [ ,c(2:11)]

##Nested for loop to bring Coverage NAÃÂ´s to beta NAÃÂ´s

for(k in 1:ncol(beta_genes)){
  for(l in 1:nrow(beta_genes)){
    if(isTRUE(is.na(cov_genes[k,l] == TRUE))){
      beta_genes[k,l] <- NA
    } 
  }
}
rm(k,l)

##Count NAÃÂ´s in dataframe

rmv.rows_beta_genes = apply(beta_genes,1, function(x){sum(is.na(x))})

##New dataframe out of beta_genes, where all rows with more than 2 NAÃÂ´s are removed

beta_genes_cleaned <- beta_genes[-which(rmv.rows_beta_genes >2),]

row_difference = nrow(genes)-nrow(beta_genes_cleaned)
sum(row_difference)
genes_deleted_percentage = row_difference/nrow(genes)*100
sum(genes_deleted_percentage)


#replace beta values for 0 and 1
beta_genes_cleaned[beta_genes_cleaned==0]<-0.00000001
beta_genes_cleaned[beta_genes_cleaned==1]<-0.99999999


#create two separate data frames for sick and healthy patients

beta_genes_healthy <- beta_genes_cleaned[,c(1:5)]
beta_genes_cancer <- beta_genes_cleaned[,c(6:10)]


#replace NAs by row means

k <- which(is.na(beta_genes_healthy), arr.ind=TRUE)
beta_genes_healthy[k] <- rowMeans(beta_genes_healthy, na.rm=TRUE)[k[,1]]

l <- which(is.na(beta_genes_cancer), arr.ind=TRUE)
beta_genes_cancer[l] <- rowMeans(beta_genes_cancer, na.rm=TRUE)[l[,1]]


#Normalisation: Transform beta-values into M-values
M_genes_healthy <- log2(beta_genes_healthy/(1-(beta_genes_healthy)))
M_genes_cancer <- log2(beta_genes_cancer/(1-(beta_genes_cancer)))

#dataset containing healthy and cancer M-values
M_genes <- cbind(M_genes_healthy, M_genes_cancer)

#Rename columns to shorter names

colnames(M_genes) <- c("1H","2H","3H","4H","5H","6CLL","7CLL","8CLL","9CLL","10CLL")

#PCA

pca_M <- prcomp(t(M_genes))

#how much variance accounts for the components
var_pca <- pca_M$sdev^2
var_pca_per <- round(var_pca/sum(var_pca)*100, 1)
plot(var_pca_per, main="variation of our data", xlab="Principal Components", ylab="Percent Variation", xlim = c(0,11), ylim=c(0,25),type = "o", pch=20)


#graph of component 1 - 9
#par(mfrow = c(3,3)) ACHTUNG GEFAHR!
#install.packages("gridExtra")
#library(gridExtra)

#graph of component 1 - 9
install.packages("ggplot2")
library(ggplot2)

##With PC1
## 1&2
pca_values2 <- data.frame(Sample=rownames(pca_M$x),
                         X=pca_M$x[,1],
                         Y=pca_M$x[,2])

ggplot(data=pca_values2, aes(x=X, y=Y, label=Sample)) +
  geom_text(aes(colour = annotation$cellTypeGroup)) +
  xlab(paste("PC1 - ", var_pca_per[1], "%", sep="")) +
  ylab(paste("PC2 - ", var_pca_per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("PCA Graph")

## 1&3
pca_values3 <- data.frame(Sample=rownames(pca_M$x),
                         X=pca_M$x[,1],
                         Y=pca_M$x[,3])

ggplot(data=pca_values3, aes(x=X, y=Y, label=Sample)) +
  geom_text(aes(colour = annotation$cellTypeGroup)) +
  xlab(paste("PC1 - ", var_pca_per[1], "%", sep="")) +
  ylab(paste("PC3 - ", var_pca_per[3], "%", sep="")) +
  theme_bw() +
  ggtitle("PCA Graph")

## 1&4
pca_values4 <- data.frame(Sample=rownames(pca_M$x),
                          X=pca_M$x[,1],
                          Y=pca_M$x[,4])

ggplot(data=pca_values4, aes(x=X, y=Y, label=Sample)) +
  geom_text(aes(colour = annotation$cellTypeGroup)) +
  xlab(paste("PC1 - ", var_pca_per[1], "%", sep="")) +
  ylab(paste("PC4 - ", var_pca_per[4], "%", sep="")) +
  theme_bw() +
  ggtitle("PCA Graph")


##With PC 2
## 2&3
pca_values23 <- data.frame(Sample=rownames(pca_M$x),
                          X=pca_M$x[,2],
                          Y=pca_M$x[,3])

ggplot(data=pca_values23, aes(x=X, y=Y, label=Sample)) +
  geom_text(aes(colour = annotation$cellTypeGroup)) +
  xlab(paste("PC2 - ", var_pca_per[2], "%", sep="")) +
  ylab(paste("PC3 - ", var_pca_per[3], "%", sep="")) +
  theme_bw() +
  ggtitle("PCA Graph")

## 2&4
pca_values24 <- data.frame(Sample=rownames(pca_M$x),
                          X=pca_M$x[,2],
                          Y=pca_M$x[,4])

ggplot(data=pca_values24, aes(x=X, y=Y, label=Sample)) +
  geom_text(aes(colour = annotation$cellTypeGroup)) +
  xlab(paste("PC2 - ", var_pca_per[2], "%", sep="")) +
  ylab(paste("PC4 - ", var_pca_per[4], "%", sep="")) +
  theme_bw() +
  ggtitle("PCA Graph")


##With PC 3
## 3&4
pca_values34 <- data.frame(Sample=rownames(pca_M$x),
                           X=pca_M$x[,3],
                           Y=pca_M$x[,4])

ggplot(data=pca_values34, aes(x=X, y=Y, label=Sample)) +
  geom_text(aes(colour = annotation$cellTypeGroup)) +
  xlab(paste("PC3 - ", var_pca_per[3], "%", sep="")) +
  ylab(paste("PC4 - ", var_pca_per[4], "%", sep="")) +
  theme_bw() +
  ggtitle("PCA Graph")



#gg plot divide by shape and colour

ggplot(data=pca_values, aes(x=X, y=Y, label=Sample)) +
  geom_point(aes(shape=annotation$cellTypeGroup, color=annotation$Predicted.Gender)) +
  xlab(paste("PC1 - ", var_pca_per[1], "%", sep="")) +
  ylab(paste("PC2 - ", var_pca_per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("PCA with PC 1 & 2 and check for sex")

#investigation for batch effekt
#PC1&2

ggplot(data=pca_values2, aes(x=X, y=Y, label=Sample)) +
  geom_point(aes(shape=annotation$cellTypeGroup, color=annotation$Predicted.Gender)) +
  xlab(paste("PC1 - ", var_pca_per[1], "%", sep="")) +
  ylab(paste("PC2 - ", var_pca_per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("PCA Graph")

ggplot(data=pca_values2, aes(x=X, y=Y, label=Sample)) +
  geom_point(aes(shape=annotation$cellTypeGroup, color=annotation$FIRST_SUBMISSION_DATE)) +
  xlab(paste("PC1 - ", var_pca_per[1], "%", sep="")) +
  ylab(paste("PC2 - ", var_pca_per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("PCA Graph")

ggplot(data=pca_values2, aes(x=X, y=Y, label=Sample)) +
  geom_point(aes(shape=annotation$cellTypeGroup, color=annotation$SAMPLE_DESC_3)) +
  xlab(paste("PC1 - ", var_pca_per[1], "%", sep="")) +
  ylab(paste("PC2 - ", var_pca_per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("PCA Graph")

ggplot(data=pca_values2, aes(x=X, y=Y, label=Sample)) +
  geom_point(aes(shape=annotation$cellTypeGroup, color=annotation$DONOR_AGE)) +
  xlab(paste("PC1 - ", var_pca_per[1], "%", sep="")) +
  ylab(paste("PC2 - ", var_pca_per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("PCA Graph")

ggplot(data=pca_values2, aes(x=X, y=Y, label=Sample)) +
  geom_point(aes(shape=annotation$cellTypeGroup, color=annotation$BIOMATERIAL_PROVIDER)) +
  xlab(paste("PC1 - ", var_pca_per[1], "%", sep="")) +
  ylab(paste("PC2 - ", var_pca_per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("PCA Graph")


#PC1&3

ggplot(data=pca_values3, aes(x=X, y=Y, label=Sample)) +
  geom_point(aes(shape=annotation$cellTypeGroup, color=annotation$Predicted.Gender)) +
  xlab(paste("PC1 - ", var_pca_per[1], "%", sep="")) +
  ylab(paste("PC3 - ", var_pca_per[3], "%", sep="")) +
  theme_bw() +
  ggtitle("PCA Graph")

ggplot(data=pca_values3, aes(x=X, y=Y, label=Sample)) +
  geom_point(aes(shape=annotation$cellTypeGroup, color=annotation$FIRST_SUBMISSION_DATE)) +
  xlab(paste("PC1 - ", var_pca_per[1], "%", sep="")) +
  ylab(paste("PC3 - ", var_pca_per[3], "%", sep="")) +
  theme_bw() +
  ggtitle("PCA Graph")

ggplot(data=pca_values3, aes(x=X, y=Y, label=Sample)) +
  geom_point(aes(shape=annotation$cellTypeGroup, color=annotation$SAMPLE_DESC_3)) +
  xlab(paste("PC1 - ", var_pca_per[1], "%", sep="")) +
  ylab(paste("PC3 - ", var_pca_per[3], "%", sep="")) +
  theme_bw() +
  ggtitle("PCA Graph")

ggplot(data=pca_values3, aes(x=X, y=Y, label=Sample)) +
  geom_point(aes(shape=annotation$cellTypeGroup, color=annotation$DONOR_AGE)) +
  xlab(paste("PC1 - ", var_pca_per[1], "%", sep="")) +
  ylab(paste("PC3 - ", var_pca_per[3], "%", sep="")) +
  theme_bw() +
  ggtitle("PCA Graph")

ggplot(data=pca_values3, aes(x=X, y=Y, label=Sample)) +
  geom_point(aes(shape=annotation$cellTypeGroup, color=annotation$BIOMATERIAL_PROVIDER)) +
  xlab(paste("PC1 - ", var_pca_per[1], "%", sep="")) +
  ylab(paste("PC3 - ", var_pca_per[3], "%", sep="")) +
  theme_bw() +
  ggtitle("PCA Graph")


#PC1&4

ggplot(data=pca_values4, aes(x=X, y=Y, label=Sample)) +
  geom_point(aes(shape=annotation$cellTypeGroup, color=annotation$Predicted.Gender)) +
  xlab(paste("PC1 - ", var_pca_per[1], "%", sep="")) +
  ylab(paste("PC4 - ", var_pca_per[4], "%", sep="")) +
  theme_bw() +
  ggtitle("PCA Graph")

ggplot(data=pca_values4, aes(x=X, y=Y, label=Sample)) +
  geom_point(aes(shape=annotation$cellTypeGroup, color=annotation$FIRST_SUBMISSION_DATE)) +
  xlab(paste("PC1 - ", var_pca_per[1], "%", sep="")) +
  ylab(paste("PC4 - ", var_pca_per[4], "%", sep="")) +
  theme_bw() +
  ggtitle("PCA Graph")

ggplot(data=pca_values4, aes(x=X, y=Y, label=Sample)) +
  geom_point(aes(shape=annotation$cellTypeGroup, color=annotation$SAMPLE_DESC_3)) +
  xlab(paste("PC1 - ", var_pca_per[1], "%", sep="")) +
  ylab(paste("PC4 - ", var_pca_per[4], "%", sep="")) +
  theme_bw() +
  ggtitle("PCA Graph")

ggplot(data=pca_values4, aes(x=X, y=Y, label=Sample)) +
  geom_point(aes(shape=annotation$cellTypeGroup, color=annotation$DONOR_AGE)) +
  xlab(paste("PC1 - ", var_pca_per[1], "%", sep="")) +
  ylab(paste("PC4 - ", var_pca_per[4], "%", sep="")) +
  theme_bw() +
  ggtitle("PCA Graph")

ggplot(data=pca_values4, aes(x=X, y=Y, label=Sample)) +
  geom_point(aes(shape=annotation$cellTypeGroup, color=annotation$BIOMATERIAL_PROVIDER)) +
  xlab(paste("PC1 - ", var_pca_per[1], "%", sep="")) +
  ylab(paste("PC4 - ", var_pca_per[4], "%", sep="")) +
  theme_bw() +
  ggtitle("PCA Graph")


#PC2&3

ggplot(data=pca_values23, aes(x=X, y=Y, label=Sample)) +
  geom_point(aes(shape=annotation$cellTypeGroup, color=annotation$Predicted.Gender)) +
  xlab(paste("PC2 - ", var_pca_per[2], "%", sep="")) +
  ylab(paste("PC3 - ", var_pca_per[3], "%", sep="")) +
  theme_bw() +
  ggtitle("PCA Graph")

ggplot(data=pca_values23, aes(x=X, y=Y, label=Sample)) +
  geom_point(aes(shape=annotation$cellTypeGroup, color=annotation$FIRST_SUBMISSION_DATE)) +
  xlab(paste("PC2 - ", var_pca_per[2], "%", sep="")) +
  ylab(paste("PC3 - ", var_pca_per[3], "%", sep="")) +
  theme_bw() +
  ggtitle("PCA Graph")

ggplot(data=pca_values23, aes(x=X, y=Y, label=Sample)) +
  geom_point(aes(shape=annotation$cellTypeGroup, color=annotation$SAMPLE_DESC_3)) +
  xlab(paste("PC2 - ", var_pca_per[2], "%", sep="")) +
  ylab(paste("PC3 - ", var_pca_per[3], "%", sep="")) +
  theme_bw() +
  ggtitle("PCA Graph")

ggplot(data=pca_values23, aes(x=X, y=Y, label=Sample)) +
  geom_point(aes(shape=annotation$cellTypeGroup, color=annotation$DONOR_AGE)) +
  xlab(paste("PC2 - ", var_pca_per[2], "%", sep="")) +
  ylab(paste("PC3 - ", var_pca_per[3], "%", sep="")) +
  theme_bw() +
  ggtitle("PCA Graph")

ggplot(data=pca_values23, aes(x=X, y=Y, label=Sample)) +
  geom_point(aes(shape=annotation$cellTypeGroup, color=annotation$BIOMATERIAL_PROVIDER)) +
  xlab(paste("PC2 - ", var_pca_per[2], "%", sep="")) +
  ylab(paste("PC3 - ", var_pca_per[3], "%", sep="")) +
  theme_bw() +
  ggtitle("PCA Graph")


#PC3&4

ggplot(data=pca_values34, aes(x=X, y=Y, label=Sample)) +
  geom_point(aes(shape=annotation$cellTypeGroup, color=annotation$Predicted.Gender)) +
  xlab(paste("PC3 - ", var_pca_per[3], "%", sep="")) +
  ylab(paste("PC4 - ", var_pca_per[4], "%", sep="")) +
  theme_bw() +
  ggtitle("PCA Graph")

ggplot(data=pca_values34, aes(x=X, y=Y, label=Sample)) +
  geom_point(aes(shape=annotation$cellTypeGroup, color=annotation$FIRST_SUBMISSION_DATE)) +
  xlab(paste("PC3 - ", var_pca_per[3], "%", sep="")) +
  ylab(paste("PC4 - ", var_pca_per[4], "%", sep="")) +
  theme_bw() +
  ggtitle("PCA Graph")

ggplot(data=pca_values34, aes(x=X, y=Y, label=Sample)) +
  geom_point(aes(shape=annotation$cellTypeGroup, color=annotation$SAMPLE_DESC_3)) +
  xlab(paste("PC3 - ", var_pca_per[3], "%", sep="")) +
  ylab(paste("PC4 - ", var_pca_per[4], "%", sep="")) +
  theme_bw() +
  ggtitle("PCA Graph")

ggplot(data=pca_values34, aes(x=X, y=Y, label=Sample)) +
  geom_point(aes(shape=annotation$cellTypeGroup, color=annotation$DONOR_AGE)) +
  xlab(paste("PC3 - ", var_pca_per[3], "%", sep="")) +
  ylab(paste("PC4 - ", var_pca_per[4], "%", sep="")) +
  theme_bw() +
  ggtitle("PCA Graph")

ggplot(data=pca_values34, aes(x=X, y=Y, label=Sample)) +
  geom_point(aes(shape=annotation$cellTypeGroup, color=annotation$BIOMATERIAL_PROVIDER)) +
  xlab(paste("PC3 - ", var_pca_per[3], "%", sep="")) +
  ylab(paste("PC4 - ", var_pca_per[4], "%", sep="")) +
  theme_bw() +
  ggtitle("PCA Graph")


#finding the most important 30 genes which have the most influence
loading_scores <- pca_M$rotation[,1]
gene_scores <- abs(loading_scores) 
ranked_gene_score <- sort(gene_scores, decreasing=TRUE)
genes_top_30 <- names(ranked_gene_score[1:30])
View(genes_top_30)

#scores with pos and neg sign
pca_M$rotation[genes_top_30,1]
View(pca_M$rotation[genes_top_30,1])

#Look for elbow in top gene variance
plot(abs(pca_M$rotation[genes_top_30,1]),type = "o", pch=20)




