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

##Mean of every gene (check for NA????????????????s first and set them to zero if NA????????????????s available)

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

##Nested for loop to bring Coverage NAs to beta NAs

for(k in 1:ncol(beta_genes)){
  for(l in 1:nrow(beta_genes)){
    if(isTRUE(is.na(cov_genes[k,l] == TRUE))){
      beta_genes[k,l] <- NA
    } 
  }
}
rm(k,l)

##Count NA????????????????s in dataframe

rmv.rows_beta_genes = apply(beta_genes,1, function(x){sum(is.na(x))})

##New dataframe out of beta_genes, where all rows with more than 2 NAs are removed

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
cat <- M_genes[1:10, 1:10]
#Rename columns to shorter names

colnames(M_genes) <- c("1H","2H","3H","4H","5H","6CLL","7CLL","8CLL","9CLL","10CLL")

#PCA

pca_M <- prcomp(t(M_genes))

#how much variance accounts for the components
var_pca <- pca_M$sdev^2
var_pca_per <- round(var_pca/sum(var_pca)*100, 1)
plot(var_pca_per, main="Variation of our data explained by PCs", xlab="Principal Components", ylab="Percentage of variation explained by PC", xlim = c(0,11), ylim=c(0,25),type = "o", pch=20)
##There is not a clear kink in the first PC, therefor we actually have to investigate on all 9 PCs, but the last ones don't explain a lot of variance anyway so we just have a look at the first 4 PCs.



#graph of component 1 - 4
#par(mfrow = c(3,3)) ACHTUNG GEFAHR! Hier muss noch ein guter code gefunden werden um grafen sch??n nebeneinander zu pr??sentieren
#install.packages("gridExtra")
#library(gridExtra)

#graph of component 1 - 4
##install.packages("ggplot2")
library(ggplot2)

##With PC1
## 1&2, den mit ins Markdown nehmen
pca_values2 <- data.frame(Sample=rownames(pca_M$x),
                          X=pca_M$x[,1],
                          Y=pca_M$x[,2])

ggplot(data=pca_values2, aes(x=X, y=Y, label=Sample)) +
  geom_text(aes(colour = annotation$DISEASE)) +
  xlab(paste("PC1 - ", var_pca_per[1], "%", sep="")) +
  ylab(paste("PC2 - ", var_pca_per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("Graph of PC 1&2")

## 1&3
pca_values3 <- data.frame(Sample=rownames(pca_M$x),
                          X=pca_M$x[,1],
                          Y=pca_M$x[,3])

ggplot(data=pca_values3, aes(x=X, y=Y, label=Sample)) +
  geom_text(aes(colour = annotation$DISEASE)) +
  xlab(paste("PC1 - ", var_pca_per[1], "%", sep="")) +
  ylab(paste("PC3 - ", var_pca_per[3], "%", sep="")) +
  theme_bw() +
  ggtitle("PCA Graph")

## 1&4
pca_values4 <- data.frame(Sample=rownames(pca_M$x),
                          X=pca_M$x[,1],
                          Y=pca_M$x[,4])

ggplot(data=pca_values4, aes(x=X, y=Y, label=Sample)) +
  geom_text(aes(colour = annotation$DISEASE)) +
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
  geom_text(aes(colour = annotation$DISEASE)) +
  xlab(paste("PC2 - ", var_pca_per[2], "%", sep="")) +
  ylab(paste("PC3 - ", var_pca_per[3], "%", sep="")) +
  theme_bw() +
  ggtitle("PCA Graph")

## 2&4
pca_values24 <- data.frame(Sample=rownames(pca_M$x),
                           X=pca_M$x[,2],
                           Y=pca_M$x[,4])

ggplot(data=pca_values24, aes(x=X, y=Y, label=Sample)) +
  geom_text(aes(colour = annotation$DISEASE)) +
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
  geom_text(aes(colour = annotation$DISEASE)) +
  xlab(paste("PC3 - ", var_pca_per[3], "%", sep="")) +
  ylab(paste("PC4 - ", var_pca_per[4], "%", sep="")) +
  theme_bw() +
  ggtitle("PCA Graph")



#investigation for batch effekt
#gg plot divide by shape and colour
#PC1&2

ggplot_1 <- ggplot(data=pca_values2, aes(x=X, y=Y, label=Sample)) +
  geom_point(aes(shape=annotation$cellTypeGroup, color=annotation$Predicted.Gender)) +
  xlab(paste("PC1 - ", var_pca_per[1], "%", sep="")) +
  ylab(paste("PC2 - ", var_pca_per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("PC1&2 check for gender")


ggplot_2 <- ggplot(data=pca_values2, aes(x=X, y=Y, label=Sample)) +
  geom_point(aes(shape=annotation$cellTypeGroup, color=annotation$SAMPLE_DESC_3)) +
  xlab(paste("PC1 - ", var_pca_per[1], "%", sep="")) +
  ylab(paste("PC2 - ", var_pca_per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("PC1&2 check for cell type/origin")

#auf das alter testen wir die signifikanz des batch effektes sp??ter nicht, denn eigentlich ist alter kein batch effekt, sonden sorgt f??r Unterschiede, die biologisch bedingt sind.
ggplot_4 <- ggplot(data=pca_values2, aes(x=X, y=Y, label=Sample)) +
  geom_point(aes(shape=annotation$cellTypeGroup, color=annotation$DONOR_AGE)) +
  xlab(paste("PC1 - ", var_pca_per[1], "%", sep="")) +
  ylab(paste("PC2 - ", var_pca_per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("PC1&2 check for age")

ggplot_3 <- ggplot(data=pca_values2, aes(x=X, y=Y, label=Sample)) +
  geom_point(aes(shape=annotation$cellTypeGroup, color=annotation$BIOMATERIAL_PROVIDER)) +
  xlab(paste("PC1 - ", var_pca_per[1], "%", sep="")) +
  ylab(paste("PC2 - ", var_pca_per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("PC1&2 check for biomaterial provider")

#install.packages("gridExtra")
library(gridExtra)
library(ggplot2)
grid.arrange(ggplot_1,ggplot_2, ncol=2)
grid.arrange(ggplot_3,ggplot_4, ncol=2)

#PC1&3

#ggplot(data=pca_values3, aes(x=X, y=Y, label=Sample)) +
# geom_point(aes(shape=annotation$cellTypeGroup, color=annotation$Predicted.Gender)) +
# xlab(paste("PC1 - ", var_pca_per[1], "%", sep="")) +
# ylab(paste("PC3 - ", var_pca_per[3], "%", sep="")) +
# theme_bw() +
#  ggtitle("PC1&3 check for gender")


#ggplot(data=pca_values3, aes(x=X, y=Y, label=Sample)) +
# geom_point(aes(shape=annotation$cellTypeGroup, color=annotation$SAMPLE_DESC_3)) +
#xlab(paste("PC1 - ", var_pca_per[1], "%", sep="")) +
#ylab(paste("PC3 - ", var_pca_per[3], "%", sep="")) +
#theme_bw() +
#ggtitle("PC1&3 check for cell type/origin")

#ggplot(data=pca_values3, aes(x=X, y=Y, label=Sample)) +
# geom_point(aes(shape=annotation$cellTypeGroup, color=annotation$DONOR_AGE)) +
#xlab(paste("PC1 - ", var_pca_per[1], "%", sep="")) +
#ylab(paste("PC3 - ", var_pca_per[3], "%", sep="")) +
#theme_bw() +
#ggtitle("PC1&3 check for age")

#ggplot(data=pca_values3, aes(x=X, y=Y, label=Sample)) +
#geom_point(aes(shape=annotation$cellTypeGroup, color=annotation$BIOMATERIAL_PROVIDER)) +
#xlab(paste("PC1 - ", var_pca_per[1], "%", sep="")) +
#ylab(paste("PC3 - ", var_pca_per[3], "%", sep="")) +
#theme_bw() +
#ggtitle("PC1&3 check for biomaterial provider")


#PC1&4

#ggplot(data=pca_values4, aes(x=X, y=Y, label=Sample)) +
#geom_point(aes(shape=annotation$cellTypeGroup, color=annotation$Predicted.Gender)) +
#xlab(paste("PC1 - ", var_pca_per[1], "%", sep="")) +
#ylab(paste("PC4 - ", var_pca_per[4], "%", sep="")) +
#theme_bw() +
#ggtitle("PC1&4 check for gender")

#ggplot(data=pca_values4, aes(x=X, y=Y, label=Sample)) +
#geom_point(aes(shape=annotation$cellTypeGroup, color=annotation$SAMPLE_DESC_3)) +
#xlab(paste("PC1 - ", var_pca_per[1], "%", sep="")) +
#ylab(paste("PC4 - ", var_pca_per[4], "%", sep="")) +
#theme_bw() +
#ggtitle("PC1&4 check for cell type/origin")

#ggplot(data=pca_values4, aes(x=X, y=Y, label=Sample)) +
#geom_point(aes(shape=annotation$cellTypeGroup, color=annotation$DONOR_AGE)) +
#xlab(paste("PC1 - ", var_pca_per[1], "%", sep="")) +
#ylab(paste("PC4 - ", var_pca_per[4], "%", sep="")) +
#theme_bw() +
#ggtitle("PC1&4 check for age")

#ggplot(data=pca_values4, aes(x=X, y=Y, label=Sample)) +
#geom_point(aes(shape=annotation$cellTypeGroup, color=annotation$BIOMATERIAL_PROVIDER)) +
#xlab(paste("PC1 - ", var_pca_per[1], "%", sep="")) +
#ylab(paste("PC4 - ", var_pca_per[4], "%", sep="")) +
#theme_bw() +
#ggtitle("PC1&4 check for biomaterial provider")


#PC2&3

#ggplot(data=pca_values23, aes(x=X, y=Y, label=Sample)) +
#geom_point(aes(shape=annotation$cellTypeGroup, color=annotation$Predicted.Gender)) +
#xlab(paste("PC2 - ", var_pca_per[2], "%", sep="")) +
#ylab(paste("PC3 - ", var_pca_per[3], "%", sep="")) +
#theme_bw() +
#ggtitle("PC2&3 check for gender")


#ggplot(data=pca_values23, aes(x=X, y=Y, label=Sample)) +
#geom_point(aes(shape=annotation$cellTypeGroup, color=annotation$SAMPLE_DESC_3)) +
#xlab(paste("PC2 - ", var_pca_per[2], "%", sep="")) +
#ylab(paste("PC3 - ", var_pca_per[3], "%", sep="")) +
#theme_bw() +
#ggtitle("PC2&3 check for cell type/origin")

#ggplot(data=pca_values23, aes(x=X, y=Y, label=Sample)) +
#geom_point(aes(shape=annotation$cellTypeGroup, color=annotation$DONOR_AGE)) +
#xlab(paste("PC2 - ", var_pca_per[2], "%", sep="")) +
#ylab(paste("PC3 - ", var_pca_per[3], "%", sep="")) +
#theme_bw() +
#ggtitle("PC2&3 check for age")

#ggplot(data=pca_values23, aes(x=X, y=Y, label=Sample)) +
#geom_point(aes(shape=annotation$cellTypeGroup, color=annotation$BIOMATERIAL_PROVIDER)) +
#xlab(paste("PC2 - ", var_pca_per[2], "%", sep="")) +
#ylab(paste("PC3 - ", var_pca_per[3], "%", sep="")) +
#theme_bw() +
#ggtitle("PC2&3 check for biomaterial provider")


#PC3&4

#ggplot(data=pca_values34, aes(x=X, y=Y, label=Sample)) +
#geom_point(aes(shape=annotation$cellTypeGroup, color=annotation$Predicted.Gender)) +
#xlab(paste("PC3 - ", var_pca_per[3], "%", sep="")) +
#ylab(paste("PC4 - ", var_pca_per[4], "%", sep="")) +
#theme_bw() +
#ggtitle("PC3&4 check for gender")


#ggplot(data=pca_values34, aes(x=X, y=Y, label=Sample)) +
#geom_point(aes(shape=annotation$cellTypeGroup, color=annotation$SAMPLE_DESC_3)) +
#xlab(paste("PC3 - ", var_pca_per[3], "%", sep="")) +
#ylab(paste("PC4 - ", var_pca_per[4], "%", sep="")) +
#theme_bw() +
#ggtitle("PC3&4 check for cell type/origin")

#ggplot(data=pca_values34, aes(x=X, y=Y, label=Sample)) +
#geom_point(aes(shape=annotation$cellTypeGroup, color=annotation$DONOR_AGE)) +
#xlab(paste("PC3 - ", var_pca_per[3], "%", sep="")) +
#ylab(paste("PC4 - ", var_pca_per[4], "%", sep="")) +
#theme_bw() +
#ggtitle("PC3&4 check for age")

#ggplot(data=pca_values34, aes(x=X, y=Y, label=Sample)) +
#geom_point(aes(shape=annotation$cellTypeGroup, color=annotation$BIOMATERIAL_PROVIDER)) +
#xlab(paste("PC3 - ", var_pca_per[3], "%", sep="")) +
#ylab(paste("PC4 - ", var_pca_per[4], "%", sep="")) +
#theme_bw() +
#ggtitle("PC3&4 check for biomaterial provider")



#checking our PC 1 to 5 for significant batch effect per category 

#create dataframes containing PC1-5 and variable of annotation we want to get p value of
#depending on the type of category (numbers, 2 categories, more than 2 different categories, aber im Markdown etwas ausf??hrlicher erkl??ren)
x <- pca_M[["x"]]
pca1_5 <- x[,1:5]
batch_kruskal <- data.frame(pca1_5, annotation$SAMPLE_DESC_3)
batch_wilcoxon <- data.frame(pca1_5, annotation$BIOMATERIAL_PROVIDER, annotation$DISEASE, annotation$Predicted.Gender)
batch_permutation <- data.frame (pca1_5, annotation$SEQ_RUNS_COUNT)

#Permutation test for numbers
#was ist mit der histogramm zeile? wir brauchen doch kein histogram an der stelle oder ?

cor.perm <- function (x, y, nperm = 1000)
{
  r.obs <- cor (x = x, y = y)
  P.par <- cor.test (x = x, y = y)$p.value
  #  r.per <- replicate (nperm, expr = cor (x = x, y = sample (y)))
  r.per <- sapply (1:nperm, FUN = function (i) cor (x = x, y = sample (y)))
  r.per <- c(r.per, r.obs)
  
  P.per <- sum (abs (r.per) >= abs (r.obs))/(nperm + 1) 
  return(list(r.obs = r.obs, P.par = P.par, P.per = P.per))
}

x <-batch_permutation$PC1
y <-batch_permutation$annotation.SEQ_RUNS_COUNT

p_PC1_SEQ <- cor.perm(x,y)

x <-batch_permutation$PC2
p_PC2_SEQ <- cor.perm(x,y)

x <-batch_permutation$PC3
p_PC3_SEQ <- cor.perm(x,y)

x <-batch_permutation$PC4
p_PC4_SEQ <- cor.perm(x,y)

x <-batch_permutation$PC5
p_PC5_SEQ <- cor.perm(x,y)

p_SEQ_RUNS_COUNT <- data.frame(p_PC1_SEQ$P.per, p_PC2_SEQ$P.per, p_PC3_SEQ$P.per, p_PC4_SEQ$P.per, p_PC5_SEQ$P.per)

#Wilcoxon test for 2 categories

pc_1_BIOMATERIAL_PROVIDER <- wilcox.test(batch_wilcoxon$PC1 ~ batch_wilcoxon$annotation.BIOMATERIAL_PROVIDER, data = batch_wilcoxon, exact = FALSE)
pc_2_BIOMATERIAL_PROVIDER <- wilcox.test(batch_wilcoxon$PC2 ~ batch_wilcoxon$annotation.BIOMATERIAL_PROVIDER, data = batch_wilcoxon, exact = FALSE)
pc_3_BIOMATERIAL_PROVIDER <- wilcox.test(batch_wilcoxon$PC3 ~ batch_wilcoxon$annotation.BIOMATERIAL_PROVIDER, data = batch_wilcoxon, exact = FALSE)
pc_4_BIOMATERIAL_PROVIDER <- wilcox.test(batch_wilcoxon$PC4 ~ batch_wilcoxon$annotation.BIOMATERIAL_PROVIDER, data = batch_wilcoxon, exact = FALSE)
pc_5_BIOMATERIAL_PROVIDER <- wilcox.test(batch_wilcoxon$PC5 ~ batch_wilcoxon$annotation.BIOMATERIAL_PROVIDER, data = batch_wilcoxon, exact = FALSE)

pc_1_BIOMATERIAL_PROVIDER <- pc_1_BIOMATERIAL_PROVIDER$p.value
pc_2_BIOMATERIAL_PROVIDER <- pc_2_BIOMATERIAL_PROVIDER$p.value
pc_3_BIOMATERIAL_PROVIDER <- pc_3_BIOMATERIAL_PROVIDER$p.value
pc_4_BIOMATERIAL_PROVIDER <- pc_4_BIOMATERIAL_PROVIDER$p.value
pc_5_BIOMATERIAL_PROVIDER <- pc_5_BIOMATERIAL_PROVIDER$p.value

pc_1_DISEASE <- wilcox.test(batch_wilcoxon$PC1 ~ batch_wilcoxon$annotation.DISEASE, data = batch_wilcoxon, exact = FALSE)
pc_2_DISEASE <- wilcox.test(batch_wilcoxon$PC2 ~ batch_wilcoxon$annotation.DISEASE, data = batch_wilcoxon, exact = FALSE)
pc_3_DISEASE <- wilcox.test(batch_wilcoxon$PC3 ~ batch_wilcoxon$annotation.DISEASE, data = batch_wilcoxon, exact = FALSE)
pc_4_DISEASE <- wilcox.test(batch_wilcoxon$PC4 ~ batch_wilcoxon$annotation.DISEASE, data = batch_wilcoxon, exact = FALSE)
pc_5_DISEASE <- wilcox.test(batch_wilcoxon$PC5 ~ batch_wilcoxon$annotation.DISEASE, data = batch_wilcoxon, exact = FALSE)

pc_1_DISEASE <- pc_1_DISEASE$p.value
pc_2_DISEASE <- pc_2_DISEASE$p.value
pc_3_DISEASE <- pc_3_DISEASE$p.value
pc_4_DISEASE <- pc_4_DISEASE$p.value
pc_5_DISEASE <- pc_5_DISEASE$p.value

pc_1_Predicted.Gender <- wilcox.test(batch_wilcoxon$PC1 ~ batch_wilcoxon$annotation.Predicted.Gender, data = batch_wilcoxon, exact = FALSE)
pc_2_Predicted.Gender <- wilcox.test(batch_wilcoxon$PC2 ~ batch_wilcoxon$annotation.Predicted.Gender, data = batch_wilcoxon, exact = FALSE)
pc_3_Predicted.Gender <- wilcox.test(batch_wilcoxon$PC3 ~ batch_wilcoxon$annotation.Predicted.Gender, data = batch_wilcoxon, exact = FALSE)
pc_4_Predicted.Gender <- wilcox.test(batch_wilcoxon$PC4 ~ batch_wilcoxon$annotation.Predicted.Gender, data = batch_wilcoxon, exact = FALSE)
pc_5_Predicted.Gender <- wilcox.test(batch_wilcoxon$PC5 ~ batch_wilcoxon$annotation.Predicted.Gender, data = batch_wilcoxon, exact = FALSE)

pc_1_Predicted.Gender <- pc_1_Predicted.Gender$p.value
pc_2_Predicted.Gender <- pc_2_Predicted.Gender$p.value
pc_3_Predicted.Gender <- pc_3_Predicted.Gender$p.value
pc_4_Predicted.Gender <- pc_4_Predicted.Gender$p.value
pc_5_Predicted.Gender <- pc_5_Predicted.Gender$p.value

p_DISEASE <- data.frame(pc_1_BIOMATERIAL_PROVIDER, pc_2_BIOMATERIAL_PROVIDER, pc_3_BIOMATERIAL_PROVIDER, pc_4_BIOMATERIAL_PROVIDER, pc_5_BIOMATERIAL_PROVIDER)
p_BIOMATERIAL_PROVIDER <- data.frame(pc_1_BIOMATERIAL_PROVIDER, pc_2_BIOMATERIAL_PROVIDER, pc_3_BIOMATERIAL_PROVIDER, pc_4_BIOMATERIAL_PROVIDER, pc_5_BIOMATERIAL_PROVIDER)
p_Predicted.Gender <- data.frame(pc_1_Predicted.Gender, pc_2_Predicted.Gender, pc_3_Predicted.Gender, pc_4_Predicted.Gender, pc_5_Predicted.Gender)


#Kruskal-Wallis for several categories
batch_kruskal <- within(batch_kruskal, {
  PC1 <- as.numeric(as.character(PC1))
})

batch_kruskal <- within(batch_kruskal, {
  PC2 <- as.numeric(as.character(PC2))
})

batch_kruskal <- within(batch_kruskal, {
  PC3 <- as.numeric(as.character(PC3))
})

batch_kruskal <- within(batch_kruskal, {
  PC4 <- as.numeric(as.character(PC4))
})

batch_kruskal <- within(batch_kruskal, {
  PC5 <- as.numeric(as.character(PC5))
})

sample_desc_3_pc1 <- kruskal.test(batch_kruskal$PC1 ~ batch_kruskal$annotation.SAMPLE_DESC_3, data = batch_kruskal)
sample_desc_3_pc2 <- kruskal.test(batch_kruskal$PC2 ~ batch_kruskal$annotation.SAMPLE_DESC_3, data = batch_kruskal)
sample_desc_3_pc3 <- kruskal.test(batch_kruskal$PC3 ~ batch_kruskal$annotation.SAMPLE_DESC_3, data = batch_kruskal)
sample_desc_3_pc4 <- kruskal.test(batch_kruskal$PC4 ~ batch_kruskal$annotation.SAMPLE_DESC_3, data = batch_kruskal)
sample_desc_3_pc5 <- kruskal.test(batch_kruskal$PC5 ~ batch_kruskal$annotation.SAMPLE_DESC_3, data = batch_kruskal)

p_SAMPLE_DESC_3 <- data.frame(sample_desc_3_pc1$p.value, sample_desc_3_pc2$p.value, sample_desc_3_pc3$p.value, sample_desc_3_pc4$p.value, sample_desc_3_pc5$p.value)




#ansatz eines verkr??ppelten nicht funktionierenden loops
##kruskal_list <- batch_kruskal[,c(1:5)]
##batch_sample_desc_3 <- for (i in 1:5){kruskal.test(kruskal_list[,i] ~ batch_kruskal$annotation.SAMPLE_DESC_3)}

#dataframe erstellen, das alle p values der Kategorien, die wir auf einen Batch Effekt untersuchen, enth??lt

p_DISEASE_t <- as.data.frame(t(p_DISEASE))
p_BIOMATERIAL_PROVIDER_t <- as.data.frame(t(p_BIOMATERIAL_PROVIDER))
p_Predicted.Gender_t <- as.data.frame(t(p_Predicted.Gender))
p_SEQ_RUNS_COUNT_t <- as.data.frame(t(p_SEQ_RUNS_COUNT))
p_SAMPLE_DESC_3_t <- as.data.frame(t( p_SAMPLE_DESC_3))

total_pvalue <- cbind(p_DISEASE_t,p_BIOMATERIAL_PROVIDER_t,p_Predicted.Gender_t,p_SEQ_RUNS_COUNT_t,p_SAMPLE_DESC_3_t)

#give columns informative names
names(total_pvalue)[1] <- "p_DISEASE"
names(total_pvalue)[2] <- "p_BIOMATERIAL"
names(total_pvalue)[3] <- "p_GENDER"
names(total_pvalue)[4] <- "p_SEQ_RUNS_COUNT"
names(total_pvalue)[5] <- "p_SAMPLE_DESC"

#rename rows
rownames(total_pvalue)[rownames(total_pvalue) == "pc_1_BIOMATERIAL_PROVIDER"] <- "PC1"
rownames(total_pvalue)[rownames(total_pvalue) == "pc_2_BIOMATERIAL_PROVIDER"] <- "PC2"
rownames(total_pvalue)[rownames(total_pvalue) == "pc_3_BIOMATERIAL_PROVIDER"] <- "PC3"
rownames(total_pvalue)[rownames(total_pvalue) == "pc_4_BIOMATERIAL_PROVIDER"] <- "PC4"
rownames(total_pvalue)[rownames(total_pvalue) == "pc_5_BIOMATERIAL_PROVIDER"] <- "PC5"

#transform dataframe into a matrix
total_pvalue <- data.matrix(total_pvalue)

#heatmap 

##preparation for heatmap
###p values under 0.1 are significant. to make visualisation easier we set values above to 10 and values equal or below to 1. So we just have to find colours for 2 values.

for (i in 1:ncol(total_pvalue)){
  for (j in 1:nrow(total_pvalue)){
    if(isTRUE(total_pvalue[i,j] > 0.1)){
      total_pvalue[i,j] <- 10
    }
    if (isTRUE(total_pvalue[i,j] <= 0.1)){
      total_pvalue [i,j] <- 1
    }
  }
}


##heatmap 
library(gplots)
heatmap.2(total_pvalue, main = "batch effect in our principal components", Colv = NA, Rowv = NA, dendrogram = "none", sepwidth = c(0.01, 0.01), sepcolor = "black", trace= "none", col = colorRampPalette(c("salmon","blue")), breaks = c(seq(0, 1, length = 2),seq(1.1, 10, length = 2)),colsep = 1:ncol(total_pvalue), rowsep = 1:nrow(total_pvalue), key=FALSE)
##Principal component 1 has a batch effect which you can see looking at the significant p values or at the heatmap
##Therefor we go on with our analysis using Principal component 2 which does't show a batch effect at all.

#create a dataframes with loading scores of the 2000 and 10000 most importat genes, which have the most influence for principal component 2.
loading_scores <- pca_M$rotation[,2]
gene_scores <- abs(loading_scores) 
ranked_gene_score <- sort(gene_scores, decreasing=TRUE)
genes_top_10000 <- names(ranked_gene_score[1:10000])
genes_top_2000 <- names(ranked_gene_score[1:2000])

#showing the most important 10000 genes in PC2 in ranked order with loading scores with pos and neg sign
pca_M$rotation[genes_top_10000,2]
View(pca_M$rotation[genes_top_10000,2])

#Look for elbow in top gene variance
plot(abs(pca_M$rotation[genes_top_10000,2]),type = "o", pch=20)

#Actually there is a kink in the ellbow plot of loading scores of the genes at round about 2000 genes. Therefor we do k-means clustering with 2000 genes.But using only 2000 genes later gives us only 3 significantly differentially methylated genes, which is not sufficient. 10000 gives us 31 genes, which is an amount we can work with. 

#build a dataframe containing all M-values from the 2000 and 10000 most important genes for principal component 2                
k_means_data2 <- M_genes[genes_top_2000,]
k_means_data <- M_genes[genes_top_10000,]
                   
#clustering the patient with k-means using a dataframe containing M values of the 2000 most important genes
k <- kmeans(x = t(k_means_data2), centers = 2)
k <-
  kmeans(
    x = t(k_means_data2),
    centers = 2,
    iter.max = 1000
  )
View(k)

#t-test (students t-test)

cluster <- k[["cluster"]]
p_value = NULL
for(i in 1:nrow(k_means_data)){
  x = k_means_data[i, 1:5]
  y = k_means_data[i, 6:10]
  t_test = t.test(x, y)
  p_value[i] = t_test$p.value
}

#create a dataframe wich contains the p-values from t-test 
pvalues <- data.frame(p_value)

#Possible combination add p values to gene names 
row_names <- row.names(k_means_data)
p_combined <- pvalues
rownames(p_combined) <- row_names

#this is just for our markdown where we have to show our clusters, we have to delete it in our code after we made the markdown
#View(beta_genes_healthy)
#cluster <- data.frame(k[["cluster"]])
#View(cluster)

#correction of p values by using holm method
#Nun haben wir p values, korrigierte p values und Gen IDs nebeneinander. Aber sobald wir die nach gr????e sortieren mit der sort Funktion, werden die Gen IDs gel??scht und man k??nnte die Werte nicht mher den Genen zuordnen, weil sie ja auch nicht mehr in der richtigen Reihenfolge sind.Hmm..
p_combined$p_adjusted = p.adjust(p_combined$p_value, method = "holm")
#p_combined = sort(p_combined$p_adjusted, decreasing = F)
p_holm <- data.frame(p_combined)

#Creating new dataframe to see our most important genes by name and order it in the correct way
symbols <- genes[rownames(p_holm),]
holm <- cbind(p_holm, "Symbols" = symbols$symbol)
holm_new <- holm[order(holm$p_value),]

#fold change calculation
#log2 fold change (use normal and not log)
#Da wir negative Werte haben in den M values, nutzen wir nicht den log2, obwohl man das normalerweise bei einem foldchange macht.
beta_10000 <- beta_genes_cleaned[genes_top_10000,]
beta_10000_log <- log2(beta_10000)
control <- beta_10000_log[, 1:5]
tumor <- beta_10000_log[, 6:10]
control_mean <- apply(control, 1, mean)
tumor_mean <- apply(tumor, 1, mean)
fold = control_mean - tumor_mean

#volcano plot (differenz der methylierung gesund und healthy (als foldchange) gegen logarythmischen p value (aber hier nicht den korrigierten, weil es sonst so scheisse aussieht mit den ganzen werten die durch die korrektur zu 1 gesetzt wurden?) )
plot(fold, -log10(p_value), ylim = c(0,3))
threshold_significant = 0.1
threshold_fold = 1
#biologically relevant
filter_by_fold = abs(fold) >= threshold_fold
#statistically relevant
filter_by_pvalue = abs(p_value) <= threshold_significant
filter_combined = filter_by_fold & filter_by_pvalue
#upmethylated
plot(fold, -log10(p_value))
points (fold[filter_combined & fold < 0],
        -log10(p_value[filter_combined & fold < 0]),
        pch = 16, col = "red")
#hypomethylated
points (fold[filter_combined & fold > 0],
        -log10(p_value[filter_combined & fold > 0]),
        pch = 16, col = "blue")

#logistic regression
##We create a dataframe, which contains M values and healthstatus of 7 random patients. Therefor the M value dataframe needs to be transformed first.

log_regression <- k_means_data[,c(1,2,4,5,8,9,10)]
log_regression <- t(log_regression)
log_regression <- data.frame(log_regression)
test_set <- k_means_data[,c(3,6,7)]
test_set <- data.frame(t(test_set))


healthstatus <- annotation$DISEASE
healthstatus <- data.frame(healthstatus)
healthstatus_regression <- healthstatus[c(1, 2, 4, 5, 8, 9, 10),]
                                      

train_set <- cbind(healthstatus_regression, log_regression)
#regression_model <- glm(healthstatus_regression ~ train_set[,1], family = "binomial", data = log_regression)
regression_model <- glm(healthstatus_regression ~ ., family = "binomial", data = train_set)
summary(regression_model)
prediction <- predict(regression_model, newdata = test_set, type = "response")
levels(train_set$healthstatus_regression)

                                      
#glm logistic regression (multicollinearity) without cross validation
k_means_data_no_cv <- t(k_means_data)
k_means_data_no_cv <- data.frame(k_means_data_no_cv)
full_data <- cbind(healthstatus, k_means_data_no_cv)
model_no_cv <- glm(healthstatus ~ ., family = "binomial", data = full_data)
prediction_no_cv <- predict(model_no_cv, newdata = full_data, type = "response")
summary(model_no_cv)
