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

names(genes)[11] <- "P1_healthy_beta"
names(genes)[12] <- "P2_healthy_beta"
names(genes)[13] <- "P3_healthy_beta"
names(genes)[14] <- "P4_healthy_beta"
names(genes)[15] <- "P5_healthy_beta"
names(genes)[16] <- "P6_CLL_beta"
names(genes)[17] <- "P7_CLL_beta"
names(genes)[18] <- "P8_CLL_beta"
names(genes)[19] <- "P9_CLL_beta"
names(genes)[20] <- "P10_CLL_beta"
names(genes)[21] <- "P1_healthy_coverage"
names(genes)[22] <- "P2_healthy_coverage"
names(genes)[23] <- "P3_healthy_coverage"
names(genes)[24] <- "P4_healthy_coverage"
names(genes)[25] <- "P5_healthy_coverage"
names(genes)[26] <- "P6_CLL_coverage"
names(genes)[27] <- "P7_CLL_coverage"
names(genes)[28] <- "P8_CLL_coverage"
names(genes)[29] <- "P9_CLL_coverage"
names(genes)[30] <- "P10_CLL_coverage"

#Create a data frame containing coverage values

cov_genes <- genes[ ,c(21:30)]

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

#we set the lower threshold to the coverage value 10 
abline(v=log10(10))



#wie kann das sein wenn wir die 0 Werte schon gelöscht haben, dass wir dann noch Werte unter 0 haben???

#Nested for loop

##Mean of every gene (check for NA´s first and set them to zero if NA´s available)

sum(is.na(cov_genes))

cov_genes[is.na(cov_genes)] <- 0

cov_genes_mean <- rowMeans(data.matrix(cov_genes))

##Set threshhold

threshhold1 <- 10

threshhold2 <- quantile(cov_genes_mean, probs = seq(0.90,0.90,0.05))


##Nested for loop 

for(i in 1:ncol(cov_genes)){
  for(j in 1:nrow(cov_genes)){
    if(isTRUE(cov_genes[i,j] <= threshhold1)){
      cov_genes[i,j] <- NA
    } 
    if(isTRUE(cov_genes[i,j] >= threshhold2)){
      cov_genes[i,j] <- NA
    }
  }
}

#separate sick and healthy patient coverage for deleting genes with 3 or more coverage values outside the threshold later
cov_genes_healthy <- cov_genes[ ,c(1:5)]
cov_genes_cancer <- cov_genes[ ,c(6:10)]

#count the NAs in the 2 dataframes per row in a new data frame

rmv.rows_healthy_cov = apply(cov_genes_healthy,1, function(x){sum(is.na(x))})
rmv.rows_cancer_cov = apply(cov_genes_cancer,1, function(x){sum(is.na(x))})

# new coverage data frame where all rows with more than 2 NAs are removed

cov_genes_cancer_cleaned = cov_genes_cancer[-which(rmv.rows_cancer_cov > 2),]
cov_genes_healthy_cleaned = cov_genes_healthy[-which(rmv.rows_healthy_cov > 2),]

#divide the beta values into a dataframe of healthy and cancer cells

beta_genes_healthy <- genes[,c(11:15)]
beta_genes_cancer <- genes[,c(16:20)]

beta_genes_healthy = beta_genes_healthy[-which(rmv.rows_healthy_cov > 2),]
beta_genes_cancer = beta_genes_cancer[-which(rmv.rows_cancer_cov > 2),]

#count the NAs in the 2 dataframes per row in a new data frame

rmv.rows_healthy = apply(beta_genes_healthy,1, function(x){sum(is.na(x))})
rmv.rows_cancer = apply(beta_genes_cancer,1, function(x){sum(is.na(x))})

# new beta value data frame where all rows with more than 2 NAs are removed

beta_genes_healthy_cleaned = beta_genes_healthy[-which(rmv.rows_healthy > 2),]
beta_genes_cancer_cleaned = beta_genes_cancer[-which(rmv.rows_cancer > 2),]

Row_Difference = nrow(beta_genes_cancer)-nrow(beta_genes_cancer_cleaned)
View(Row_Difference)
genes_deleted_cancer_percentage = Row_Difference/nrow(beta_genes_cancer)*100
sum(genes_deleted_cancer_percentage)

Row_Difference_2 = nrow(beta_genes_healthy)-nrow(beta_genes_healthy_cleaned)
View(Row_Difference_2)
genes_deleted_healthy_percentage = Row_Difference_2/nrow(beta_genes_healthy)*100
sum(genes_deleted_healthy_percentage)

#Normalisation: Transform beta-values into M-values

M_genes_h <- log2(beta_genes_healthy/(1-(beta_genes_healthy)))
M_genes_c <- log2(beta_genes_cancer/(1-(beta_genes_cancer)))


