#Nested for loop

##Load dataset

input_data <- readRDS(file ="CLL-Bcells_list.RDS.gz")

##Add subgroup "genes"

genes <- input_data$genes

##Make a dataset with coverage values only

cov_genes <- genes[,c(21:30)]

##Mean of every gene (check for NA´s first and set them to zero if NA´s available)

sum(is.na(cov_genes))

cov_genes[is.na(cov_genes)] <- 0

cov_genes_mean <- rowMeans(data.matrix(cov_genes))

##Set threshhold

threshhold1 <- 10

threshhold2 <- quantile(cov_genes_mean, probs = seq(0.95,0.95,0.05))


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
