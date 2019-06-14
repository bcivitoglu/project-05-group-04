#Nested for loop

##Load dataset

input_data <- readRDS(file ="CLL-Bcells_list.RDS.gz")

##Add subgroup "genes"

genes <- input_data$genes

##Make a dataset with coverage values only

cov_genes <- genes[,c(1,21:30)]

##Remove ChrX and ChrY

cov_genes_new <- cov_genes[-which(cov_genes =="chrX"),]
cov_genes <- cov_genes_new[-which(cov_genes_new =="chrY"),]
rm(cov_genes_new)

##Mean of every gene (check for NAÂ´s first and set them to zero if NAÂ´s available)

sum(is.na(cov_genes))

cov_genes[is.na(cov_genes)] <- 0

cov_genes_mean <- rowMeans(data.matrix(cov_genes))

##Set threshhold

threshhold1 <- 10

threshhold2 <- quantile(cov_genes_mean, probs = seq(0.95,0.95,0.05))

##Forloop umgehen? klappt nicht

cov_genes[cov_genes<=threshhold1] <- NA
cov_genes[cov_genes>=threshhold2] <- NA

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
rm(i,j,threshhold1,threshhold2)

##Dataframe for beta values(+Chromosomes)

beta_genes <- genes[,c(1,11:20)]

##Also remove chrX and chrY

beta_genes_new <- beta_genes[-which(beta_genes =="chrX"),]
beta_genes <- beta_genes_new[-which(beta_genes_new =="chrY"),]
rm(beta_genes_new)


##Nested for loop to bring Coverage NA´s to beta NA´s

for(k in 1:ncol(beta_genes)){
  for(l in 1:nrow(beta_genes)){
    if(isTRUE(is.na(cov_genes[k,l] == TRUE))){
      beta_genes[k,l] <- NA
    } 
  }
}
rm(k,l)

##Count NA´s in dataframe

rmv.rows_beta_genes = apply(beta_genes,1, function(x){sum(is.na(x))})

##New dataframe out of beta_genes, where all rows with more than 2 NA´s are removed

beta_genes_cleaned <- beta_genes[-which(rmv.rows_beta_genes >0),]


##Make a dataframe with cleaned beta values and the associated chromosome



