#create dataframe with only beta values

beta.values.genes <- genes[,c(11,12,13,14,15,16,17,18,19,20)]
View(beta.values.genes)

#find NAs in beta.values.genes
lapply(beta.values.genes, function (x) which(is.na(x)))

#create data frame/list? with all NAs
NA.beta.genes <- lapply(beta.values.genes, function (x) which(is.na(x)))
View(NA.beta.genes$Patient1_healthy_beta)

#count number of NAs per row
rowSums(is.na(beta.values.genes))

#do the same for only healthy patients
beta.healthy.genes <- genes[,c(11,12,13,14,15)]
rowSums(is.na(beta.healthy.genes))

beta.healthy.genes.counts <- rowSums(is.na(beta.healthy.genes))
View(beta.healthy.genes.counts)

#create data frame with 3 or more NAs so the corresponding genes can be removed later
remove.NAs <- (beta.healthy.genes.counts > 2)
View(remove.NAs)

#remove all beta values with 3 or more NAs
beta.healthy <- beta.healthy.genes.counts[-which((beta.healthy.genes.counts > 2) == TRUE) ]
View(beta.healthy)

