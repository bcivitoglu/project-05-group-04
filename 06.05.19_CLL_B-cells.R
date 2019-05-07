#load data
input.data <- readRDS(file = "CLL-Bcells_list.RDS")

annotation <- read.csv("sample_annotation.csv")


#divide data into subgroups

genes <- input.data$genes

promoters <- input.data$promoters

cpgislands <- imput.data$cpgislands


#rename columns containing patient data

colnames(genes)[11] <- "Patient1_healthy_beta"
colnames(genes)[12] <- "Patient2_healthy_beta"
colnames(genes)[13] <- "Patient3_healthy_beta"
colnames(genes)[14] <- "Patient4_healthy_beta"
colnames(genes)[15] <- "Patient5_healthy_beta"
colnames(genes)[16] <- "Patient6_CLL_beta"
colnames(genes)[17] <- "Patient7_CLL_beta"
colnames(genes)[18] <- "Patient8_CLL_beta"
colnames(genes)[19] <- "Patient9_CLL_beta"
colnames(genes)[20] <- "Patient10_CLL_beta"
colnames(genes)[21] <- "Patient1_healthy_coverage"
colnames(genes)[22] <- "Patient2_healthy_coverage"
colnames(genes)[23] <- "Patient3_healthy_coverage"
colnames(genes)[24] <- "Patient4_healthy_coverage"
colnames(genes)[25] <- "Patient5_healthy_coverage"
colnames(genes)[26] <- "Patient6_CLL_coverage"
colnames(genes)[27] <- "Patient7_CLL_coverage"
colnames(genes)[28] <- "Patient8_CLL_coverage"
colnames(genes)[29] <- "Patient9_CLL_coverage"
colnames(genes)[30] <- "Patient10_CLL_coverage"