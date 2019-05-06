
#methods to rename columns

colnames(dataset)[colnames(dataset) == "old_name" ] <- "new_name"

names(dataset)[columnnumber] <- "new_name"

#rename our patients

names(genes)[11] <- "Patent1_healthy_beta"
names(genes)[12] <- "Patent2_healthy_beta"
names(genes)[13] <- "Patent3_healthy_beta"
names(genes)[14] <- "Patent4_healthy_beta"
names(genes)[15] <- "Patent5_healthy_beta"
names(genes)[16] <- "Patent6_CLL_beta"
names(genes)[17] <- "Patent7_CLL_beta"
names(genes)[18] <- "Patent8_CLL_beta"
names(genes)[19] <- "Patent9_CLL_beta"
names(genes)[20] <- "Patent10_CLL_beta"
names(genes)[21] <- "Patent1_healthy_coverage"
names(genes)[22] <- "Patent2_healthy_coverage"
names(genes)[23] <- "Patent3_healthy_coverage"
names(genes)[24] <- "Patent4_healthy_coverage"
names(genes)[25] <- "Patent5_healthy_coverage"
names(genes)[26] <- "Patent6_CLL_coverage"
names(genes)[27] <- "Patent7_CLL_coverage"
names(genes)[28] <- "Patent8_CLL_coverage"
names(genes)[29] <- "Patent9_CLL_coverage"
names(genes)[30] <- "Patent10_CLL_coverage"