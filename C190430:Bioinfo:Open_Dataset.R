#Create variable that contains our DNA-methylation dataset.
input_data <- readRDS(file = "CLL-Bcells_list.RDS.gz")

#We want to have a look at the data this set contains.
view(input_data)

#Create a variable that contains the sample annotation
annotation <- read.csv("sample_annotation.csv")

#To make sure which file on my computer offers the environment for my R project:
getwd()

#Another way to find the right the right datasat (to upload it in your R project) is to write the whole path:
annotation <- read.csv('/Users/carlottabrueggen/Desktop/MoBi/4. FS/Bioinfo/CLLvsBcells/sample_annotation.csv')

#Open the dataset annotation
View(annotation)

#For the start it's sufficient to have a look at just one type of DNA sequence, which is in our case the "genes".
#Later, when the code is written, you can use it for all type of DNA sequence, which is besided "genes" also "promoters" and "CpGislands".
genes <- input_data$genes

#This already looks much neatlier arranged:
View(genes)

#To get in averview at how big our dataset we work on is:
dim(genes)

[1] 56235    30


