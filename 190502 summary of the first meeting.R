#how to load our data and create new variables

annotation <- read.csv("sample_annotation.csv")

input_data <- readRDS(file = "CLL-Bcells_list.RDS.gz")

#how to get an overview of the data

View(annotation)

View(input_data)

#how to check your current working directory

getwd()

#Our input_data can be divided into 4 subgroups (tiling, genes, promotors, cpgislands).
#In our project we are going to work only with 3 of them so we have to create new variables.

genes <- input_data$genes

promotors <- input_data$promoters

cpgislands <- input_data$cpgislands


#how to see the dimension of the data (number of columns and rows)

dim()

