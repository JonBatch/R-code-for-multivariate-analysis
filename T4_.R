setwd("E:/R-code-for-multivariate-analysis")

install.packages("vegan")
install.packages("cluster")
install.packages("pvclust")

library("vegan")
library("cluster")
library("pvclust")
library("stats")


source("E:/R-code-for-multivariate-analysis/biostats.R")

envdata<- read.csv('MAHA_environment.csv', header=TRUE, row.names = 1)
speabu<- read.csv('MAHA_speciesabu.csv', header=TRUE, row.names=1)

envdata.tra<-data.stand(envdata, method='standardize', margin = 'column',plot=F)
mean(envdata.tra[,1])

site.eucd<-vegdist(envdata.tra, method = 'euclidean')
?hclust

sitecl.ave<-hclust(site.eucd,method="average")

