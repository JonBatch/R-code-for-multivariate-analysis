
install.packages("vegan")
install.packages("MASS")
install.packages("caret")
installed.packages("e1071")
library("vegan")
library("e1071")
library("MASS")
library("caret")
source('D:/GIT/R-code-for-multivariate-analysis/biostats.R')

envdata<-read.csv('MAHA_environment.csv', header = TRUE, row.names = 1)
sitegroup<- read.csv("MAHA_groups.csv", header = TRUE, row.names = 1)
envdata.tran<-data.trans(envdata, method = 'log',plot = F)

env.lda<-lda(envdata.tran, sitegroup[,1])
env.lda
