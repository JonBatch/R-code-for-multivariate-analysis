install.packages("vegan")
library("vegan")
install.packages("pastecs")
library("pastecs")

setwd("D:/GIT/R-code-for-multivariate-analysis")
source("biostats.R")

envdata<-read.csv('MAHA_environment.csv', header=TRUE, row.names = 1)
spedata<-read.csv('MAHA_speciesabu.csv', header = TRUE, row.names=1)

str(envdata)
stat.desc(envdata)
stat.desc(spedata)

testdata<- replace.missing(spedata)
testdata<- drop.var(spedata, pct.missing = 5)

foa.plots(spedata)

testdata<-drop.var(spedata,min.fo=5)
str(spedata)
str(testdata)

testdata<-drop.var(spedata,max.po = 95)
testdata<-drop.var(spedata,min.cv = 5)

