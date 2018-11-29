
install.packages("vegan")
library("vegan")
install.packages("pastecs")
library("pastecs")

setwd("D:/GIT/R-code-for-multivariate-analysis")
source("biostats.R")

envdata<-read.csv('MAHA_environment.csv', header=TRUE, row.names = 1)
spedata<-read.csv('MAHA_speciesabu.csv', header = TRUE, row.names=1)

FSD <-read.csv('Squirl Depth.csv', header = TRUE, row.names = 1)
FSO <-read.csv('Squirl Open.csv', header = TRUE, row.names = 1)

str(envdata)
stat.desc(envdata)
stat.desc(spedata)

str(FSD)
stat.desc(FSD)
stat.desc(FSO)

testdata<- replace.missing(spedata)
testdata<- drop.var(spedata, pct.missing = 5)

foa.plots(spedata)
foa.plots(FSD)



testdata<-drop.var(spedata,min.fo=5)
str(spedata)
str(testdata)

testdata<-drop.var(spedata,max.po = 95)
testdata<-drop.var(spedata,min.cv = 5)

ecdf.plots(spedata)
hist.plots(spedata)
box.plots(spedata)
qqnorm.plots(spedata)
uv.plots(spedata)


data.trans(spedata,method='log')
data.trans(spedata,method = 'power',exp = .5)
data.trans(spedata,method='power',exp = 0)
data.trans(spedata,method = 'asin')


data.stand(spedata,method='total', margin = 'row')

uv.outliers(envdata,id='Sinuosity:BasinAre',var='Elev',sd.limit = 1)
mv.outliers(envdata,method='euclidean', sd.limit = 1)


###My DATA###

