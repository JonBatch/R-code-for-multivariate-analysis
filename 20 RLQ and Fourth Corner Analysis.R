install.packages('ade4')
install.packages('vegan')
install.packages('dummies')

library('vegan')
library('dummies')
library('ade4')
source('biostats.R')

speabu<-read.csv('MAHA_speciesabu.csv',header = TRUE, row.names = 1)
envdata<-read.csv('MAHA_environment.csv',header=TRUE, row.names = 1)
spetra<-read.csv('MAHA_speciestrait.csv', header=TRUE, row.names = 1)
str(spetra)

spetra.tran <- dummy.data.frame(spetra)
spetra.tran

spetra.tran<-data.stand(spetra.tran,method='range',margin='column')
speabu.tran<-data.stand(speabu,method ='total',margin = 'row', plot=F)
envdata.tran<-data.trans(envdata,method = 'log',plot=F)

L.species<-dudi.coa(speabu.tran,scannf=FALSE)
Q.trait<-dudi.pca(spetra.tran,row.w = L.species$cw,scannf=FALSE)
R.env<-dudi.pca(envdata.tran,row.w=L.species$lw,scannf =FALSE)
rlq.MAHA<-rlq(R.env, L.species, Q.trait,scannf = FALSE)
