

setwd("D:/GIT/R-code-for-multivariate-analysis")
source("biostats.R")
install.packages("simba")
install.packages("vegan")
install.packages("ecodist")
install.packages("cluster")
library("simba")
library("vegan")
library("ecodist")
library("cluster")

envdata<-read.csv('MAHA_environment.csv', header = TRUE, row.names = 1)
speabu <- read.csv('MAHA_speciesabu.csv', header=TRUE, row.names = 1)
spetrait<-read.csv('MAHA_speciestrait.csv', header=TRUE, row.names=1)

speocc<-data.trans(speabu,method = 'power', exp=0, plot = F)

sp.jac<-sim(speocc,method="jaccard")
sp.jac
sp.sim<-sim(speocc,method = "simplematching")
sp.sor<-sim(speocc,method = "soerensen")
plot(sp.jac,sp.sim, xlab="Jaccard's coefficient", ylab="Simple Matching coefficient")
abline(0,1,col="darkgray")

sp.bray <- vegdist(speabu, method="bray")
sp.bray

install.packages("gclus")
library("gclus")
source(coldiss.R)

coldiss(sp.bray,nc=4, byrank=FALSE, diag=TRUE)

sp.jacd<-vegdist(speocc,method = "jaccard")
plot(sp.jac, 1-sp.jacd)
abline(0,1,col="darkgray")

speabu.norm<-decostand(speabu, method= 'nor')
sp.chordd<-dist(speabu.norm)

speabu.hel<-decostand(speabu, method='hel')
sp.held<-dist(speabu.hel)

str(spetrait)
sptr.gower <-daisy(spetrait,metric = "gower")
coldiss(sptr.gower,nc=4,byrank=FALSE,diag=TRUE)
env.euc<-vegdist(envdata, method="euclidean")
env.euc
