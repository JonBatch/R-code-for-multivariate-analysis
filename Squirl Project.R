install.packages('vegan')
install.packages('dplyr')
install.packages('factoextra')
install.packages('FactoMineR')
install.packages('cluster')
library('vegan')
library('dplyr')
library('FactoMineR')
library('factoextra')
library('cluster')

########Coralation between trapping years##########

SPECIESCOUNT <- read_csv("SPECIESCOUNT.csv")
cor.test(SPECIESCOUNT$'GLSA2016', SPECIESCOUNT$'GLSA2015', method = "spearman")
cor.test(SPECIESCOUNT$'GLSA2014', SPECIESCOUNT$'GLSA2015', method = "spearman")
cor.test(SPECIESCOUNT$'GLSA2016', SPECIESCOUNT$'GLSA2014', method = "spearman")
cor.test(SPECIESCOUNT$'NEFU2016', SPECIESCOUNT$'NEFU2015', method = "spearman")
cor.test(SPECIESCOUNT$'NEFU2014', SPECIESCOUNT$'NEFU2015', method = "spearman")
cor.test(SPECIESCOUNT$'NEFU2016', SPECIESCOUNT$'NEFU2014', method = "spearman")
cor.test(SPECIESCOUNT$'NECI2016', SPECIESCOUNT$'NECI2015', method = "spearman")
cor.test(SPECIESCOUNT$'NECI2014', SPECIESCOUNT$'NECI2015', method = "spearman")
cor.test(SPECIESCOUNT$'NECI2016', SPECIESCOUNT$'NECI2014', method = "spearman")
cor(SPECIESCOUNT, use = "pairwise.complete.obs", method = "spearman")

#################NMDS to identify the highest loadings#################

PlotOpen<-read.csv('Squirl Open.csv', header = TRUE, row.names = 1) 
PlotDepth<-read.csv('Squirl Depth.csv', header=TRUE, row.names = 1)
Groups<-read.csv('Squirl groups.csv', header=TRUE, row.names = 1)

?daisy #vegdist can't use gower to deal with missing values
Ddist<-daisy(PlotDepth, 'gower')
Odist<-daisy(PlotOpen,'gower')

  
DepthNMDS<-metaMDS(Ddist, k=2, autotransform = FALSE, trymax = 200)
DepthNMDS
nmds.scree(Ddist,k=10, autotransform = FALSE, trymax = 30)
nmds.monte(Ddist, k=3, autotransform = FALSE, trymax = 20)
stressplot(DepthNMDS)
plot(DepthNMDS)

vec.D<-envfit(DepthNMDS$points,PlotDepth,perm=2000, na.rm = TRUE)
vec.D #D38-D47 have the highest loadings on NMDS1 (<0.97) D42 is highest at 0.99974
ordiplot(DepthNMDS, choices = c(1,2), display = "sites")
plot(vec.D, p.max = 0.01, col = "blue")

#subset for D38-D47 
DepthSub <- PlotDepth[,c(54:63)]
#NMDS on Subset and switching from gowers to Euclidian
DepthNMDSub<-metaMDS(DepthSub, distance = "euc", k=2, autotransform = FALSE, trymax = 200)


OpenNMDS<-metaMDS(Odist, k=2, autotransform = FALSE, trymax = 200)
OpenNMDS
nmds.scree(Odist,k=10, autotransform = FALSE, trymax = 30)
nmds.monte(Odist, k=3, autotransform = FALSE, trymax = 20)
stressplot(OpenNMDS)
plot(OpenNMDS)

vec.O<-envfit(OpenNMDS$points,PlotOpen,perm=2000)
vec.O$vectors #O38-O47 have the highest loadings on NMDS1 (<0.97) O42 is highest at 0.99975
ordiplot(OpenNMDS, choices = c(1,2), display = "sites")
plot(vec.O, p.max = 0.01, col = "blue")

#subset for O38-O47
OpenSub <- PlotOpen[,c(54:63)]
#NMDS on Subset and switching from gowers to Euclidian
OpenNMDSub<-metaMDS(OpenSub, distance = "euc", k=2, autotransform = FALSE, trymax = 200)

###############PROCRUSTIES###################
protest(OpenNMDS,DepthNMDS)
crust<-procrustes(OpenNMDS,DepthNMDS)
plot(crust,kind=1,ar.col="red",len=0.2)
plot(crust,kind=2)

protest(OpenNMDSub,DepthNMDSub)
crust<-procrustes(OpenNMDSub,DepthNMDSub)
plot(crust,kind=1,ar.col="red",len=0.2)
plot(crust,kind=2)


########Mantel#######
OMAN<-vegdist(OpenSub,'euc')
DMAN<-vegdist(DepthSub,'euc')
se.man<-mantel(OMAN,DMAN,method='pearson')
se.man

OMAN<-daisy(PlotOpen,'gower')
DMAN<-daisy(PlotDepth,'gower')
se.man<-mantel(OMAN,DMAN,method='pearson')
se.man

#########ANOSIM#########

?vegdist
Opendist<-vegdist(OpenSub,"euc")
Depthdist<-vegdist(DepthSub,"euc")

O.anosim<-anosim(Opendist,Groups[,2])
summary(O.anosim)
plot.anosim(O.anosim)


D.anosim<-anosim(Depthdist,Groups[,2])
summary(D.anosim)
plot.anosim(D.anosim)




data<-cbind(OpenSub,Groups)
newdata<-data[which(data$GLSA.RANK=='VHIGH' | data$GLSA.RANK=='HIGH'),]
newdata$GLSA.RANK<-factor(newdata$GLSA.RANK)
OpenDD<-vegdist(newdata[,1:10],"euc")
OO.anosim<-anosim(OpenDD,newdata$GLSA.RANK)
summary(OO.anosim)

data<-cbind(OpenSub,Groups)
newdata<-data[which(data$GLSA.RANK=='VHIGH' | data$GLSA.RANK=='MID'),]
newdata$GLSA.RANK<-factor(newdata$GLSA.RANK)
OpenDD<-vegdist(newdata[,1:10],"euc")
OO.anosim<-anosim(OpenDD,newdata$GLSA.RANK)
summary(OO.anosim)

data<-cbind(OpenSub,Groups)
newdata<-data[which(data$GLSA.RANK=='VHIGH' | data$GLSA.RANK=='LOW'),]
newdata$GLSA.RANK<-factor(newdata$GLSA.RANK)
OpenDD<-vegdist(newdata[,1:20],"euc")
OO.anosim<-anosim(OpenDD,newdata$GLSA.RANK)
summary(OO.anosim)

data<-cbind(OpenSub,Groups)
newdata<-data[which(data$GLSA.RANK=='HIGH' | data$GLSA.RANK=='MID'),]
newdata$GLSA.RANK<-factor(newdata$GLSA.RANK)
OpenDD<-vegdist(newdata[,1:20],"euc")
OO.anosim<-anosim(OpenDD,newdata$GLSA.RANK)
summary(OO.anosim)

data<-cbind(OpenSub,Groups)
newdata<-data[which(data$GLSA.RANK=='HIGH' | data$GLSA.RANK=='LOW'),]
newdata$GLSA.RANK<-factor(newdata$GLSA.RANK)
OpenDD<-vegdist(newdata[,1:20],"euc")
OO.anosim<-anosim(OpenDD,newdata$GLSA.RANK)
summary(OO.anosim)

data<-cbind(OpenSub,Groups)
newdata<-data[which(data$GLSA.RANK=='MID' | data$GLSA.RANK=='LOW'),]
newdata$GLSA.RANK<-factor(newdata$GLSA.RANK)
OpenDD<-vegdist(newdata[,1:20],"euc")
OO.anosim<-anosim(OpenDD,newdata$GLSA.RANK)
summary(OO.anosim)

#####Depth Anosim#########
data<-cbind(DepthSub,Groups)
newdata<-data[which(data$GLSA.RANK=='VHIGH' | data$GLSA.RANK=='HIGH'),]
newdata$GLSA.RANK<-factor(newdata$GLSA.RANK)
DepthDD<-vegdist(newdata[,1:20],"euc")
OO.anosim<-anosim(DepthDD,newdata$GLSA.RANK)
summary(OO.anosim)

data<-cbind(DepthSub,Groups)
newdata<-data[which(data$GLSA.RANK=='VHIGH' | data$GLSA.RANK=='MID'),]
newdata$GLSA.RANK<-factor(newdata$GLSA.RANK)
DepthDD<-vegdist(newdata[,1:20],"euc")
OO.anosim<-anosim(DepthDD,newdata$GLSA.RANK)
summary(OO.anosim)

data<-cbind(DepthSub,Groups)
newdata<-data[which(data$GLSA.RANK=='VHIGH' | data$GLSA.RANK=='LOW'),]
newdata$GLSA.RANK<-factor(newdata$GLSA.RANK)
DepthDD<-vegdist(newdata[,1:20],"euc")
OO.anosim<-anosim(DepthDD,newdata$GLSA.RANK)
summary(OO.anosim)

data<-cbind(DepthSub,Groups)
newdata<-data[which(data$GLSA.RANK=='HIGH' | data$GLSA.RANK=='MID'),]
newdata$GLSA.RANK<-factor(newdata$GLSA.RANK)
DepthDD<-vegdist(newdata[,1:20],"euc")
OO.anosim<-anosim(DepthDD,newdata$GLSA.RANK)
summary(OO.anosim)

data<-cbind(DepthSub,Groups)
newdata<-data[which(data$GLSA.RANK=='HIGH' | data$GLSA.RANK=='LOW'),]
newdata$GLSA.RANK<-factor(newdata$GLSA.RANK)
DepthDD<-vegdist(newdata[,1:20],"euc")
OO.anosim<-anosim(DepthDD,newdata$GLSA.RANK)
summary(OO.anosim)

data<-cbind(DepthSub,Groups)
newdata<-data[which(data$GLSA.RANK=='MID' | data$GLSA.RANK=='LOW'),]
newdata$GLSA.RANK<-factor(newdata$GLSA.RANK)
DepthDD<-vegdist(newdata[,1:20],"euc")
OO.anosim<-anosim(DepthDD,newdata$GLSA.RANK)
summary(OO.anosim)