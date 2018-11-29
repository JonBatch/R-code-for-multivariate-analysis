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



PlotOpen<-read.csv('Squirl Open.csv', header = TRUE, row.names = 1) 
PlotDepth<-read.csv('Squirl DepthTRUNK.csv',header=TRUE, row.names = 1)
PlotDepthNA<-read.csv('Squirl Depth.csv', header=TRUE, row.names = 1)
Groups<-read.csv('Squirl groups.csv', header=TRUE, row.names = 1)

OpenPCA<-prcomp(PlotOpen,scale. = TRUE)
summary(OpenPCA)
Open.eig<-OpenPCA$sdev^2
Open.eig
pca.eigenval(OpenPCA)

#can't use PCA on Open due to missing values
DepthPCA<-prcomp(PlotDepth,scale. = TRUE)
summary(DepthPCA)
Depth.eig<-DepthPCA$sdev^2
Depth.eig
pca.eigenval(DepthPCA)



screeplot(OpenPCA, bstick=TRUE)
ordi.monte(PlotOpen, ord='pca',dim = 5)
screeplot(DepthPCA, bstick=TRUE)
ordi.monte(PlotDepth, ord='pca',dim = 5)


OpenPCA$rotation
pca.eigenvec(OpenPCA, dim = 5,digits = 3,cutoff = 0.1)
pca.structure(OpenPCA,PlotOpen,dim = 5, cutoff = 0.4)
Openscores<-OpenPCA$x
Openscores


biplot(OpenPCA)
OpenPCA2<-PCA(PlotOpen,graph = F)
fviz_screeplot(OpenPCA2, addlabels = TRUE, ylim=c(0,100))
fviz_pca_var(OpenPCA2, col.var= "black")
fviz_pca_var(OpenPCA2, col.var= "contrib", gradient.cols= c("#00AFBB","#E7B800","#FC4E07"), repel = TRUE )
fviz_contrib(OpenPCA2,choice = "var",axes = 1, top = 20)
 fviz_pca(OpenPCA2, label="var", habillage=Groups$GLSA.RANK, palette=c("orange","green","yellow","red"), addEllipses=TRUE)

biplot(DepthPCA)
DepthPCA2<-PCA(PlotDepth,graph = F)
fviz_screeplot(DepthPCA2, addlabels = TRUE, ylim=c(0,100))
fviz_pca_var(DepthPCA2, col.var= "black")
fviz_pca_var(DepthPCA2, col.var= "contrib", gradient.cols= c("#00AFBB","#E7B800","#FC4E07"), repel = TRUE )
fviz_contrib(DepthPCA2,choice = "var",axes = 1, top = 20)
fviz_pca(DepthPCA2, label="var", habillage=Groups$GLSA.RANK, palette=c("orange","green","yellow","red"), addEllipses=TRUE)
?fviz
#subset for O47-O66 and new PCA
OpenSub <- PlotOpen[,c(51:70)]
OpenPCA3<-PCA(OpenSub,graph = F)
fviz_pca(OpenPCA3, label="var", habillage=Groups$GLSA.RANK, palette=c("orange","green","yellow","red"), addEllipses=TRUE)
OpenSubPCA<-prcomp(OpenSub, scale. = TRUE)
OpenLocation<-OpenSubPCA$x
OpenLocation

#subset for D37-D56 and new PCA
DepthSub <- PlotDepth[,c(39:58)]
DepthPCA3<-PCA(DepthSub,graph = F)
fviz_pca(DepthPCA3, label="var", habillage=Groups$GLSA.RANK, palette=c("orange","green","yellow","red"), addEllipses=TRUE)
DepthSubPCA<-prcomp(DepthSub, scale. = TRUE)
DepthLocation<-DepthSubPCA$x
DepthLocation

#subset for D41-D45 and new PCA
Depthopt <- PlotDepth[,c(49:54)]


#########Ploting PCAlocations#########
#boxplots
?boxplot
BoxOpen<-cbind(OpenLocation,Groups)
o <- ordered(BoxOpen$names, levels = c("Bonanza", "Chintiminy", "Three", "Wildcat", "Erikson", "Easy", "Bull", "Farmer", "Buzzard", "Schooner", "Beehave", "Savage", "Trail", "East", "Potnu", "Ferngully"))
boxplot(PC1~o, data=BoxOpen,notch=FALSE)
boxplot(MDS1~o, data=BoxOpen,notch=FALSE)


BoxDepth<-cbind(DepthLocation,Groups)
o <- ordered(BoxOpen$names, levels = c("Bonanza", "Chintiminy", "Three", "Wildcat", "Erikson", "Easy", "Bull", "Farmer", "Buzzard", "Schooner", "Beehave", "Savage", "Trail", "East", "Potnu", "Ferngully"))
boxplot(PC1~o, data=BoxDepth,notch=FALSE)




write.csv(BoxOpen, file = "BoxOpen.csv")
 
#########Moving into PCoA################

?vegdist
Opendist<-vegdist(PlotOpen,"altGower")
Depthdist<-vegdist(PlotDepth,"altGower")

OpenPCOA<-cmdscale(Opendist,eig = TRUE,add = T)
plot(OpenPCOA$eig[1:35]/sum(OpenPCOA$eig)*100,type="b",lwd=2,col="blue",xlab=
       "Principal Component from PCoA", ylab="% variation explained", main="%
     variation explained by PCoA (blue) vs. random expectation (red)")
lines(bstick(35)*100,type = "b",lwd=2,col="red")

OpenPCOA<-cmdscale(Opendist,eig = TRUE,add = T)
plot(OpenPCOA$eig[1:35]/sum(OpenPCOA$eig)*100,type="b",lwd=2,col="blue",xlab=
       "Principal Component from PCoA", ylab="% variation explained", main="%
     variation explained by PCoA (blue) vs. random expectation (red)")
lines(bstick(35)*100,type = "b",lwd=2,col="red")

ordiplot(OpenPCOA, choices = c(1, 2), type="text", display="sites", xlab="PC
         1 (15%)", ylab="PC 2 (11%)")
vec.sp<-envfit(scores(OpenPCOA), PlotOpen, perm=1000)
vec.sp


#################Depth NMDS#################
DepthNMDS<-metaMDS(PlotDepth, distance = 'bray', k=2, autotransform = FALSE, trymax = 200)
DepthNMDSub<-metaMDS(DepthSub, distance = 'bray', k=2, autotransform = FALSE, trymax = 200)
DepthNMDopt<-metaMDS(Depthopt, distance = 'bray', k=2, autotransform = FALSE, trymax = 200)
DepthNMDS
nmds.scree(PlotDepth,distance = 'bray',k=10, autotransform = FALSE, trymax = 30)
nmds.monte(PlotDepth, distance = 'bray', k=3, autotransform = FALSE, trymax = 20)
stressplot(DepthNMDS)
stressplot(DepthNMDSub)
plot(DepthNMDS)
plot(DepthNMDopt)
?metaMDS


plot(DepthNMDS, type = 'n')
text(DepthNMDS, labels =Squirl_groups$`GLSA RANK`)
ordiellipse(DepthNMDS, Squirl_groups$`GLSA RANK`, show.groups = "VHIGH", col = c('red'))
ordiellipse(DepthNMDS, Squirl_groups$`GLSA RANK`, show.groups = "HIGH", col = c('Purple'))
ordiellipse(DepthNMDS, Squirl_groups$`GLSA RANK`, show.groups = "MID", col = c('Blue'))
ordiellipse(DepthNMDS, Squirl_groups$`GLSA RANK`, show.groups = "LOW", col = c('Green'))

plot(DepthNMDSub, type="n", display = "sites")
text(DepthNMDSub, labels =Squirl_groups$`GLSA RANK`)
ordiellipse(DepthNMDSub, Squirl_groups$`GLSA RANK`, show.groups = "VHIGH", col = c('red'))
ordiellipse(DepthNMDSub, Squirl_groups$`GLSA RANK`, show.groups = "HIGH", col = c('Purple'))
ordiellipse(DepthNMDSub, Squirl_groups$`GLSA RANK`, show.groups = "MID", col = c('Blue'))
ordiellipse(DepthNMDSub, Squirl_groups$`GLSA RANK`, show.groups = "LOW", col = c('Green'))

vec.D<-envfit(DepthNMDS$points,PlotDepth,perm=2000)
vec.Dsub<-envfit(DepthNMDSub$points,DepthSub,perm=2000)
vec.D
vec.Dsub

ordiplot(DepthNMDS, choices = c(1,2), display = "sites")
plot(vec.D, p.max = 0.01, col = "blue")
ordiplot(DepthNMDSub, choices = c(1,2), display = "sites")
plot(vec.Dsub, p.max = 0.01, col = "blue")

#funtion to plot points and elipses
NMDS = data.frame(MDS1 = DepthNMDSub$points[,1], MDS2 = DepthNMDSub$points[,2],group=Squirl_groups$`GLSA RANK`)
ord<-ordiellipse(DepthNMDSub, Squirl_groups$`GLSA RANK`, display = "sites", 
                 kind = "se", conf = 0.95, label = T)
df_ell <- data.frame()
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}
for(g in levels(NMDS$group)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$group==g,],
                                                   veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale)))
                                ,group=g))
}
ggplot(data = NMDS, aes(MDS1, MDS2)) + geom_point(aes(color = group), size=3) +
  geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2,colour=group), size=2, linetype=1)

#########Ploting DEPTH NMDS Point locations#########
DepthNMDSub$points


BoxDepth<-cbind(DepthLocation,Groups)
o <- ordered(BoxOpen$names, levels = c("Bonanza", "Chintiminy", "Three", "Wildcat", "Erikson", "Easy", "Bull", "Farmer", "Buzzard", "Schooner", "Beehave", "Savage", "Trail", "East", "Potnu", "Ferngully"))
gr <- ordered(BoxOpen$GLSA.RANK, levels = c("LOW", "MID", "HIGH", "VHIGH"))
nr<-  ordered(BoxOpen$NECI.RANK, levels = c("LOW", "MID", "HIGH", "VHIGH"))
boxplot(DMDS1~o, data=BoxOpen,notch=FALSE)
boxplot(DMDS1~gr, data = BoxOpen)
boxplot(D45~o, data = DepthSub)
boxplot(D45~gr, data = DepthSub)
plot(Groups$GLSA, DepthSub$D45)
plot(Groups$NECI, DepthSub$D45)



#################OPEN NMDS#################
OpenNMDS<-metaMDS(PlotOpen, distance = 'bray', k=2, autotransform = FALSE, trymax = 200)
OpenNMDSub<-metaMDS(OpenSub, distance = 'bray', k=2, autotransform = FALSE, trymax = 200)
OpenNMDS
nmds.scree(PlotOpen,distance = 'bray',k=10, autotransform = FALSE, trymax = 30)
nmds.monte(PlotOpen, distance = 'bray', k=3, autotransform = FALSE, trymax = 20)
stressplot(OpenNMDS)
stressplot(OpenNMDSub)
plot(OpenNMDS)
plot(OpenNMDSub)


plot(OpenNMDS, type = 'n')
text(OpenNMDS, labels =Squirl_groups$`GLSA RANK`)
ordiellipse(OpenNMDS, Squirl_groups$`GLSA RANK`, show.groups = "VHIGH", col = c('red'))
ordiellipse(OpenNMDS, Squirl_groups$`GLSA RANK`, show.groups = "HIGH", col = c('Purple'))
ordiellipse(OpenNMDS, Squirl_groups$`GLSA RANK`, show.groups = "MID", col = c('Blue'))
ordiellipse(OpenNMDS, Squirl_groups$`GLSA RANK`, show.groups = "LOW", col = c('Green'))

plot(OpenNMDSub, type = 'n')
text(OpenNMDSub, labels =Squirl_groups$`GLSA RANK`)
ordiellipse(OpenNMDSub, Squirl_groups$`GLSA RANK`, show.groups = "VHIGH", col = c('red'))
ordiellipse(OpenNMDSub, Squirl_groups$`GLSA RANK`, show.groups = "HIGH", col = c('Purple'))
ordiellipse(OpenNMDSub, Squirl_groups$`GLSA RANK`, show.groups = "MID", col = c('Blue'))
ordiellipse(OpenNMDSub, Squirl_groups$`GLSA RANK`, show.groups = "LOW", col = c('Green'))

vec.O<-envfit(OpenNMDS$points,PlotOpen,perm=2000)
vec.Osub<-envfit(OpenNMDSub$points,OpenSub,perm=2000)
vec.O
vec.Osub
ordiplot(OpenNMDS, choices = c(1,2), display = "sites")
plot(vec.O, p.max = 0.01, col = "blue")
OpenNMDSub$points

#funtion to plot points and elipses
NMDS = data.frame(MDS1 = OpenNMDSub$points[,1], MDS2 = OpenNMDSub$points[,2],group=Squirl_groups$`GLSA RANK`)
ord<-ordiellipse(OpenNMDSub, Squirl_groups$`GLSA RANK`, display = "sites", 
                 kind = "se", conf = 0.95, label = T)
df_ell <- data.frame()
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}
for(g in levels(NMDS$group)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$group==g,],
                                                   veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale)))
                                ,group=g))
}
ggplot(data = NMDS, aes(MDS1, MDS2)) + geom_point(aes(color = group), size=3) +
  geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2,colour=group), size=2, linetype=1)

#########Ploting OPEN NMDS Point locations#########
OpenNMDSub$points

BoxOpen<-cbind(OpenLocation,Groups)
boxplot(PC1~o, data=BoxOpen,notch=FALSE)
boxplot(OMDS1~o, data=BoxOpen,notch=FALSE)


o <- ordered(BoxOpen$names, levels = c("Bonanza", "Chintiminy", "Three", "Wildcat", "Erikson", "Easy", "Bull", "Farmer", "Buzzard", "Schooner", "Beehave", "Savage", "Trail", "East", "Potnu", "Ferngully"))
gr <- ordered(BoxOpen$GLSA.RANK, levels = c("LOW", "MID", "HIGH", "VHIGH"))
nr<-  ordered(BoxOpen$NECI.RANK, levels = c("LOW", "MID", "HIGH", "VHIGH"))
boxplot(OMDS1~o, data=BoxOpen,notch=FALSE)
boxplot(OMDS1~gr, data = BoxOpen)
boxplot(O38~o, data = OpenSub)
boxplot(O38~gr, data = OpenSub)
plot(Groups$GLSA, OpenSub$O38)
plot(Groups$NECI, OpenSub$O38)



###############PROCRUSTIES###################
protest(OpenNMDS,DepthNMDS)
crust<-procrustes(OpenNMDS,DepthNMDS)
plot(crust,kind=1,ar.col="red",len=0.2)
plot(crust,kind=2)

########Mantel#######
OMAN<-vegdist(PlotOpen,'euc')
DMAN<-daisy(PlotDepthNA,'gower')

se.man<-mantel(OMAN,DMAN,method='pearson')
se.man


 


