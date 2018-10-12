install.packages("vegan")
install.packages("pastecs")
source("biostats.R")
library(vegan)
library(pastecs)


envdata <- read.csv('MAHA_environment.csv', header=TRUE, row.names = 1)
spedata <- read.csv('MAHA_speciesabu.csv', header = TRUE,row.names=1)

str(envdata)
stat.desc(envdata)
stat.desc(spedata)
testdata <- replace.missing(spedata)
testdatamean <- replace.missing(spedata, method = "mean")
testdatadrop <- drop.var(spedata, pct.missing = 5)
foa.plots(spedata)

testdata <- drop.var(spedata,min.po = 5)
testdata <- drop.var(spedata,min.fo = 5)
testdata <- drop.var(spedata,max.po=95)
testdata <- drop.var(spedata,min.cv = 5)

str(spedata)
str(testdata)

ecdf.plots(spedata)
hist.plots(spedata)
box.plots(spedata)
qqnorm.plots(spedata)
uv.plots(spedata)


data.trans(spedata,method = 'log')
data.trans(spedata.method='power', exp=.5)
