library(plotrix)
require(plotrix)
library(ggplot2)
mydata=read.table('data1')
boxplot(mydata, use.cols = TRUE, col =2:8, xlab="", names=c("CMAP","GPCP","TRMM","ERAI","MERRA","TMEAN"),pch=0, cex.main=2,cex.lab=1, frame.plot=TRUE)
title(ylab=(expression(mm ~ day^{-1})),cex.lab=1.25,line=2)
mydata1=read.table('cammodelsdat')
points(1:6,mydata1[1,],pch=1)
points(1:6,mydata1[2,],pch=2)
points(1:6,mydata1[3,],pch=4)
