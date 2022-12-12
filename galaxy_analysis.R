# ####################################################### ####
# Expansion and evolution of the R programming language   ####
# Author:    Timothy L Staples                            ####
# ####################################################### ####

# Script purpose ####

# Conduct all statistical analyses and produce figures.

# This script makes extensive use of R-Studio's code folding syntax.
# Alt + x to collapse all folds on Windows/Linux: Cmd + Alt + x on Mac.

# SET DETAILS ####

rm(list=ls())
setwd("PATH TO THIS FILE")
sapply(paste0("./functions/", list.files("./functions")), source)

package.loader(c("parallel", "vegan", "mgcv", "lme4", "shape"))

# DATA PREP ####

funTableSub <- read.csv(unz("./outputs/commonFunctionLong.csv.zip",
                            "commonFunctionLong.csv"))

funTableSub$functionCat <- "other"
funTableSub$functionCat[funTableSub$basePacks] = "base"
funTableSub$functionCat[funTableSub$tidyverse] = "tidyverse"
funTableSub$functionCat[funTableSub$tidyExtend] = "tidyverse"
funTableSub$functionCat <- as.factor(funTableSub$functionCat)

funTableSub$urID <- paste0(funTableSub$userId, ".", funTableSub$repoId)

# # Proportionalize function counts
scriptN <- as.data.frame(tapply(funTableSub$count,funTableSub$ursID, sum, na.rm=TRUE))
scriptN$ursID <- rownames(scriptN)
colnames(scriptN)[1] = "scriptN"

funTableSub <- funTableSub[,colnames(funTableSub) != "scriptN"]

funTableSub <- merge(funTableSub, scriptN,
                     all.x=TRUE, all.y=FALSE, sort=FALSE)
funTableSub$prop <- funTableSub$count / funTableSub$scriptN

# GENERATE REPO-LEVEL RA DATA ####

# proportionalize means to a sum of 1 for each repo

funTableRepo <- do.call("rbind", lapply(split(funTableSub[,c("urID", "count", "function.")], 
                                              f=funTableSub$urID),
                                        function(x){
                                          
                                          # average RA across scripts
                                          data.frame(function.= sort(unique(x$function.)),
                                                                     count = tapply(x$count,
                                                                                   x$function.,
                                                                                   sum, na.rm=TRUE),
                                                     urID=x$urID[1])


                                        }))

# merge one brings across package data by function
subFun <- funTableSub[!duplicated(funTableSub$function.),
                      c("package", "function.", "functionCat")]
funTableRepo$package <- subFun$package[match(funTableRepo$function., subFun$function.)]
funTableRepo$functionCat <- subFun$functionCat[match(funTableRepo$function., subFun$function.)]

# merge two brings across time & user details for repo

funTableRepo <- merge(funTableRepo,
                      funTableSub[!duplicated(funTableSub$urID),
                                  c("urID", "userId", "monthsSinceJan10Created", "monthsSinceJan10")],
                      by.x="urID", by.y="urID", all.x=TRUE, all.y=FALSE, sort=FALSE)

funTableRepo$funAsCat <- as.factor(funTableRepo$function.)

funTableRepo$updateLag <- funTableRepo$monthsSinceJan10 - funTableRepo$monthsSinceJan10Created

# Prop data ####

repoCount <- as.data.frame(tapply(funTableRepo$count,
                                  funTableRepo$urID, sum))
repoCount$urID <- rownames(repoCount)
colnames(repoCount)[1] = "repoCount"

funTableRepo <- merge(x=funTableRepo, y=repoCount,
                      by.x="urID", by.y="urID", all.x=TRUE,
                      all.y=FALSE, sort=FALSE)

funTableRepo$prop <- funTableRepo$count / funTableRepo$repoCount

funTableSub$funAsCat <- as.factor(funTableSub$function.)

# Compositional change over time ####
#               nMDS ####

tempProp <- tapply(funTableRepo$prop,
                   list(funTableRepo$monthsSinceJan10Created,
                        funTableRepo$funAsCat),
                   sum, na.rm=TRUE)
tempProp[is.na(tempProp)] = 0
tempProp <- prop.table(tempProp, 1)

funOrd <- metaMDS(tempProp, distance="bray")
plot(funOrd, display="sites", type = "t")
points(funOrd$species, col="red", cex=0.5)

#               time RDA ####

monthVect <- as.numeric(rownames(funMonthProp))

funMonthSqrt <- sqrt(funMonthProp)
funMonthSqrt <- prop.table(funMonthSqrt, margin=1)

funrda <- dbrda(funMonthSqrt ~ monthVect, distance="bray")

#               Function changes over time ####

funOccur <- as.matrix(table(funTableSub$monthsSinceJan10, funTableSub$function.))
monthCount <- as.vector(table(funTableSub$monthsSinceJan10[!duplicated(funTableSub$ursID)]))
centreMonth <- max(as.numeric(rownames(funOccur)))
monthRange <- as.numeric(rownames(funOccur))

funTrends <- as.data.frame(t(sapply(1:ncol(funOccur), function(n){
  print(n)
  x <- funOccur[,n]
  xF <- monthCount - x
  monthVect <- monthRange - min(monthRange) - (centreMonth - min(monthRange))
  
  tempGlm <- glm(cbind(x, xF) ~ monthVect, family=binomial)
  return(as.vector(summary(tempGlm)$coefficients)[1:4])
  
})))
colnames(funTrends) = c("(Intercept)", "monthVect", "IntSE", "monthSE")
funTrends$function. <- colnames(funOccur)

funTop <- merge(funTrends, funTableSub[!duplicated(funTableSub$function.), 
                                       c("function.", "package", "tidyverse", "tidyExtend", "basePacks")],
                by.x="function.", by.y="function.", all.x=TRUE, all.y=FALSE, sort=FALSE)

funTop$cat = as.factor(paste0(funTop$tidyverse, ":", funTop$tidyExtend, ":", funTop$basePacks))
levels(funTop$cat) = c("other", "base", "tidyextend", "tidyverse")


# now do the same at the package level
packTableSub <- funTableSub[!duplicated(paste0(funTableSub$ursID,":", funTableSub$package)),]

packOccur <- as.matrix(table(packTableSub$monthsSinceJan10, packTableSub$package))
monthCount <- as.vector(table(packTableSub$monthsSinceJan10[!duplicated(packTableSub$ursID)]))

packTrends <- as.data.frame(t(sapply(1:ncol(packOccur), function(n){
  print(n)
  x <- packOccur[,n]
  xF <- monthCount - x
  monthVect <- monthRange - min(monthRange) - (centreMonth - min(monthRange))
  
  tempGlm <- glm(cbind(x, xF) ~ monthVect, family=binomial)
  return(as.vector(summary(tempGlm)$coefficients)[1:4])
  
})))
colnames(packTrends) = c("(Intercept)", "monthVect", "IntSE", "monthSE")
packTrends$package <- colnames(packOccur)

packTop <- merge(packTrends, funTableSub[!duplicated(funTableSub$package), 
                                         c("package", "tidyverse", "tidyExtend", "basePacks")],
                 by.x="package", by.y="package", all.x=TRUE, all.y=FALSE, sort=FALSE)

packTop$cat = as.factor(paste0(packTop$tidyverse, ":", packTop$tidyExtend, ":", packTop$basePacks))
levels(packTop$cat) = c("other", "base", "tidyextend", "tidyverse")

#               Plot ####

# vectors to include
library(viridisLite)
yearVect <- as.numeric(rownames(funOrd$points)) %/% 12
yearVect <- yearVect - min(yearVect) + 1
yearVect <- c(1, yearVect)
yearBase <- viridis(max(yearVect), option="C", end=0.8)
yearCol <- yearBase[yearVect]

yearPoints <- funOrd$points[as.numeric(rownames(funOrd$points)) %% 12 == 0, ]

library(shape)
pdf(date.wrap("./plots/2013-2020comp",".pdf"), height=12, width=4.5, useDingbats = FALSE)

split.screen(rbind(c(0.125,0.96,0.69,0.99),
                   c(0.125,0.96,0.37,0.65),
                   c(0.125,0.96,0.04,0.32),
                   c(0.125,0.96,0.37,0.65),
                   c(0.125,0.96,0.04,0.32)))

screen(1)
par(mar=c(0,0,0,0), ps=8, tcl=-0.25, las=1, mgp=c(3,0,0))
plot(funOrd, display="sites", xlim=c(-0.4,0.5), type="n", axes=FALSE)
axis(side=1)
axis(side=2, mgp=c(3,0.5,0))
mtext(side=1, line=0.75, text="nMDS 1")
mtext(side=2, line=1.75, text="nMDS 2", las=0)
box()

# convex hulls for 2 target years
baseCols <- col2rgb(yearBase)/255

sapply(3:(ncol(baseCols)+2), function(n){
  
  yearRows <- ((n*12)+1):((n+1)*12)
  
  if(sum(rownames(funOrd$points) %in% yearRows)==0){return(NULL)}
  sub <- funOrd$points[rownames(funOrd$points) %in% yearRows,]
  
  subHull <- sub[chull(sub),]
  polygon(subHull, border=NA, 
          col=rgb(baseCols[1,n-2], baseCols[2,n-2], baseCols[3,n-2], 0.35))
  
  text(x=colMeans(subHull)[1], colMeans(subHull)[2], labels = 2010 + (n),
       cex=1, col=yearBase[n-3], font=2)
  
})

sapply(2:nrow(funOrd$points), function(n){
  lines(funOrd$points[c(n-1,n),], lwd=1, col = yearCol[n])
})
points(funOrd$points, pch=21, cex=1, bg=yearCol)
# text(funOrd$points, labels=c("J","F","M","A","M","J","J","A","S","O","N","D"),
#      col="white")
text(x=relative.axis.point(0.02, "x"), y=relative.axis.point(0.975, "y"),
     adj=0, labels="(A)", font=2)

text(x=relative.axis.point(0.02, "x"), y=relative.axis.point(0.02, "y"),
     adj=0, labels=paste0("Stress = ", round(funOrd$stress, 3)))
text(x=relative.axis.point(0.875, "x"), y=relative.axis.point(0.02, "y"),
     adj=1, labels=expression("Variance explained by time = "))
text(x=relative.axis.point(0.98, "x"), y=relative.axis.point(0.02, "y"),
     adj=1, labels=paste0(sprintf("%.2f", round(summary(funrda)$cont[[1]][2,1]*100, 2)), "%"))
close.screen(1)

screen(2)
par(mar=c(0,0,0,0), ps=8, tcl=-0.25, las=1)

plot(plogis(funTop$`(Intercept)`) ~ funTop$monthVect, type="n", yaxs="i", xaxs="i",
     axes=FALSE, xlab="", ylab="", xlim=c(-0.1,0.175), ylim=log10(c(8e-7, 1)))
#plot(log10(mainTop$prop2020) ~ log10(mainTop$proprat), type="n", yaxs="i", xaxs="i",
#     axes=FALSE, xlab="", ylab="", xlim=c(-3,3), ylim=c(-6,-0.1))
mainXs <- par("usr")[1:2]
mainYs <- par("usr")[3:4]

axis(side=1, mgp=c(3,0,0))

axis(side=2, at=log10(c(1e-5, 0.0001,0.001,0.01,0.1,1)), 
     labels=c(1e-5, 0.0001,0.001,0.01,0.1,1),
     mgp=c(3,0.5,0))
axis(side=2, at=log10(c(seq(1e-6, 1e-5, 1e-6),
                        seq(1e-5,1e-4,1e-5),
                        seq(1e-4,1e-3,1e-4),
                        seq(1e-3, 1e-2, 1e-3),
                        seq(1e-2, 0.1, 1e-2),
                        seq(0.1,1,0.1))), tcl=-0.125, labels=NA)

mtext(side=2, line=2.15, 
      text="Probability of function occurrence in December 2021", las=0)
mtext(side=1, line=0.75, text="Slope of change in function occurrence probability over time")

abline(v=log10(1), lty="31", col="grey60")

funOut <- funTop[funTop$monthVect >= par("usr")[2],]
funIn <- funTop[funTop$monthVect < par("usr")[2],]

# base package points
baseFun <- funIn[funIn$basePacks,]
baseHull <- chull(baseFun[,c("(Intercept)", "monthVect")])
basecol <- col2rgb("#2369bd")/255

otherFun <- funIn[!funIn$tidyverse & !funIn$basePacks & !funIn$tidyExtend,]
otherHull <- chull(otherFun[, c("(Intercept)", "monthVect")])
othercol <- col2rgb("grey75")/255

tidyFun <- funIn[funIn$tidyverse | funIn$tidyExtend,]
tidyHull <- chull(tidyFun[, c("(Intercept)", "monthVect")])
tidycol <- col2rgb("#d12239")/255

# y-densities
densScale = 0.05

baseYDens <- density(log10(plogis(baseFun$`(Intercept)`)),
                     from=min(log10(plogis(baseFun$`(Intercept)`))),
                     to=max(log10(plogis(baseFun$`(Intercept)`))))
tidyYDens <- density(log10(plogis(tidyFun$`(Intercept)`)),
                     from=min(log10(plogis(tidyFun$`(Intercept)`))),
                     to=max(log10(plogis(tidyFun$`(Intercept)`))))
otherYDens <- density(log10(plogis(otherFun$`(Intercept)`)),
                      from=min(log10(plogis(otherFun$`(Intercept)`))),
                      to=max(log10(plogis(otherFun$`(Intercept)`))))

polygon(y = c(otherYDens$x[1], otherYDens$x, rev(otherYDens$x)[1]),
        x = c(par('usr')[1], par('usr')[1] + otherYDens$y * densScale, par('usr')[1]),
        border="grey75", col=rgb(othercol[1], othercol[2],othercol[3], 0.2))

polygon(y = c(tidyYDens$x[1], tidyYDens$x, rev(tidyYDens$x)[1]),
        x = c(par('usr')[1], par('usr')[1] + tidyYDens$y * densScale, par('usr')[1]),
        border="#d12239", col=rgb(tidycol[1], tidycol[2],tidycol[3], 0.2))

polygon(y = c(baseYDens$x[1], baseYDens$x, rev(baseYDens$x)[1]),
        x = c(par('usr')[1], par('usr')[1] + baseYDens$y * densScale, par('usr')[1]),
        border="#2369bd", col=rgb(basecol[1], basecol[2],basecol[3], 0.2))

# x-densities

densScale = 0.025

baseXDens <- density(baseFun$monthVect,
                     from=min(baseFun$monthVect),
                     to=max(baseFun$monthVect))
tidyXDens <- density(tidyFun$monthVect,
                     from=min(tidyFun$monthVect),
                     to=max(tidyFun$monthVect))
otherXDens <- density(otherFun$monthVect,
                      from=min(otherFun$monthVect),
                      to=max(otherFun$monthVect))

polygon(x = c(otherXDens$x[1], otherXDens$x, rev(otherXDens$x)[1]),
        y = c(par('usr')[3], par('usr')[3] + otherXDens$y * densScale, par('usr')[3]),
        border="grey75", 
        col=rgb(othercol[1], othercol[2],othercol[3], 0.2))

polygon(x = c(tidyXDens$x[1], tidyXDens$x, rev(tidyXDens$x)[1]),
        y = c(par('usr')[3], par('usr')[3] + tidyXDens$y * densScale, par('usr')[3]),
        border="#d12239", col=rgb(tidycol[1], tidycol[2],tidycol[3], 0.2))

polygon(x = c(baseXDens$x[1], baseXDens$x, rev(baseXDens$x)[1]),
        y = c(par('usr')[3], par('usr')[3] + baseXDens$y * densScale, par('usr')[3]),
        border="#2369bd", col=rgb(basecol[1], basecol[2],basecol[3], 0.2))

# Polygons

# with(otherFun[otherHull,],
#      polygon(y=log10(plogis(`(Intercept)`)), x=monthVect,
#              col=rgb(othercol[1], othercol[2], othercol[3], 0.2), border=NA))
# with(baseFun[baseHull,],
#      polygon(y=log10(plogis(`(Intercept)`)), x=monthVect,
#              col=rgb(basecol[1], basecol[2],basecol[3], 0.2), border=NA))
# with(tidyExtend[extendHull,],
#      polygon(y=log10(plogis(`(Intercept)`)), x=monthVect,
#              col=rgb(extendcol[1], extendcol[2],extendcol[3], 0.2), border="#d12239"))
# with(tidyFun[tidyHull,],
#      polygon(y=log10(plogis(`(Intercept)`)), x=monthVect,
#              col=rgb(tidycol[1], tidycol[2],tidycol[3], 0.2), border=NA))

with(otherFun,
     points(y=log10(plogis(`(Intercept)`)), x=monthVect, 
            pch=16, col="grey75", cex=0.5))
with(baseFun,
     points(y=log10(plogis(`(Intercept)`)), x=monthVect, 
            pch=21, bg="#2369bd", lwd=0.5, cex=0.5))
with(tidyFun,
     points(y=log10(plogis(`(Intercept)`)), x=monthVect, 
            pch=21, bg="#d12239", lwd=0.5, cex=0.5))

mtext(side=3, line=0.1, at=par("usr")[1], text="(B) Change in function use",
      adj=0, font=2)

legend(x=relative.axis.point(0.65, "x"), y=relative.axis.point(1.02, "y"),
       pch=c(21,21,21, 16), col=c("black","black", "grey75"), pt.bg=c("#2369bd", "#d12239", NA),
       legend = c("Base packages", "Tidyverse", "Other packages"), pt.lwd=0.5,
       bty="n", y.intersp=0.65, x.intersp=0.75)

# text(x=relative.axis.point(c(0.02), "x"), y=relative.axis.point(0.02, "y"),
#      labels=c("Used less over time"),
#      adj=0)
# 
# text(x=relative.axis.point(c(0.98), "x"),y=relative.axis.point(0.02, "y"),
#      labels="Used more over time",
#      adj=1)
# 
# Arrows(x0=relative.axis.point(0.3, "x"),
#        y0 = relative.axis.point(0.02, "y"),
#        y1 = relative.axis.point(0.02, "y"),
#        x1 = relative.axis.point(0.68, "x"),
#        arr.type="triangle", arr.width=0.1, arr.length=0.1)

box()
close.screen(2)

# FUNCTION TEXT
screen(4)
par(mar=c(0,0,0,0), ps=8, tcl=-0.25, las=1)

plot(plogis(funTop$`(Intercept)`) ~ funTop$monthVect, type="n", yaxs="i", xaxs="i",
     axes=FALSE, xlab="", ylab="", xlim=c(-0.1,0.175), ylim=log10(c(8e-7, 1)))

textFuns <- funTop[plogis(funTop$`(Intercept)`) > 0.05 |
                     abs(funTop$monthVect) > 0.02,]

# if we just want particular functions
#textFuns <- funTop[funTop$function. %in% c("merge", "left_join", "right_join"),]

text(x=textFuns$monthVect,
     y=log10(plogis(textFuns$`(Intercept)`)),
     labels=textFuns$function., cex=0.75,
     col=c("grey75", "#2369bd", "#d12239", "#d12239")[as.factor(textFuns$cat)])

with(funOut,
     Arrows(x0=relative.axis.point(0.98, "x"),
            x1=relative.axis.point(0.995, "x"),
            y0=log10(plogis(`(Intercept)`)),
            y1=log10(plogis(`(Intercept)`)),
            col=c("grey75", "#2369bd", "#d12239", "#d12239")[as.factor(cat)],
            arr.type="triangle", arr.length=0.1, arr.width=0.1))

with(funOut,
     text(x=relative.axis.point(0.98, "x"),
          y=log10(plogis(`(Intercept)`)),
          labels=function., cex=0.75, adj=1,
          col=c("grey75", "#2369bd", "#d12239", "#d12239")[as.factor(cat)]))

close.screen(4)

# package-levels
screen(3)
par(mar=c(0,0,0,0), ps=8, tcl=-0.25, las=1)

plot(plogis(packTop$`(Intercept)`) ~ packTop$monthVect, type="n", yaxs="i", xaxs="i",
     axes=FALSE, xlab="", ylab="", xlim=c(-0.05,0.11), ylim=log10(c(8e-5, 1.1)))
mainXs <- par("usr")[1:2]
mainYs <- par("usr")[3:4]

axis(side=1, mgp=c(3,0,0))

axis(side=2, at=log10(c(1e-5, 0.0001,0.001,0.01,0.1,1)), labels=c(1e-5, 0.0001,0.001,0.01,0.1,1),
     mgp=c(3,0.5,0))
axis(side=2, at=log10(c(seq(1e-6, 1e-5, 1e-6),
                        seq(1e-5,1e-4,1e-5),
                        seq(1e-4,1e-3,1e-4),
                        seq(1e-3, 1e-2, 1e-3),
                        seq(1e-2, 0.1, 1e-2),
                        seq(0.1,1,0.1))), tcl=-0.125, labels=NA)

mtext(side=2, line=2.15, 
      text="Probability of package occurrence in December 2021", las=0)
mtext(side=1, line=1, text="Slope of change in package occurrence probability over time")

abline(v=log10(1), lty="31", col="grey60")

funOut <- packTop[packTop$monthVect >= par("usr")[2],]
funIn <- packTop[packTop$monthVect < par("usr")[2],]

# base package points
baseFun <- funIn[funIn$basePacks,]
baseHull <- chull(cbind(log10(plogis(baseFun[,"(Intercept)"])), 
                        baseFun[,"monthVect"]))
basecol <- col2rgb("#2369bd")/255

otherFun <- funIn[!funIn$tidyverse & !funIn$basePacks & !funIn$tidyExtend,]
otherHull <- chull(cbind(log10(plogis(otherFun[,"(Intercept)"])), 
                         otherFun[,"monthVect"]))
othercol <- col2rgb("grey75")/255

tidyFun <- funIn[funIn$tidyverse | funIn$tidyExtend,]
tidyHull <- chull(cbind(log10(plogis(tidyFun[,"(Intercept)"])), 
                        tidyFun[,"monthVect"]))
tidycol <- col2rgb("#d12239")/255

# y-densities
densScale = 0.0125

baseYDens <- density(log10(plogis(baseFun$`(Intercept)`)),
                     from=min(log10(plogis(baseFun$`(Intercept)`))),
                     to=max(log10(plogis(baseFun$`(Intercept)`))))
tidyYDens <- density(log10(plogis(tidyFun$`(Intercept)`)),
                     from=min(log10(plogis(tidyFun$`(Intercept)`))),
                     to=max(log10(plogis(tidyFun$`(Intercept)`))))
otherYDens <- density(log10(plogis(otherFun$`(Intercept)`)),
                      from=min(log10(plogis(otherFun$`(Intercept)`))),
                      to=max(log10(plogis(otherFun$`(Intercept)`))))

polygon(y = c(otherYDens$x[1], otherYDens$x, rev(otherYDens$x)[1]),
        x = c(par('usr')[1], par('usr')[1] + otherYDens$y * densScale, par('usr')[1]),
        border="grey75", col=rgb(othercol[1], othercol[2],othercol[3], 0.2))

polygon(y = c(tidyYDens$x[1], tidyYDens$x, rev(tidyYDens$x)[1]),
        x = c(par('usr')[1], par('usr')[1] + tidyYDens$y * densScale, par('usr')[1]),
        border="#d12239", col=rgb(tidycol[1], tidycol[2],tidycol[3], 0.2))

polygon(y = c(baseYDens$x[1], baseYDens$x, rev(baseYDens$x)[1]),
        x = c(par('usr')[1], par('usr')[1] + baseYDens$y * densScale, par('usr')[1]),
        border="#2369bd", col=rgb(basecol[1], basecol[2],basecol[3], 0.2))

# x-densities

densScale = 0.0125

baseXDens <- density(baseFun$monthVect,
                     from=min(baseFun$monthVect),
                     to=max(baseFun$monthVect))
tidyXDens <- density(tidyFun$monthVect,
                     from=min(tidyFun$monthVect),
                     to=max(tidyFun$monthVect))
otherXDens <- density(otherFun$monthVect,
                      from=min(otherFun$monthVect),
                      to=max(otherFun$monthVect))

polygon(x = c(otherXDens$x[1], otherXDens$x, rev(otherXDens$x)[1]),
        y = c(par('usr')[3], par('usr')[3] + otherXDens$y * densScale, par('usr')[3]),
        border="grey75", 
        col=rgb(othercol[1], othercol[2],othercol[3], 0.2))

polygon(x = c(tidyXDens$x[1], tidyXDens$x, rev(tidyXDens$x)[1]),
        y = c(par('usr')[3], par('usr')[3] + tidyXDens$y * densScale, par('usr')[3]),
        border="#d12239", col=rgb(tidycol[1], tidycol[2],tidycol[3], 0.2))

polygon(x = c(baseXDens$x[1], baseXDens$x, rev(baseXDens$x)[1]),
        y = c(par('usr')[3], par('usr')[3] + baseXDens$y * densScale, par('usr')[3]),
        border="#2369bd", col=rgb(basecol[1], basecol[2],basecol[3], 0.2))

# Polygons

# with(otherFun[otherHull,],
#      polygon(y=log10(plogis(`(Intercept)`)), x=monthVect,
#              col=rgb(othercol[1], othercol[2], othercol[3], 0.2), border=NA))
# with(baseFun[baseHull,],
#      polygon(y=log10(plogis(`(Intercept)`)), x=monthVect,
#              col=rgb(basecol[1], basecol[2],basecol[3], 0.2), border=NA))
# with(tidyExtend[extendHull,],
#      polygon(y=log10(plogis(`(Intercept)`)), x=monthVect,
#              col=rgb(extendcol[1], extendcol[2],extendcol[3], 0.2), border="#d12239"))
# with(tidyFun[tidyHull,],
#      polygon(y=log10(plogis(`(Intercept)`)), x=monthVect,
#              col=rgb(tidycol[1], tidycol[2],tidycol[3], 0.2), border=NA))

with(otherFun,
     points(y=log10(plogis(`(Intercept)`)), x=monthVect, 
            pch=16, col="grey75", cex=0.5))
with(baseFun,
     points(y=log10(plogis(`(Intercept)`)), x=monthVect, 
            pch=21, bg="#2369bd", lwd=0.5, cex=0.5))
with(tidyFun,
     points(y=log10(plogis(`(Intercept)`)), x=monthVect, 
            pch=21, bg="#d12239", lwd=0.5, cex=0.5))

mtext(side=3, line=0.1, at=par("usr")[1], text="(C) Change in package use",
      adj=0, font=2)

legend(x=relative.axis.point(0.65, "x"), y=relative.axis.point(1.02, "y"),
       pch=c(21,21,21, 16), col=c("black","black","grey75"), pt.bg=c("#2369bd", "#d12239",NA),
       legend = c("Base packages", "Tidyverse","Other packages"), pt.lwd=0.5,
       bty="n", y.intersp=0.65, x.intersp=0.75)

# text(x=relative.axis.point(c(0.02), "x"), y=relative.axis.point(0.02, "y"),
#      labels=c("Used less over time"),
#      adj=0)
# 
# text(x=relative.axis.point(c(0.98), "x"),y=relative.axis.point(0.02, "y"),
#      labels="Used more over time",
#      adj=1)
# 
# Arrows(x0=relative.axis.point(0.3, "x"),
#        y0 = relative.axis.point(0.02, "y"),
#        y1 = relative.axis.point(0.02, "y"),
#        x1 = relative.axis.point(0.68, "x"),
#        arr.type="triangle", arr.width=0.1, arr.length=0.1)

box()

close.screen(3)

screen(5)
par(mar=c(0,0,0,0), ps=8, tcl=-0.25, las=1)

plot(plogis(packTop$`(Intercept)`) ~ packTop$monthVect, type="n", yaxs="i", xaxs="i",
     axes=FALSE, xlab="", ylab="", xlim=c(-0.05,0.11), ylim=log10(c(8e-5, 1.1)))

textFuns <- packTop
text(x=textFuns$monthVect,
     y=log10(plogis(textFuns$`(Intercept)`)),
     labels=textFuns$package, cex=0.75,
     col=c("grey75", "#2369bd", "#d12239", "#d12239")[as.factor(textFuns$cat)])
close.screen(5)

close.screen(all.screens=TRUE)

dev.off()

#               Summary stats ####

# how many functions used more than 10% in Dec 2021?
table(plogis(funTop$`(Intercept)`) >= 0.1) / nrow(funTop)
table(funTop$basePacks[plogis(funTop$`(Intercept)`) >= 0.1])

funTop$cat[funTop$tidyExtend] = "tidyverse"
funTop <- droplevels(funTop)

# % of functions increasing versus decreasing
funTop$sigIncr <- (funTop$monthVect - 1.96 * funTop$monthSE) > 0
funTop$sigDecr <- (funTop$monthVect + 1.96 * funTop$monthSE) < 0

# relationship between occurrence probability and slope
intSlopelm <- lm(log(plogis(funTop$`(Intercept)`)) ~ log(abs(funTop$monthVect)))
summary(intSlopelm)
plot(log(plogis(funTop$`(Intercept)`)) ~ log(abs(funTop$monthVect)))

#               Trend summary plot ####

sigChange <- cbind(table(funTop$sigDecr, funTop$cat)[2,],
                   table(!funTop$sigIncr & !funTop$sigDecr, funTop$cat)[2,],
                   table(funTop$sigIncr, funTop$cat)[2,])

sigProp <- as.data.frame.matrix(sigChange) / t(t(table(funTop$cat)))

catCentres <- table(funTop$cat)

es <- 0.01
meanChange <- cbind(table(exp(funTop$monthVect) <= 1-es,
                          funTop$cat)[2,],
                    table(exp(funTop$monthVect) > 1-es & exp(funTop$monthVect) < 1+es,
                          funTop$cat)[2,],
                    table(exp(funTop$monthVect) >= 1+es,
                          funTop$cat)[2,])

meanProp <- as.data.frame.matrix(meanChange) / t(t(table(funTop$cat)))                    

# package stats

packTop$cat[packTop$tidyExtend] = "tidyverse"
packTop <- droplevels(packTop)

packTop$sigIncr <- (packTop$monthVect - 1.96 * packTop$monthSE) > 0
packTop$sigDecr <- (packTop$monthVect + 1.96 * packTop$monthSE) < 0

pSigChange <- cbind(table(packTop$sigDecr, packTop$cat)[2,],
                    table(!packTop$sigIncr & !packTop$sigDecr, packTop$cat)[2,],
                    table(packTop$sigIncr, packTop$cat)[2,])

pSigProp <- as.data.frame.matrix(pSigChange) / t(t(table(packTop$cat)))

pMeanChange <- cbind(table(exp(packTop$monthVect) <= 1-es,
                           packTop$cat)[2,],
                    table(exp(packTop$monthVect) > 1-es & exp(packTop$monthVect) < 1+es,
                          packTop$cat)[2,],
                    table(exp(packTop$monthVect) >= 1+es,
                          packTop$cat)[2,])

pMeanProp <- as.data.frame.matrix(pMeanChange) / t(t(table(packTop$cat)))                    

# PLOT ####

boxWidth <- 0.175
offset = 0.35
xLims <- c(-0.6,1.1)

lightCols <- c("#bdd2eb", "#f1bcc3", "#ececec")[c(3,1,2)]
# lightCols <- rgb(c(basecol[1], tidycol[1], othercol[1]),
#                  c(basecol[2], tidycol[2], othercol[2]),
#                  c(basecol[3], tidycol[3], othercol[3]),
#                  0.3)[c(3,1,2)]
catCols <- rgb(c(basecol[1], tidycol[1], othercol[1]),
               c(basecol[2], tidycol[2], othercol[2]),
               c(basecol[3], tidycol[3], othercol[3]),
               1)[c(3,1,2)]

pdf("./plots/incrDec.pdf", height=7, width=6, useDingbats = FALSE)

par(mfrow=c(2,1), mar=c(0,0,0,0), oma=c(2.5,4,0.5,0.5), ps=10, tcl=-0.25, mgp=c(3,0.5,0))

plot(x=NULL, y=NULL, xlim=xLims, ylim=c(0.6,3.7), axes=FALSE, xlab="", ylab="")
segments(x0=0, x1=0,
         y0=par("usr")[3], 
         y1=relative.axis.point(0.925,"y"),#3+boxWidth + offset*0.5,
         lty="31")

axis(side=1,at=seq(-0.5,1,0.1), labels=NA)

# sigNegs
rect(xleft= - sigProp[,1] -sigProp[,2] * 0.5,
     xright= -sigProp[,2] * 0.5,
     ybottom= c(1,3,2) - boxWidth + offset * 0.5,
     ytop = c(1,3,2) + boxWidth + offset * 0.5,
     col=lightCols)

text(x = colMeans(rbind(-sigProp[,1] -sigProp[,2] * 0.5, -sigProp[,2] * 0.5)),
     y = c(1,3,2) + offset * 0.5,
     labels=sigChange[,1], col=catCols)

# sig Pos
rect(xright= sigProp[,3] + sigProp[,2] * 0.5,
     xleft= sigProp[,2] * 0.5,
     ybottom= c(1,3,2) - boxWidth + offset * 0.5,
     ytop = c(1,3,2) + boxWidth + offset * 0.5, col=lightCols)
text(x = colMeans(rbind(sigProp[,3] + sigProp[,2] * 0.5, sigProp[,2] * 0.5)),
     y = c(1,3,2)+ offset * 0.5,
     labels=sigChange[,3], col=catCols)

# sig Stable
rect(xleft= -sigProp[,2] * 0.5,
     xright= sigProp[,2] * 0.5,
     ybottom= c(1,3,2) - boxWidth + offset * 0.5,
     ytop = c(1,3,2) + boxWidth + offset * 0.5, col="white")
text(x = 0,
     y = c(1,3,2) + offset * 0.5,
     labels=sigChange[,2], col=catCols)

# mean Negs
rect(xleft= - sigProp[,1]  -sigProp[,2] * 0.5,
     xright= -sigProp[,1]  -sigProp[,2] * 0.5 + meanProp[,1],
     ybottom= c(1,3,2) - boxWidth - offset * 0.5,
     ytop = c(1,3,2) + boxWidth - offset * 0.5, col=catCols)
text(x = colMeans(rbind(- sigProp[,1]  -sigProp[,2] * 0.5, 
                  -sigProp[,1]  -sigProp[,2] * 0.5 + meanProp[,1])),
     y = c(1,3,2) - offset * 0.5,
     labels=meanChange[,1], col=lightCols)

# mean Pos
rect(xright= sigProp[,3] + sigProp[,2] * 0.5,
     xleft= sigProp[,3]  + sigProp[,2] * 0.5 - meanProp[,3],
     ybottom= c(1,3,2) - boxWidth - offset * 0.5,
     ytop = c(1,3,2) + boxWidth - offset * 0.5, col=catCols)
text(x = colMeans(rbind(sigProp[,3] + sigProp[,2] * 0.5, 
                        sigProp[,3]  + sigProp[,2] * 0.5 - meanProp[,3])),
     y = c(1,3,2) - offset * 0.5,
     labels=meanChange[,3], col=lightCols)

# mean stable text
rect(xleft= -sigProp[,1]  -sigProp[,2] * 0.5 + meanProp[,1],
     xright= sigProp[,3]  + sigProp[,2] * 0.5 - meanProp[,3],
     ybottom= c(1,3,2) - boxWidth - offset * 0.5,
     ytop = c(1,3,2) + boxWidth - offset * 0.5, col="white")
text(x = colMeans(rbind(-sigProp[,1]  -sigProp[,2] * 0.5 + meanProp[,1],
                        sigProp[,3]  + sigProp[,2] * 0.5 - meanProp[,3])),
     y = c(1,3,2) - offset * 0.5,
     labels=meanChange[,2], col=lightCols)

# labels 
text(x = - sigProp[,1]  -sigProp[,2] * 0.5,
     y = c(1,3,2) + offset * 0.5,
     labels="-", pos=2, offset=0.25)
text(x = sigProp[,3] + sigProp[,2] * 0.5,
     y = c(1,3,2) + offset * 0.5,
     labels="+", pos=4, offset=0.25)

text(x = - sigProp[,1]  -sigProp[,2] * 0.5,
     y = c(1,3,2) - offset * 0.5,
     labels="< -1%", pos=2, offset=0.25)
text(x = sigProp[,3] + sigProp[,2] * 0.5,
     y = c(1,3,2) - offset * 0.5,
     labels="> 1%", pos=4, offset=0.25)

mtext(side=2,
     at = c(1,3,2),
     text=c("Other", "Base", "tidyverse"), line=0.25, las=1,
     col=catCols, font=2)
text(x=relative.axis.point(0.01, "x"),
     y=relative.axis.point(0.95, "y"),
     labels="(A) Temporal function trends", font=2, adj=0)

box()

# PACKAGES

plot(x=NULL, y=NULL, xlim=xLims, ylim=c(0.6,3.7), axes=FALSE, xlab="", ylab="")

segments(x0=0, x1=0,
         y0=par("usr")[3],
         y1=relative.axis.point(0.925,"y"),#3+boxWidth + offset*0.5,
         lty="31")

axis(side=1,at=seq(-0.5,1,0.1), labels=NA)
axis(side=1,at=seq(-0.4,1,0.2), labels=abs(seq(-0.4,1,0.2)), 
     mgp=c(3,0.2,0))
mtext(side=1, line=1.25, text="Proportion of total", at=0)

# sigNegs
rect(xleft= - pSigProp[,1] -pSigProp[,2] * 0.5,
     xright= -pSigProp[,2] * 0.5,
     ybottom= c(1,3,2) - boxWidth + offset * 0.5,
     ytop = c(1,3,2) + boxWidth + offset * 0.5,
     col=lightCols)

text(x = colMeans(rbind(-pSigProp[,1] -pSigProp[,2] * 0.5, -pSigProp[,2] * 0.5)),
     y = c(1,3,2) + offset * 0.5,
     labels=pSigChange[,1], col=catCols)

# sig Pos
rect(xright= pSigProp[,3] + pSigProp[,2] * 0.5,
     xleft= pSigProp[,2] * 0.5,
     ybottom= c(1,3,2) - boxWidth + offset * 0.5,
     ytop = c(1,3,2) + boxWidth + offset * 0.5, col=lightCols)
text(x = colMeans(rbind(pSigProp[,3] + pSigProp[,2] * 0.5, pSigProp[,2] * 0.5)),
     y = c(1,3,2)+ offset * 0.5,
     labels=pSigChange[,3], col=catCols)

# sig Stable
rect(xleft= -pSigProp[,2] * 0.5,
     xright= pSigProp[,2] * 0.5,
     ybottom= c(1,3,2) - boxWidth + offset * 0.5,
     ytop = c(1,3,2) + boxWidth + offset * 0.5, col="white")
text(x = 0,
     y = c(1,3,2) + offset * 0.5,
     labels=pSigChange[,2], col=catCols)

# mean Negs
rect(xleft= - pSigProp[,1]  -pSigProp[,2] * 0.5,
     xright= -pSigProp[,1]  -pSigProp[,2] * 0.5 + pMeanProp[,1],
     ybottom= c(1,3,2) - boxWidth - offset * 0.5,
     ytop = c(1,3,2) + boxWidth - offset * 0.5, col=catCols)
text(x = colMeans(rbind(- pSigProp[,1]  -pSigProp[,2] * 0.5, 
                        -pSigProp[,1]  -pSigProp[,2] * 0.5 + pMeanProp[,1])),
     y = c(1,3,2) - offset * 0.5,
     labels=pMeanChange[,1], col=lightCols)

# mean Pos
rect(xright= pSigProp[,3] + pSigProp[,2] * 0.5,
     xleft= pSigProp[,3]  + pSigProp[,2] * 0.5 - pMeanProp[,3],
     ybottom= c(1,3,2) - boxWidth - offset * 0.5,
     ytop = c(1,3,2) + boxWidth - offset * 0.5, col=catCols)
text(x = colMeans(rbind(pSigProp[,3] + pSigProp[,2] * 0.5, 
                        pSigProp[,3]  + pSigProp[,2] * 0.5 - pMeanProp[,3])),
     y = c(1,3,2) - offset * 0.5,
     labels=pMeanChange[,3], col=lightCols)

# mean stable text
rect(xleft= -pSigProp[,1]  -pSigProp[,2] * 0.5 + pMeanProp[,1],
     xright= pSigProp[,3]  + pSigProp[,2] * 0.5 - pMeanProp[,3],
     ybottom= c(1,3,2) - boxWidth - offset * 0.5,
     ytop = c(1,3,2) + boxWidth - offset * 0.5, col="white")
text(x = colMeans(rbind(-pSigProp[,1]  -pSigProp[,2] * 0.5 + pMeanProp[,1],
                        pSigProp[,3]  + pSigProp[,2] * 0.5 - pMeanProp[,3])),
     y = c(1,3,2) - offset * 0.5,
     labels=pMeanChange[,2], col=lightCols)

# labels 
text(x = - pSigProp[,1]  -pSigProp[,2] * 0.5,
     y = c(1,3,2) + offset * 0.5,
     labels="-", pos=2, offset=0.25)
text(x = pSigProp[,3] + pSigProp[,2] * 0.5,
     y = c(1,3,2) + offset * 0.5,
     labels="+", pos=4, offset=0.25)

text(x = - pSigProp[,1]  -pSigProp[,2] * 0.5,
     y = c(1,3,2) - offset * 0.5,
     labels="< -1%", pos=2, offset=0.25)
text(x = pSigProp[,3] + pSigProp[,2] * 0.5,
     y = c(1,3,2) - offset * 0.5,
     labels="> 1%", pos=4, offset=0.25)

mtext(side=2,
      at = c(1,3,2),
      text=c("Other", "Base", "tidyverse"), line=0.25, las=1,
      col=catCols, font=2)

text(x=relative.axis.point(0.01, "x"),
     y=relative.axis.point(0.95, "y"),
     labels="(B) Temporal package trends", font=2, adj=0)
box()

dev.off()



# Categorical Compositional change ####

funCatProp <- lapply(split(funTableRepo, f=funTableRepo$functionCat), function(x){
  
  tempMonth <- tapply(x$prop,
                      list(x$monthsSinceJan10Created,
                           x$funAsCat),
                      sum, na.rm=TRUE)
  tempMonth[is.na(tempMonth)] = 0
  tempMonth <- tempMonth[ ,colSums(tempMonth)>0]
  tempMonth <- prop.table(tempMonth, 1)
  return(tempMonth)
  
})

#               nMDS ####

funCatMDS <- lapply(funCatProp, function(x){
  metaMDS(x, distance="bray")
})

#               Year RDA ####

catRDA <- lapply(funCatProp, function(x){
  
  dbrda(x ~ monthVect)
  
})

sapply(catRDA, function(x){
  summary(x)$cont[[1]][2,1]
})

#               Plot ####

pdf(date.wrap("./plots/nmdsCat", ".pdf"), height=10, width=3.65)

mdsLims <- list(c(-0.5,0.45),
                c(-1.5,0.9),
                c(-1.1,0.4))

library(viridisLite)
yearVect <- (as.numeric(rownames(funCatMDS[[1]]$points))-1) %/% 12
yearVect <- yearVect - min(yearVect) + 1
#yearVect <- c(1, yearVect)

par(mfrow=c(3,1), oma=c(2,2,1,1), mar=c(1,1,0.25,1), ps=10, las=1, tcl=-0.25)

sapply(c(1:3), function(n){
  
  plot(x=NULL, y=NULL, xlim=mdsLims[[n]], ylim=mdsLims[[n]], axes=FALSE, xlab="", ylab="")  
  
  if(n==1){yearBase <- colorRampPalette(c("#FFCE1D", "#2369bd"))(max(yearVect))}
  if(n==2){yearBase <- colorRampPalette(c("grey80", "black"))(max(yearVect))}
  if(n==3){yearBase <- colorRampPalette(c("#A2D022", "#D1223A"))(max(yearVect))} 
  yearCol <- yearBase[yearVect]
  
  yearPoints <- funCatMDS[[n]]$points[as.numeric(rownames(funCatMDS[[n]]$points)) %% 12 == 0, ]
  
  axis(side=1, mgp=c(3,0.1,0))
  axis(side=2, mgp=c(3,0.5,0))
  
  if(n %in% 3){
    mtext(side=1, line=1, text="nMDS 1", cex=0.8)
  }
  
  mtext(side=2, line=1.75, text="nMDS 2", las=0, cex=0.8)
  
  box()
  
  # convex hulls for 2 target years
  baseCols <- col2rgb(yearBase)/255
  
  tempOrd <- funCatMDS[[n]]
  
  sapply(4:(ncol(baseCols)+4), function(n1){
    print(n1)
    yearRows <- ((n1*12)+1):((n1+1)*12)
    
    if(sum(rownames(tempOrd$points) %in% yearRows)==0){return(NULL)}
    sub <- tempOrd$points[rownames(tempOrd$points) %in% yearRows,]
    
    subHull <- sub[chull(sub),]
    polygon(subHull, border=NA, 
            col=rgb(baseCols[1,n1-3], baseCols[2,n1-3], baseCols[3,n1-3], 0.35))
    
    text(x=colMeans(subHull)[1], colMeans(subHull)[2], labels = 2010 + (n1),
           cex=1, col=yearBase[n1-3], font=2)
    
  })
  
  sapply(2:nrow(tempOrd$points), function(n1){
    lines(tempOrd$points[c(n1-1,n1),], lwd=1, col = yearCol[n1])    
  })
  
  points(tempOrd$points, pch=21, cex=1, bg=yearCol)
  
  text(x=relative.axis.point(0.02, "x"), y=relative.axis.point(0.95, "y"),
       adj=0, labels=paste0("(",LETTERS[n],") ",
                            c("Base functions",
                              "Other functions",
                              "tidyverse functions")[n],
                            " (n = ",
                            ncol(funCatProp[[n]]), ")"), font=2)
  
  text(x=relative.axis.point(0.02, "x"), y=relative.axis.point(0.035, "y"),
       adj=0, labels=paste0("Stress = ", round(tempOrd$stress, 3)))
  
  text(x=relative.axis.point(0.87, "x"), y=relative.axis.point(0.035, "y"),
       adj=1, labels=expression("Time var = "))
  text(x=relative.axis.point(0.98, "x"), y=relative.axis.point(0.035, "y"),
       adj=1, labels=paste0(sprintf("%.2f", round(summary(catRDA[[n]])$cont[[1]][2,1]*100, 2)), "%"))
})

dev.off()

# Function diversity ####
#               Gamma diversity ####

# month by month hill numbers
funHill <- do.call("rbind", apply(funMonthProp, 1, function(x){
  data.frame(h0 = hillCalc(as.vector(x[x>0]), l=1),
             h1 = hillCalc(as.vector(x[x>0]), l=1e-6),
             h2 = hillCalc(as.vector(x[x>0]), l=-1))
}))
funHill$month <- as.numeric(rownames(funMonthProp))

# model over time
h0gam <- gam(h0 ~ s(month, bs="tp"), data=funHill, family=nb)
h1gam <- gam(h1 ~ s(month, bs="tp"), data=funHill)

#               Per-repo alpha diversity ####

funOccurS <- tapply(funTableRepo$prop,
                    funTableRepo$urID,
                    function(x){
                      
                      c(h0 = hillCalc(x, l=1),
                        h1 = hillCalc(x, l=1e-6),
                        h2 = hillCalc(x, l=-1))
                    })
funOccurSDf <- as.data.frame(do.call("rbind", funOccurS))
funOccurSDf$urID <- rownames(funOccurSDf)

funOccurS <- merge(funOccurSDf, funTableRepo[!duplicated(funTableRepo$urID), 
                                          c("urID", "monthsSinceJan10Created",
                                            "userId", "repoCount")],
                   by.x="urID", by.y = "urID", all.x=TRUE, all.y=FALSE,
                   sort=FALSE)

h0gamL <- gam(h0 ~ s(monthsSinceJan10Created, bs = "tp") + log(repoCount), 
              data=funOccurS, family=nb)
h1gamL <- gam(h1 ~ s(monthsSinceJan10Created, bs = "tp") + log(repoCount), data=funOccurS)

#               Turnover beta diversity ####

funDiss <- as.matrix(vegdist(funMonthProp))
funTurn <- diag(funDiss[-1,-ncol(funDiss)])
funBase <- funDiss[1,-1]

funHill$funTurn <- c(NA, funTurn)
funHill$funBase <- c(NA, funBase)

Turngam <- gam(funTurn ~ s(month, bs="tp"), data=funHill)
Basegam <- gam(funBase ~ s(month, bs="tp"), data=funHill)

#               Save models ####

saveRDS(list(list(funHill, h0gam, h1gam),
             list(funOccurS, h0gamL, h1gamL),
             list(funHill, Turngam, Basegam)),
        "./outputs/hillDiversityModels.rds")

divList <- readRDS("./outputs/hillDiversityModels.rds")
funHill <- divList[[1]][[1]]
h0gam <- divList[[1]][[2]]
h1gam <- divList[[1]][[3]]
funOccurS <- divList[[2]][[1]]
h0gamL <- divList[[2]][[2]]
h1gamL <- divList[[2]][[3]]

colnames(funOccurS)[colnames(funOccurS)=="month"] = "monthsSinceJan10Created"
colnames(funOccurS)[colnames(funOccurS)=="scriptCount"] = "repoCount"

#               Plot ####

# gamma diversity
pred.df <- data.frame(month = seq(min(funHill$month), max(funHill$month), 1))
h0P <- cbind(pred.df, as.data.frame(predict(h0gam, newdata=pred.df, se.fit=TRUE)))
h1P <- cbind(pred.df, as.data.frame(predict(h1gam, newdata=pred.df, se.fit=TRUE)))
TurnP <- cbind(pred.df, as.data.frame(predict(Turngam, newdata=pred.df, se.fit=TRUE)))
BaseP <- cbind(pred.df, as.data.frame(predict(Basegam, newdata=pred.df, se.fit=TRUE)))

pred.dfL <- data.frame(monthsSinceJan10Created = seq(min(funOccurS$monthsSinceJan10Created), 
                                                     max(funOccurS$monthsSinceJan10Created), 1),
                       repoCountL = mean(log(funOccurS$repoCount)),
                       repoCount = exp(mean(log(funOccurS$repoCount))))

h0PL <- cbind(pred.dfL, as.data.frame(predict(h0gamL, newdata=pred.dfL, se.fit=TRUE)))
h1PL <- cbind(pred.dfL, as.data.frame(predict(h1gamL, newdata=pred.dfL, se.fit=TRUE)))

colnames(h0PL)[1] = "month"
colnames(h1PL)[1] = "month"

colnames(funOccurS)[colnames(funOccurS) == "monthsSinceJan10Created"] = "month"

pdf(date.wrap("./plots/function_diversity",".pdf"), 
        height=4.5, width=8, useDingbats=FALSE)

par(mfcol=c(2,3), mar=c(0,3,0,1), oma=c(2,1,2,0), ps=10, tcl=-0.25)

ylims <- list(c(1225,1700), c(750,1400), 
              c(27,32), c(21,26), 
              c(0.34,0.45), c(0.1,0.45))

sapply(1:4, function(n){
  print(n)
  temp.pred <- list(h0P, h1P, h0PL, h1PL, BaseP, TurnP)[[n]]
  temp.raw <- list(funHill[,c("h0", "month")], 
                   funHill[,c("h1", "month")], 
                   funOccurS[,c("h0", "month")], 
                   funOccurS[,c("h1", "month")],
                   funHill[,c("funBase", "month")],
                   funHill[,c("funTurn", "month")])[[n]]
  
  if(n %in% 3:4){
    temp.raw <- cbind(tapply(temp.raw[,1], temp.raw[,2], mean),
                      sort(unique(temp.raw[,2])))
  }
  
  plot(x=NULL, y=NULL, xlim=range(funHill$month) + c(0, 20), ylim=ylims[[n]], axes=FALSE,
       ylab="", xlab="")
  
  if(n %in% c(1,2)){
    axis(side=2, las=1, mgp=c(3,0.5,0),
         at=seq(500,2000,100),
         labels=format(seq(500,2000,100), big.mark=","))
  } else {
    axis(side=2, las=1, mgp=c(3,0.5,0))
  }
  
  mtext(side=2, line=ifelse(n %in% 1:2, 2.75, 2.5), 
        text=c("Unweighted (Hill q = 0)",
               "Rare functions down-weighted (Hill q = 1)",
               "Unweighted (Hill q = 0)",
               "Rare functions down-weighted (Hill q = 1)",
               "Dissimilarity from Jan 2014",
               "Dissimilarity from previous month")[n], las=0, cex=0.8)
  
  axis(side=1, at=seq(37,145,12), labels=NA)
  
  if(n %in% c(2,4,6)){
    par(xpd=NA)
    text(x=seq(49,145,12),
         y=relative.axis.point(-0.04,"y"),
         labels=2014:2022, srt=30, adj=1)
    par(xpd=FALSE)
    
  }
  
  if(n %in% c(1,3,5)){
    mtext(side=3, text=c("Global function diversity","",
                         "Per-script function diversity","",
                         "Change in function use")[n],
          font=2, cex=0.8)
  }
  
  yearVect <- temp.raw[,2] %/% 12
  yearVect <- yearVect - min(yearVect) + 1
  yearBase <- viridis(max(yearVect), option="C", end=0.8)
  yearCol <- yearBase[yearVect]
  
  if(n == 1){targets = seq(-0.5, 0.5, 0.05)}
  if(n == 2){targets = seq(-0.5, 1, 0.1)}
  if(n %in% c(3:4)){targets = seq(-0.1,0.15,0.02)}
  if(n == 5){targets = seq(-0.8,1,0.05)}
  if(n == 6){targets = seq(-0.8,1,0.1)}
  
  sapply(targets, function(x){
    
    if(n %in% c(1,3)){
      xpos = exp(temp.pred$fit[1]) + x * exp(temp.pred$fit[1])
    } else {
      xpos = temp.pred$fit[1] + x * temp.pred$fit[1]
    }
    
    segments(y0=xpos,
             y1=xpos,
             x0=relative.axis.point(0, "x"),
             x1=relative.axis.point(0.85, "x"), col="grey", lty="31")
    text(x=relative.axis.point(0.85, "x"),
         y=xpos,
         labels=paste0(x*100,"%"), pos=4, col="grey", offset=0.25)
  })
  
  if(n %in% c(1:2, 5:6)){
    points(x=temp.raw[,2], y=temp.raw[,1], pch=16, col=yearCol, cex=0.75)
    
    segments(x0=temp.raw[-nrow(temp.raw),2], x1=temp.raw[-1,2],
             y0=temp.raw[-nrow(temp.raw),1], y1=temp.raw[-1,1],
             pch=16, col=yearCol, lwd=1)
    
    rect(xleft=par("usr")[1], xright=par("usr")[2], 
         ybottom=par("usr")[3], ytop=par("usr")[4],
         border=NA, col=rgb(1,1,1,0.3))
  }
  
  if(n %in% c(1,3)){
    polygon(x=c(temp.pred[,1], rev(temp.pred[,1])),
            y=exp(c(temp.pred[,"fit"] + 1.96*temp.pred[,"se.fit"], rev(temp.pred[,"fit"] - 1.96 * temp.pred[,"se.fit"]))),
            border=NA, col=rgb(0.5,0.5,0.5,0.5))
    lines(exp(temp.pred[,"fit"]) ~ temp.pred[,"month"], lwd=2)
    
  } else {
    polygon(x=c(temp.pred[,"month"], rev(temp.pred[,"month"])),
            y=c(temp.pred[,"fit"] + 1.96*temp.pred[,"se.fit"], rev(temp.pred[,"fit"] - 1.96 * temp.pred[,"se.fit"])),
            border=NA, col=rgb(0.5,0.5,0.5,0.5))
    lines(temp.pred[,"fit"] ~ temp.pred[,"month"], lwd=2)
  }
  
  text(x=relative.axis.point(0.05, "x"),
       y=relative.axis.point(0.95, "y"),
       labels=paste0("(", LETTERS[n],")"), font=2)
  
  box()
  
})

dev.off()

#               Are tidyverse models more rich? ####

hasTidy <- as.data.frame(tapply(funTableRepo$functionCat == "tidyverse", funTableRepo$urID, sum))
colnames(hasTidy) = "tidyCount"
hasTidy$urID <- rownames(hasTidy)
summary(hasTidy$tidyCount > 0)

funOccurS$hasTidy <- as.factor((hasTidy$tidyCount > 0)[match(hasTidy$urID, funOccurS$urID)])

tidyRich <- glmer(h0 ~ hasTidy + log(repoCount) + (1|userId), family=poisson, data=funOccurS)
tidyRich0 <- glm(h0 ~ hasTidy + log(repoCount), family=poisson, data=funOccurS)
tidyRich1 <- lm(log(h1) ~ hasTidy + log(repoCount), data=funOccurS)
performance(tidyRich0)
summary(tidyRich0)

# predictions
predDf <- data.frame(hasTidy = c("FALSE", "TRUE"),
           repoCount = mean(log(funOccurS$repoCount)))
predDf <- cbind(predDf,
                as.data.frame(predict(tidyRich0, newdata=predDf, se.fit=TRUE)),
                as.data.frame(predict(tidyRich1, newdata=predDf, se.fit=TRUE)))

pdf("./plots/hasTidy.pdf", height=4, width=3, useDingbats = FALSE)

par(mar=c(2,3,0.5,0.5), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)
plot(x=NULL, y=NULL, xlim=c(0.5,2.5), ylim=c(1.5,2.25), xaxt="n", yaxt="n", xlab="", ylab="")

axis(side=1, at=(1:2) - 0.25, labels = rep("Hill 0", 2), mgp=c(3,0.1,0)) 
axis(side=1, at=(1:2) + 0.25, labels = rep("Hill 1", 2), mgp=c(3,0.1,0))
mtext(side=1, at=c(1,2), line=1, text=c("No tidyverse", "Tidyverse"))
mtext(side=2, line=1.5, las=0, text="Function richness (n = )")

axis(side=2, at=log(c(1:15,10,100)), labels = c(1:15,10,100), tcl=-0.25)

axis(side=2, at=log(c(seq(1,10,1),
                      seq(10,100,10),
                      seq(100,1000,100))), 
     labels=NA, tcl=-0.125)

segments(x0=1:2 - 0.25, x1=1:2 - 0.25,
         y0=predDf[,3] + 1.96 * predDf[,4],
         y1=predDf[,3] - 1.96 * predDf[,4])

segments(x0=1:2 + 0.25, x1=1:2 + 0.25,
         y0=predDf[,6] + 1.96 * predDf[,7],
         y1=predDf[,6] - 1.96 * predDf[,7])

segments(x0=0.75, x1=1.75, y0=predDf[1,3], y1=predDf[1,3], lty="31", col="grey50")
segments(x0=1.75, x1=1.75, y0=predDf[1,3], y1=predDf[2,3], lty="31", col="grey50")
segments(x0=1.25, x1=2.25, y0=predDf[1,6], y1=predDf[1,6], lty="31", col="grey50")
segments(x0=2.25, x1=2.25, y0=predDf[1,6], y1=predDf[2,6], lty="31", col="grey50")

text(x=1.75, y=mean(predDf[,3]), pos=4,
     labels=paste0(round(exp(predDf[2,3]) / exp(predDf[1,3]),2), "x"))
text(x=2.25, y=mean(predDf[,6]), pos=2,
     labels=paste0(round(exp(predDf[2,6]) / exp(predDf[1,6]),2), "x"))

points(y=predDf[,3], x=(1:2) - 0.25, pch=21, bg=c("#2369bd", "#d12239"))
points(y=predDf[,6], x=(1:2) + 0.25, pch=21, bg=c("#2369bd", "#d12239"))

dev.off()

# Category diversity ####
#               Prop data ####

# month by month hill numbers
funHillCat <- do.call("rbind", lapply(1:length(funCatProp), function(n){
  
  y<-funCatProp[[n]]
  
  temp <- do.call("rbind", apply(y, 1, function(x){
  data.frame(h0 = hillCalc(as.vector(x[x>0]), l=1),
             h1 = hillCalc(as.vector(x[x>0]), l=1e-6),
             h2 = hillCalc(as.vector(x[x>0]), l=-1))
  }))
  
  temp$month <- as.numeric(rownames(y))# wide form data 
  temp$functionCat <- levels(funTableRepo$functionCat)[n]
  return(temp)
}))

# per repo diversity

funOccurSCat <- tapply(funTableRepo$prop,
                    list(as.factor(funTableRepo$urID),
                         funTableRepo$functionCat),
                    function(x){
                      
                      c(h0 = hillCalc(x, l=1),
                        h1 = hillCalc(x, l=1e-6),
                        h2 = hillCalc(x, l=-1))
                      
                    }, simplify=FALSE)

funOccurSCatDf <- do.call("rbind", lapply(1:nrow(funOccurSCat), function(n){
  
  x <- funOccurSCat[n,]
  temp <- as.data.frame(t(sapply(x, function(x1){
    if(is.null(x1)){return(c(0,0,0))} else{return(x1)}
  })))
  colnames(temp) = c("h0", "h1", "h2")
  temp$functionCat <- rownames(temp)
  temp$urID <- rownames(funOccurSCat)[n]
  rownames(temp)<-NULL
  return(temp)
  
}))


funOccurSCatDf <- merge(funOccurSCatDf, funTableRepo[!duplicated(funTableRepo$urID), 
                                             c("urID", "monthsSinceJan10Created",
                                               "userId", "repoCount")],
                   by.x="urID", by.y = "urID", all.x=TRUE, all.y=FALSE,
                   sort=FALSE)

#               Gamma diversity ####
  
funHillCat$functionCat <- as.factor(funHillCat$functionCat)

# model over time
h0CatGam <- gam(h0 ~ s(month, bs="tp", by=functionCat) + functionCat, data=funHillCat, family=nb)
h1CatGam <- gam(h1 ~ s(month, bs="tp", by=functionCat) + functionCat, data=funHillCat)

#               Alpha diversity ####

funOccurSCatDf$functionCat <- as.factor(funOccurSCatDf$functionCat)

# Two step process for repo diversity
funOccurSCatDf$h1Bin <- ifelse(funOccurSCatDf$h1 > 0, 1, 0)

h1CatgamBin <- gam(h1Bin ~ s(monthsSinceJan10Created, bs = "tp", by=functionCat) + functionCat + log(repoCount), 
                 data=funOccurSCatDf, family=binomial)

saveRDS(h1CatgamBin, "./outputs/CatAlphaBinaryModel.rds")

funOccurNonZero <- droplevels(funOccurSCatDf[funOccurSCatDf$h1Bin > 0, ])

h1CatgamNonZero <- gam(h1 ~ s(monthsSinceJan10Created, bs = "tp", by=functionCat) + functionCat + log(repoCount), 
                       data=funOccurNonZero, family=Gamma)

saveRDS(h1CatgamNonZero, "./outputs/CatAlphaGammaModel.rds")

#               Model predictions ####

# read in alpha models
h1CatgamBin <- readRDS("./outputs/CatAlphaBinaryModel.rds")
h1CatgamNonZero <- readRDS("./outputs/CatAlphaGammaModel.rds")

monthSeq <- seq(min(funOccurSCatDf$monthsSinceJan10Created), 
                max(funOccurSCatDf$monthsSinceJan10Created), 1)

# gamma predictions

pred.df <- data.frame(month = rep(monthSeq, 3),
                      functionCat = rep(levels(as.factor(funOccurSCatDf$functionCat)), each=length(monthSeq)))
hoCatPred <- cbind(pred.df, as.data.frame(predict(h0CatGam, newdata=pred.df, se.fit=TRUE)))
h1CatPred <- cbind(pred.df, as.data.frame(predict(h1CatGam, newdata=pred.df, se.fit=TRUE)))
funOccurNonZero <- droplevels(funOccurSCatDf[funOccurSCatDf$h1Bin > 0, ])

# alpha predictions

pred.dfL <- data.frame(monthsSinceJan10Created = rep(monthSeq, 3),
                       repoCountL = mean(log(funOccurNonZero$repoCount)),
                       repoCount = exp(mean(log(funOccurNonZero$repoCount))),
                       functionCat=rep(levels(as.factor(funOccurSCatDf$functionCat)), each=length(monthSeq)))
h1CatBinPL <- cbind(pred.dfL, as.data.frame(predict(h1CatgamBin, newdata=pred.dfL, se.fit=TRUE)))
h1CatGammaPL <- cbind(pred.dfL, as.data.frame(predict(h1CatgamNonZero, newdata=pred.dfL, se.fit=TRUE)))
colnames(h1CatBinPL)[1] = "month"
colnames(h1CatBinPL)[1] = "month"

#               Plot ####

pdf(date.wrap("./plots/functionDiversityByCat", ".pdf"), 
    height=4.5, width=8, useDingbats=FALSE)

par(mfcol=c(2,3), mar=c(0,3,0,1), oma=c(2,1,2,0), ps=10, tcl=-0.25)

ylims <- list(log(c(140,650)), log(c(70,500)), 
              c(0,1.1), c(1,11.5))

sapply(1:4, function(n){
  
  print(paste0("--", n))
  temp.pred <- list(hoCatPred, h1CatPred, h1CatBinPL, h1CatGammaPL)[[n]]
  temp.raw <- list(funHillCat[,c("h0", "month", "functionCat")], 
                   funHillCat[,c("h1", "month", "functionCat")], 
                   funOccurSCatDf[,c("h1", "monthsSinceJan10Created", "functionCat")], 
                   funOccurSCatDf[,c("h1", "monthsSinceJan10Created", "functionCat")])[[n]]
  
  colnames(temp.pred)[colnames(temp.pred)=="functionCat"] = "cat"
  colnames(temp.raw)[colnames(temp.raw)=="functionCat"] = "cat"
  
  colnames(temp.raw)[2] = "month"
  
  plot(x=NULL, y=NULL, xlim=range(funHillCat$month) + c(0, 25), ylim=ylims[[n]], axes=FALSE,
       ylab="", xlab="")
  
  if(n < 3){
    
    axis(side=2, las=1, mgp=c(3,0.5,0),
         at=log(c(seq(1,20,1),
                  seq(20,200,10),
                  seq(200,1000,100),
                  seq(1000,10000,1000))),
         labels=NA, tcl=-0.125)
    
    if(n==1){
      
      axis(side=2, las=1, mgp=c(3,0.5,0),
           at=log(seq(100,2000,100)),
           labels=format(seq(100,2000,100), big.mark=","))
    
      } else {
      
        axis(side=2, las=1, mgp=c(3,0.5,0),
           at=log(c(1:5,
                    seq(10,100,10),
                    seq(100,2000,100))),
           labels=format(c(1:5,
                           seq(10,100,10),
                           seq(100,2000,100)), big.mark=","))
      
    }
    } else {
      axis(side=2, las=1, mgp=c(3,0.5,0))
    }
    
    
  mtext(side=2, line=2.5, 
        text=c("Unweighted (Hill q = 0)",
               "Rare functions down-weighted (Hill q = 1)",
               "Probability of occurrence",
               "Rare functions down-weighted (Hill q = 1)")[n], las=0, cex=0.8)
  
  axis(side=1, at=seq(49,145,12), labels=NA)
  
  if(n %in% c(2,4,6)){
    par(xpd=NA)
    text(x=seq(49,145,12),
         y=relative.axis.point(-0.04,"y"),
         labels=2014:2022, srt=30, adj=1)
    par(xpd=FALSE)
    
  }
  
  if(n %in% c(1,3,5)){
    mtext(side=3, text=c("Global function diversity","",
                         "Per-script function diversity","")[n],
          font=2, cex=0.8)
  }
  
  predSplit <- split(temp.pred, f=temp.pred$cat)
  rawSplit <- split(temp.raw, f=temp.raw$cat)
  
  sapply(1:3, function(n1){
    
    print(n1)
    subPred <- predSplit[[n1]]
    subRaw <- rawSplit[[n1]]
    
    pointCol <- col2rgb(c("#2369bd", "grey75", "#d12239")[n1])/255
    faintCol <- col2rgb(c("#8eb8ea", "grey90", "#ef9da7")[n1])/255
    linecol <- col2rgb(c("black", "grey75", "#d12239")[n1])/255
    
    if(!n %in% c(3:4)){
      
    tempPoints <- log(subRaw[,1])
    
    points(tempPoints ~ subRaw[,2],
           pch=21, 
           bg=rgb(faintCol[1], faintCol[2], faintCol[3], 1),
           col=rgb(linecol[1], linecol[2], linecol[3], 0.5), cex=0.7)
    
    segments(x0=subRaw[-nrow(subRaw),2], x1=subRaw[-1,2],
             y0=tempPoints[-nrow(subRaw)], y1=tempPoints[-1],
             pch=16, col=rgb(faintCol[1], faintCol[2], faintCol[3], 1), lwd=1)
    }
    
    # log link
    if(n == 1){
    
      polygon(x=c(subPred$month, rev(subPred$month)),
              y=c(subPred$fit + 1.96 * subPred$se.fit,
                  rev(subPred$fit - 1.96 * subPred$se.fit)),
              border=NA, col=rgb(pointCol[1], pointCol[2], pointCol[3], 0.5)) 
      
    if(n1 == 3){pointCol <- linecol}
    lines(subPred$fit ~ subPred$month,
          col=rgb(pointCol[1], pointCol[2], pointCol[3]), lwd=2,
          lty=c("solid", "solid", "solid")[n1])
    
    propChange <- exp(rev(subPred$fit)[1]) / exp(subPred$fit[1])
    
    text(x=rev(subPred$month)[1],
         y=rev(subPred$fit)[1],
         pos=4,
         labels=paste0(c("B", "O", "T")[n1], " (", sprintf("%.1f", round(propChange,1)),"x)"),
         font=2,
         col=rgb(pointCol[1], pointCol[2], pointCol[3]))
    }
    
    # logit link
    if(n == 3){
      
      polygon(x=c(subPred$month, rev(subPred$month)),
              y=c(plogis(subPred$fit + 1.96 * subPred$se.fit),
                  rev(plogis(subPred$fit - 1.96 * subPred$se.fit))),
              border=NA, col=rgb(pointCol[1], pointCol[2], pointCol[3], 0.5)) 
      
      if(n1 == 3){pointCol <- linecol}
      lines(plogis(subPred$fit) ~ subPred$month,
            col=rgb(pointCol[1], pointCol[2], pointCol[3]), lwd=2,
            lty=c("solid", "solid", "solid")[n1])
      
      propChange <- plogis(rev(subPred$fit)[1]) / plogis(subPred$fit[1])
      text(x=rev(subPred$month)[1],
           y=rev(plogis(subPred$fit))[1],
           pos=4,
           labels=paste0(c("B", "O", "T")[n1], " (", sprintf("%.1f", round(propChange,1)),"x)"),
           font=2,
           col=rgb(pointCol[1], pointCol[2], pointCol[3]))
      
    }
    
    # inverse link
    if(n == 4){
      
      polygon(x=c(subPred$month, rev(subPred$month)),
              y=1/(c(subPred$fit + 1.96 * subPred$se.fit,
                  rev(subPred$fit - 1.96 * subPred$se.fit))),
              border=NA, col=rgb(pointCol[1], pointCol[2], pointCol[3], 0.3)) 
      
      if(n1 == 3){pointCol <- linecol}
      lines(1/(subPred$fit) ~ subPred$month,
            col=rgb(pointCol[1], pointCol[2], pointCol[3]), lwd=2)
      
      propChange <- (1/rev(subPred$fit)[1]) / (1/(subPred$fit[1]))
      text(x=rev(subPred$month)[1],
           y=rev(1/(subPred$fit))[1],
           pos=4,
           labels=paste0(c("B", "O", "T")[n1], " (", sprintf("%.1f", round(propChange,1)),"x)"),
           font=2,
           col=rgb(pointCol[1], pointCol[2], pointCol[3]))
      
    }
    
    # identity link
    if(n == 2){
      
      polygon(x=c(subPred$month, rev(subPred$month)),
              y=log(c(subPred$fit + 1.96 * subPred$se.fit,
                     rev(subPred$fit - 1.96 * subPred$se.fit))),
              border=NA, col=rgb(pointCol[1], pointCol[2], pointCol[3], 0.3)) 
      
      if(n1 == 3){pointCol <- linecol}
      lines(log(subPred$fit) ~ subPred$month,
            col=rgb(pointCol[1], pointCol[2], pointCol[3]), lwd=2,
            lty=c("solid", "solid", "solid")[n])
      
      propChange <- rev(subPred$fit)[1] / (subPred$fit[1])
      text(x=rev(subPred$month)[1],
           y=rev(log(subPred$fit))[1],
           pos=4,
           labels=paste0(c("B", "O", "T")[n1], " (", sprintf("%.1f", round(propChange,1)),"x)"),
           font=2,
           col=rgb(pointCol[1], pointCol[2], pointCol[3]))
      
    }
    
})
  
  text(x=relative.axis.point(0.05, "x"),
       y=relative.axis.point(0.95, "y"),
       labels=paste0("(", LETTERS[n],")"), font=2)
  
  box()
  
})

dev.off()

# function frequency distributions ####
pdf("./plots/wordFreqDist.pdf", height=3.5, width=7)

split.screen(rbind(c(0.1,0.65,0.15,0.99),
                   c(0.75,0.99,0.15,0.99)))

screen(1)
par(mar=c(0,0,0,0), ps=10, tcl=-0.25, las=1, mgp=c(3,0.5,0))

plot(x=NULL, y=NULL, xlim=c(-15,-4), ylim=c(-0.02,0.75), yaxs="i", xaxt="n",
     xlab="", ylab="")
mtext(side=2, line=2, text="Density", las=0)

xTicks <- sapply(as.numeric(paste0("1e", (-10:-1))),
                 function(n){seq(n, n*10, n)})
axis(side=1, at=log(as.vector(xTicks)), tcl=-0.125, labels=NA)
axis(side=1, at=log(as.numeric(paste0("1e", (-10:-1)))), 
     labels=as.numeric(paste0("1e", (-10:-1))),
     mgp=c(3,0.2,0))
mtext(side=1, line=1.5, text="Function relative abundance")

sapply(1:nrow(funMonthProp), function(n){
  lX <- log(funMonthProp[n,][funMonthProp[n,]>0])
  dlX <- density(lX, from=min(lX), to=max(lX))
  lines(x=dlX$x, y=dlX$y, col=yearCol[n], lwd=0.25)
  
})

# find normal distribution that best describes data
funPropVect <- as.vector(funMonthProp)
funPropVect <- funPropVect[funPropVect>0]

funlm <- lm(log(funPropVect) ~ 1)
summary(funlm)
funPred <- predict(funlm, interval="prediction", level=0.67)[1,]
funSeq <- seq(min(log(funPropVect)),max(log(funPropVect)),len=200)

funQuan <- dnorm(x=funSeq, mean=funPred[1], sd= funPred[3] - funPred[1])

lines(funQuan[funQuan > 1e-5] ~ funSeq[funQuan > 1e-5], lwd=2, lty="31")

text(x=relative.axis.point(0.02, "x"),
     y=relative.axis.point(0.965, "y"),
     labels="(A)", font=2, adj=0)

close.screen(1)

screen(2)
par(mar=c(0,0,0,0), ps=10, tcl=-0.25, las=1, mgp=c(3,0.5,0))

a <- qqnorm(residuals(funlm), plot.it=FALSE)

plot(a, pch=16, cex=0.5, xlab="", ylab="", xaxt="n", col="grey65")


axis(side=1, mgp=c(3,0.2,0))
mtext(side=1, line=1.5, text="Theoretical quantiles")
mtext(side=2, line=2, text="Observed data quantiles", las=0)

b <- qqline(residuals(funlm), col="black", lwd=2, lty="31")

text(x=relative.axis.point(0.02, "x"),
     y=relative.axis.point(0.965, "y"),
     labels="(B)", font=2, adj=0)
close.screen(2)


dev.off()

# Function extinction/birth ####

# first occurrence
firstOccur <- min(as.numeric(rownames(funMonthProp))) - 1 + apply(funMonthProp, 2, function(x){which(x>0)[1]})

# last occurrence
lastOccur <- max(as.numeric(rownames(funMonthProp))) + 1 - apply(funMonthProp, 2, function(x){which(rev(x)>0)[1]})

plot(firstOccur ~ lastOccur)

summary(firstOccur == min(firstOccur) & lastOccur == max(lastOccur))

# how many zero occurrences
zeroOccur <-  apply(funMonthProp, 2, function(x){sum(x == 0)})
table(zeroOccur)

extStats <- data.frame(function. = colnames(funMonthProp),
                       firstOccur = firstOccur,
                       lastOccur = lastOccur,
                       zeroOccur = zeroOccur,
                       averageProp = colMeans(funMonthProp))

table(extStats$lastOccur)

extStats[extStats$lastOccur < 140,]

head(extStats, 50)

plot(funMonthProp[,"left_join"], ylim=c(0,0.002))

# ####
# Vinettes of change ####

#               Import functions ####

ylims <- log(c(1e-5, 1))
mainYlims=c(0,1)
funSub <- funTableRepo


importFuns <- functionVignetteRepo(c("read.csv", "read.table", "read_csv", "read.delim", "read_excel",
                                 "read_xlsx", "read.xlsx", "read.csv2", "read_tsv", "read_delim",
                                 "read.xls", "fread"),
                 droplevels(funTableRepo),
                 "ImportFunctions2",
                 ylims=ylims,
                 mainYlims=mainYlims)

saveRDS(importFuns,"./outputs/importFunsFull.rds")


#               Merge functions ####

joinFuns <- functionVignetteRepo(c("merge", "inner_join", "left_join", "right_join", 
                   "full_join", "semi_join", "anti_join", "merge.data.table",
                   "join", "join_all"),
                 funTableRepo,
                 "JoinFunctions",
                 ylims=ylims,
                 mainYlims=mainYlims)

#               reshape functions #####

reshapeFuns <- functionVignetteRepo(c("reshape", "pivot_longer", "pivot_wider",
                                      "melt", "cast", "gather" ,"spread", "xtabs",
                                      "dcast", "acast", "stack", "unstack"),
                                 funTableRepo,
                                 "reshapeFuns",
                                 ylims=ylims,
                                 mainYlims=mainYlims)


#             models ####

modelFuns <- functionVignetteRepo(c("lm", "glm", "aov", "lmer", "gam", "gamm", "rq", "naiveBayes", "rpart", "auto.arima",
                                    "ctree", "glmer.nb", "lda", "caret",
                                   "glmer", "brms", "randomForest", "train", "stan", "glmnet", "knn", "t.test"),
                                 funTableRepo,
                                 "modelFunctions",
                                 ylims=ylims,
                                 mainYlims=mainYlims)

saveRDS(joinFuns, "./outputs/joinFuns.rds")

#           plotting ####

modelFuns <- functionVignetteRepo(c("plot", "ggplot"),
                                  funTableRepo,
                                  "plotFunctions",
                                  ylims=log(c(0.1,0.5)),
                                  mainYlims=c(0,0.8))

sort(table(funTableSub$function.[funTableSub$package=="dplyr"]), decreasing=TRUE)

modelFuns <- functionVignetteRepo(c("ddply", "ldply", "mapvalues", "llply", "join_all",
                                    "filter", "select", "group_by", "left_join", "bind_rows", "distinct", "funs", "case_when"),
                                  funTableRepo,
                                  "plyrFunctions",
                                  ylims=log(c(1e-6,0.5)),
                                  mainYlims=c(0,1))
