functionVignette <- function(importFuns, 
                             funSub, 
                             plotName, 
                             ylims, mainYlims){

funSub$monthsSinceJan10CreatedFact <- as.factor(funSub$monthsSinceJan10Created)
funSub$urIDF <- as.factor(funSub$urID)

scriptCount <- table(funSub$urIDF[!duplicated(funSub$ursID)])
  
funSub <- funSub[funSub$function. %in% importFuns,]

importDf <- as.data.frame(table(funSub$urIDF[!duplicated(funSub$ursID)]))
importDf$total <- as.vector(scriptCount)
importDf <- merge(x=importDf, y=funSub[!duplicated(funSub$urID), 
                                            c("urID", "userId", "monthsSinceJan10Created")],
                  by.x="Var1", by.y="urID", all.x=TRUE, all.y=FALSE, sort=FALSE)
colnames(importDf) = c("urID", "count", "total", "userId", "monthNum")
importDf$fail <- importDf$total - importDf$count

importList <- lapply(importFuns, function(x){
  print(x)

  subDf <- funSub[funSub$function.==x,]
  subTable <- as.data.frame(table(subDf$urIDF))
  subTable$total <- scriptCount
  
  subTable <- merge(x=subTable, y=funSub[!duplicated(funSub$urID), 
                                              c("urID", "userId", "monthsSinceJan10Created")],
                    by.x="Var1", by.y="urID", all.x=TRUE, all.y=FALSE, sort=FALSE)
  colnames(subTable) = c("urID", "count", "total", "userId", "monthNum")
  subTable$fail <- subTable$total - subTable$count
  
  return(subTable)
})
  
  #                  Models over time & plot ####

monthTotals <- table(funSub$monthsSinceJan10Created[!duplicated(funSub$urID)])
monthNum <- as.numeric(names(monthTotals))

pdf(date.wrap(paste0("./plots/",plotName),".pdf"), 
    height=7, width=9.5, useDingbats = FALSE)

split.screen(rbind(c(0.1,0.5,0.65,0.99),
                   c(0.1,0.5,0.1,0.65),
                   c(0.1,0.5,0.675,0.775)))

screen(1)
par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)
plot(x=NULL, y=NULL, xlim=c(24,165), ylim=mainYlims,
     axes=FALSE, xlab="", ylab="")
axis(side=2)
axis(side=1, at=seq(49,12*12+1, 12),
     labels=NA)
mtext(side=2, line=2.25, text="Probability of any function\noccurring in a repository", las=0)

tempGam <- gam(cbind(count, fail) ~ s(monthNum, k = 6), data=importDf, family=binomial, method="REML")

library(gamm4)
tempGamm <- gamm4(cbind(count, fail) ~ s(monthNum), data=importDf, family=binomial,
                  random = ~(1|userId))
summary(tempGamm$gam)

pred.df <- data.frame(monthNum=seq(min(monthNum), max(monthNum), len=200))
pred.df <- cbind(pred.df,
             as.data.frame(predict(tempGamm$gam, newdata=pred.df, se.fit=TRUE)))

polygon(x=c(pred.df$monthNum, rev(pred.df$monthNum)),
        y=plogis(c(pred.df$fit + 1.96 * pred.df$se.fit,
                       rev(pred.df$fit - 1.96 * pred.df$se.fit))),
        border=NA, col=rgb(0, 0, 0, 0.35))
lines(plogis(pred.df$fit) ~ pred.df$monthNum, lwd=2)
text(x=pred.df$monthNum[nrow(pred.df)],
     y=plogis(pred.df$fit)[nrow(pred.df)],
     labels=paste0(round(plogis(pred.df$fit)[nrow(pred.df)]/
                           plogis(pred.df$fit)[1],2),"x"),
     pos=4, cex=0.8)
text(x=relative.axis.point(0.025,"x"),
     y=relative.axis.point(0.935, "y"),
     labels=c("(A) Function category probability over time"), font=2, adj=0)

segments(x0=min(monthNum), x1=min(monthNum),
         y0=relative.axis.point(0.25,"y"),
         y1=plogis(pred.df$fit)[1], lty="31")

segments(x0=max(monthNum), x1=max(monthNum),
         y0=relative.axis.point(0.25,"y"),
         y1=rev(plogis(pred.df$fit))[1], lty="31")

box()
close.screen(1)

screen(2)
par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)
plot(x=NULL, y=NULL, xlim=c(24,165), ylim=ylims,
     axes=FALSE, xlab="", ylab="")

axis(side=2, at=log(c(0.1,0.01,1e-3,1e-4)),
     labels=c(0.1,0.01,1e-3,1e-4))
mtext(side=2, line=3, text="Probability of function occurring in a repository", las=0)
mtext(side=1, line=2, text="Time of repository creation", las=0)

axis(side=1, at=seq(49,12*12+1, 12),
     labels=NA)
par(xpd=NA)
text(x=seq(49,12*12+1, 12), y=relative.axis.point(-0.025, "y"),
     labels=2014:2022, srt=30, adj=c(1))
par(xpd=FALSE)

axis(side=2, at=log(c(seq(0.1,1,0.1),
                      seq(0.01,0.1,0.01),
                      seq(1e-3,0.01,1e-3),
                      seq(1e-4, 1e-3, 1e-4))),
     labels=NA, tcl=-0.125)

importMods <- lapply(1:length(importList), function(n){
print(n)
  
  x <- importList[[n]]
  # crop data to months observed
  monthCount <- tapply(x$count, x$monthNum, sum)
  monthRange <- as.numeric(c(names(monthCount)[which(monthCount > 0)[1]],
                  names(monthCount)[rev(which(monthCount > 0))[1]]))
  
  if(sum(is.na(monthRange))==2 | diff(monthRange) < 10){return(NULL)}
  
  x <- x[x$monthNum >= monthRange[1] &
         x$monthNum <= monthRange[2],]              
  
  #tempGam <- gam(cbind(count, fail) ~ s(monthNum), data=x, family=binomial, method="REML")
  tempGamm <- gamm4(cbind(count, fail) ~ s(monthNum), random = ~(1|userId),
                  data=x, family=binomial)
  
  pred.df <- data.frame(monthNum=seq(monthRange[1], monthRange[2], 1))
  
  return(cbind(pred.df,
               as.data.frame(predict(tempGamm$gam, newdata=pred.df, se.fit=TRUE))))
})

sapply(1:length(importMods), function(n){
  
  print(n)
  subDf <- importMods[[n]]
  
  if(is.null(subDf)){return(NULL)}
  
  changeRatio <- mean(plogis(subDf$fit)[nrow(subDf) - (0:9)])/mean(plogis(subDf$fit)[1:10])
  
  funCat <- funSub$functionCat[funSub$function. == importFuns[n]][1]
  tempCol <- c("#2369bd", "grey75", "#d12239")[as.numeric(funCat)] 
  backCol <- col2rgb(c("#2369bd", "grey75", "#d12239")[as.numeric(funCat)])/255
  lineStyle <- "solid"
  
  if(is.na(tempCol[1])){tempCol=c(0,0,0)}
  
  polygon(x=c(subDf$monthNum, rev(subDf$monthNum)),
          y=log(plogis(c(subDf$fit + 1.96 * subDf$se.fit,
                         rev(subDf$fit - 1.96 * subDf$se.fit)))),
          border=NA, col=rgb(backCol[1], backCol[2], backCol[3], 0.35))
  
  lines(log(plogis(subDf$fit)) ~ subDf$monthNum, 
        col=tempCol, lwd=2)
  text(x=subDf$monthNum[1],
       y=log(plogis(subDf$fit)[1]), cex=0.8,
       labels=importFuns[n], col=tempCol, pos=2)
  text(x=rev(subDf$monthNum)[1],
       y=log(plogis(rev(subDf$fit)[1])),
       labels=paste0(round(changeRatio,2),"x"),
       col=tempCol, pos=4, cex=0.8)
})

box()
text(x=relative.axis.point(0.025,"x"),
     y=relative.axis.point(0.96, "y"),
     labels=c("(B) Function occurrence probability over time"), font=2, adj=0)

close.screen(2)

screen(3)
par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)

total <- tapply(funSub$count,
                list(funSub$monthsSinceJan10Created,
                     funSub$functionCat),
                sum, na.rm=TRUE)
total[is.na(total)] = 0
tempDf <- prop.table(total, 1)

cumulDf <- apply(tempDf, 1, cumsum)
cumulDf <- rbind(rep(0, ncol(cumulDf)), cumulDf)

plot(x=NULL, y=NULL, xlim=c(24,165), ylim=c(0,1),
     axes=FALSE, xlab="", ylab="", yaxs="i")

axis(side=2, at=c(0,0.5,1), line=-4, lty="31")

sapply(2:4, function(n){
  
  polygon(x=c(monthNum, rev(monthNum)),
          y=c(cumulDf[n,], rev(cumulDf[n-1,])),
          col=c("#2369bd", "grey75", "#d12239")[n-1])
  
})

sapply(seq(0,1,0.25), function(x){
segments(x0=min(monthNum), x1=max(monthNum),
         y0=x, y1=x,col=rgb(0,0,0,0.5))
})
rect(xleft=min(monthNum), xright=max(monthNum),
     ybottom=par("usr")[3], ytop=par("usr")[4],
     col=rgb(0,0,0,0))

close.screen(3)

close.screen(all.screens = TRUE)

dev.off()

return(list(allData = importDf,
            funData = importList,
            allModel = tempGamm,
            funModels = importMods))

}