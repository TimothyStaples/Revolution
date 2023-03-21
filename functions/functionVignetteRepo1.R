functionVignetteRepo1 <- function(importFuns, 
                             funSub, 
                             plotName, 
                             ylims, mainYlims){

funSub$monthsSinceJan10CreatedFact <- as.factor(funSub$monthsSinceJan10Created)
funSub$urIDF <- as.factor(funSub$urID)

repoCount <- as.data.frame(table(funSub$monthsSinceJan10CreatedFact[!duplicated(funSub$urID)]))

funSub <- funSub[funSub$function. %in% importFuns,]

importDf <- as.data.frame(table(funSub$monthsSinceJan10CreatedFact[!duplicated(funSub$urID)]))

importDf <- cbind(importDf, repoCount[,2])
colnames(importDf) = c("monthNum", "count", "total")
importDf$fail <- importDf$total - importDf$count
importDf$monthNum <- as.numeric(as.character(importDf$monthNum))

importList <- lapply(importFuns, function(x){
  print(x)

  subDf <- funSub[funSub$function.==x,]
  subTable <- as.data.frame(table(subDf$monthsSinceJan10CreatedFact))
  
  subTable <- cbind(subTable, repoCount[,2])
  colnames(subTable) = c("monthNum", "count", "total")
  subTable$fail <- subTable$total - subTable$count
  subTable$monthNum <- as.numeric(as.character(subTable$monthNum))
  
  return(subTable)
})
  
  #                  Models over time & plot ####

monthTotals <- table(funSub$monthsSinceJan10Created[!duplicated(funSub$urID)])
monthNum <- as.numeric(names(monthTotals))

plot(x=NULL, y=NULL, xlim=c(24,165), ylim=ylims,
     axes=FALSE, xlab="", ylab="")

axis(side=2, at=log(c(1, 0.1,0.01,1e-3,1e-4)),
     labels=c(1, 0.1,0.01,1e-3,1e-4))
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
  # monthCount <- tapply(x$count, x$monthNum, sum)
  monthRange <- range(x$monthNum)
  # monthRange <- as.numeric(c(names(monthCount)[which(monthCount > 0)[1]],
  #                 names(monthCount)[rev(which(monthCount > 0))[1]]))
  # 
  # if(sum(is.na(monthRange))==2 | diff(monthRange) < 10){return(NULL)}
  # 
  # x <- x[x$monthNum >= monthRange[1] &
  #        x$monthNum <= monthRange[2],]              
  
  tempGam <- gam(cbind(count, fail) ~ s(monthNum), data=x, family=binomial, method="REML")
  
  pred.df <- data.frame(monthNum=seq(monthRange[1], monthRange[2], 1))
  
  return(cbind(pred.df,
               as.data.frame(predict(tempGam, newdata=pred.df, se.fit=TRUE))))
})

sapply(1:length(importMods), function(n){
  
  print(n)
  subDf <- importMods[[n]]
  
  if(is.null(subDf)){return(NULL)}
  
  changeRatio <- mean(plogis(subDf$fit)[nrow(subDf) - (0:11)])/mean(plogis(subDf$fit)[1:12])

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

return(list(allData = importDf,
            funData = importList,
            funModels = importMods))

}