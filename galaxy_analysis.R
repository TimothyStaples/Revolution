# ########################################################## ####
# Diversification and change of the R programming language   ####
# Author:    Timothy L Staples                               ####
# ########################################################## ####

# Script purpose ####

# Conduct all statistical analyses and produce figures.

# This script makes extensive use of R-Studio's code folding syntax.
# Alt + x to collapse all folds on Windows/Linux: Cmd + Alt + x on Mac.

# SET DETAILS ####

rm(list=ls())
setwd("PATH TO THIS FILE")
sapply(paste0("./functions/", list.files("./functions")), source)

# this function will install packages if not installed, or load them to the library
# if they are
package.loader(c("parallel", "vegan", "mgcv", "lme4", "shape", "DHARMa", "performance", "emmeans",
                 "viridisLite", "hilldiv"))

# DATA PREP ####

# This file is available for download via Dryad (https://doi.org/10.5061/dryad.h18931zrg)
funTableSub <- read.csv("./outputs/commonFunctionLong.csv")

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

funTableRepo$funAsCat <- as.factor(funTableRepo$function.)

funMonthProp <- tapply(funTableRepo$count,
                       list(funTableRepo$monthsSinceJan10Created,
                            funTableRepo$funAsCat),
                       sum, na.rm=TRUE)
funMonthProp[is.na(funMonthProp)] = 0
funMonthProp <- prop.table(funMonthProp, 1)


# This is Table S3
write.csv(cbind(tapply(funTableRepo$count,
                       funTableRepo$monthsSinceJan10Created, sum),
                funMonthProp), 
          "./outputs/monthlyFunctionProp.csv")

# use proportions for actual analyses (counts are simpler for supp table data)
funMonthProp <- tapply(funTableRepo$prop,
                       list(funTableRepo$monthsSinceJan10Created,
                            funTableRepo$funAsCat),
                       sum, na.rm=TRUE)
funMonthProp[is.na(funMonthProp)] = 0
funMonthProp <- prop.table(funMonthProp, 1)

packMonthProp <- tapply(funTableRepo$count,
                       list(funTableRepo$monthsSinceJan10Created,
                            as.factor(funTableRepo$package)),
                       sum, na.rm=TRUE)
packMonthProp[is.na(packMonthProp)] = 0

packMonthProp <- cbind(rowSums(packMonthProp),
                       packMonthProp)
packMonthProp[,1] == funMonthProp[,1]

# This is Table S4
write.csv(packMonthProp, "./outputs/monthlyPackageProp.csv")

# function category proportion data
funCatProp <- list(funMonthProp[,colnames(funMonthProp) %in% 
                                  unique(funTableSub$function.[funTableSub$functionCat=="base"])],
                   funMonthProp[,colnames(funMonthProp) %in% 
                                  unique(funTableSub$function.[funTableSub$functionCat=="other"])],
                   funMonthProp[,colnames(funMonthProp) %in% 
                                  unique(funTableSub$function.[funTableSub$functionCat=="tidyverse"])])

# Function diversity (Models for Fig 1 and Fig S2) ####
#               Gamma diversity ####

# month by month hill numbers
funHill <- do.call("rbind", apply(funMonthProp, 1, function(x){
  data.frame(h0 = hill_div(as.vector(x[x>0]), 0),
             h1 = hill_div(as.vector(x[x>0]), 1),
             h2 = hill_div(as.vector(x[x>0]), 2))
}))
funHill$month <- as.numeric(rownames(funMonthProp))
funHill$monthCount <- table(funTableSub$monthsSinceJan10Created[!duplicated(funTableSub$urID)])

# model over time
h0gam <- gam(h0 ~ s(month, bs="tp") + log(monthCount), data=funHill, family=nb)
h1gam <- gam(log(h1) ~ s(month, bs="tp") + log(monthCount), data=funHill)

#               Per-repo alpha diversity ####

funOccurS <- tapply(funTableRepo$prop,
                    funTableRepo$urID,
                    function(x){
                      
                      c(h0 = hill_div(x, 0),
                        h1 = hill_div(x, 1),
                        h2 = hill_div(x, 2))
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
h1gamL <- gam(log(h1) ~ s(monthsSinceJan10Created, bs = "tp") + log(repoCount), data=funOccurS)

#               Save models ####

saveRDS(list(list(funHill, h0gam, h1gam),
             list(funOccurS, h0gamL, h1gamL)),
        "./outputs/hillDiversityModels.rds")

# Category diversity ####
#               Prop data ####

# month by month hill numbers
funHillCat <- do.call("rbind", lapply(1:length(funCatProp), function(n){
  
  y<-funCatProp[[n]]
  
  temp <- do.call("rbind", apply(y, 1, function(x){
    data.frame(h0 = hill_div(as.vector(x), 0),
               h1 = hill_div(as.vector(x), 1),
               h2 = hill_div(as.vector(x), 2))
  }))
  
  temp$month <- as.numeric(rownames(y))# wide form data 
  temp$functionCat <- levels(funTableRepo$functionCat)[n]
  return(temp)
}))
monthCount <- as.data.frame(table(funTableSub$monthsSinceJan10Created[!duplicated(funTableSub$urID)]))
colnames(monthCount) <- c("month", "monthCount")
funHillCat <- merge(funHillCat, monthCount,
                    by.x="month", by.y="month", all.x=TRUE, all.y=FALSE, sort=FALSE)

# per repo diversity

funOccurSCat <- tapply(funTableRepo$prop,
                       list(as.factor(funTableRepo$urID),
                            funTableRepo$functionCat),
                       function(x){
                         
                         c(h0 = hill_div(x[x>0], 0),
                           h1 = hill_div(x[x>0], 1),
                           h2 = hill_div(x[x>0], 2))
                         
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
h0CatGam <- gam(h0 ~ s(month, bs="tp", by=functionCat) + functionCat + log(monthCount), data=funHillCat, family=nb)
h1CatGam <- gam(log(h1) ~ s(month, bs="tp", by=functionCat) + functionCat + log(monthCount), data=funHillCat)

#               Alpha diversity ####

funOccurSCatDf$functionCat <- as.factor(funOccurSCatDf$functionCat)

# Two step process for repo diversity
funOccurSCatDf$h1Bin <- ifelse(funOccurSCatDf$h1 > 0, 1, 0)

h1CatgamBin <- gam(h1Bin ~ s(monthsSinceJan10Created, bs = "tp", by=functionCat) + functionCat + log(repoCount), 
                   data=funOccurSCatDf, family=binomial)

saveRDS(h1CatgamBin, "./outputs/CatAlphaBinaryModel.rds")

funOccurNonZero <- droplevels(funOccurSCatDf[funOccurSCatDf$h1Bin > 0, ])

h1CatgamNonZero <- gam(log(h1) ~ s(monthsSinceJan10Created, bs = "tp", by=functionCat) + functionCat + log(repoCount), 
                       data=funOccurNonZero)

saveRDS(h1CatgamNonZero, "./outputs/CatAlphaGammaModel.rds")

#               Model predictions ####

colnames(funOccurS)[colnames(funOccurS)=="month"] = "monthsSinceJan10Created"
colnames(funOccurS)[colnames(funOccurS)=="scriptCount"] = "repoCount"

monthSeq <- seq(min(funOccurSCatDf$monthsSinceJan10Created), 
                max(funOccurSCatDf$monthsSinceJan10Created), 1)

#                       gamma predictions (overall) ####

pred.df <- data.frame(month = seq(min(funHill$month), max(funHill$month), 1),
                      monthCount = exp(mean(log(funHill$monthCount))))
h0P <- cbind(pred.df, as.data.frame(predict(h0gam, newdata=pred.df, se.fit=TRUE)))
h1P <- cbind(pred.df, as.data.frame(predict(h1gam, newdata=pred.df, se.fit=TRUE)))

#                       gamma predictions (for categories) ####

pred.df <- data.frame(month = rep(monthSeq, 3),
                      monthCount = exp(mean(log(funHillCat$monthCount))),
                      functionCat = rep(levels(as.factor(funOccurSCatDf$functionCat)), each=length(monthSeq)))
hoCatPred <- cbind(pred.df, as.data.frame(predict(h0CatGam, newdata=pred.df, se.fit=TRUE)))
h1CatPred <- cbind(pred.df, as.data.frame(predict(h1CatGam, newdata=pred.df, se.fit=TRUE)))
funOccurNonZero <- droplevels(funOccurSCatDf[funOccurSCatDf$h1Bin > 0, ])

#                       alpha predictions (overall) ####

pred.dfL <- data.frame(monthsSinceJan10Created = seq(min(funOccurS$monthsSinceJan10Created), 
                                                     max(funOccurS$monthsSinceJan10Created), 1),
                       repoCountL = mean(log(funOccurS$repoCount)),
                       repoCount = exp(mean(log(funOccurS$repoCount))))

h0PL <- cbind(pred.dfL, as.data.frame(predict(h0gamL, newdata=pred.dfL, se.fit=TRUE)))
h1PL <- cbind(pred.dfL, as.data.frame(predict(h1gamL, newdata=pred.dfL, se.fit=TRUE)))

colnames(h0PL)[1] = "month"
colnames(h1PL)[1] = "month"

#                       alpha predictions (categories) ####

pred.dfL <- data.frame(monthsSinceJan10Created = rep(monthSeq, 3),
                       repoCountL = mean(log(funOccurNonZero$repoCount)),
                       repoCount = exp(mean(log(funOccurNonZero$repoCount))),
                       functionCat=rep(levels(as.factor(funOccurSCatDf$functionCat)), each=length(monthSeq)))
h1CatBinPL <- cbind(pred.dfL, as.data.frame(predict(h1CatgamBin, newdata=pred.dfL, se.fit=TRUE)))
h1CatGammaPL <- cbind(pred.dfL, as.data.frame(predict(h1CatgamNonZero, newdata=pred.dfL, se.fit=TRUE)))
colnames(h1CatBinPL)[1] = "month"
colnames(h1CatBinPL)[1] = "month"

#               H0 plot (Fig 1)####

colnames(funOccurS)[colnames(funOccurS) == "monthsSinceJan10Created"] = "month"

pdf(date.wrap("./plots/function_diversityRevised",".pdf"), 
    height=4.5, width=8, useDingbats=FALSE)

par(mfcol=c(2,3), mar=c(0,3,0,1), oma=c(2,1,2,0), ps=10, tcl=-0.25)

ylims <- list(c(600,1250), c(50,600), 
              c(20,22), c(0,1.1))

sapply(1:4, function(n){
  
  print(n)
  temp.pred <- list(h0P, hoCatPred, h0PL, h1CatBinPL)[[n]]
  temp.raw <- list(funHill[,c("h0", "month")], 
                   funHillCat[,c("h0", "month", "functionCat")], 
                   funOccurS[,c("h0", "month")], 
                   funOccurSCatDf[,c("h1", "monthsSinceJan10Created", "functionCat")]
  )[[n]]
  
  if(n == 3){
    temp.raw <- cbind(tapply(temp.raw[,1], temp.raw[,2], mean),
                      sort(unique(temp.raw[,2])))
  }
  
  if(n %in% c(2,4)){
    colnames(temp.pred)[colnames(temp.pred)=="functionCat"] = "cat"
    colnames(temp.raw)[colnames(temp.raw)=="functionCat"] = "cat"
    
    colnames(temp.raw)[2] = "month"
  }
  
  plot(x=NULL, y=NULL, xlim=range(funHill$month) + c(0, 20), ylim=ylims[[n]], axes=FALSE,
       ylab="", xlab="")
  
  if(n ==1){
    axis(side=2, las=1, mgp=c(3,0.5,0),
         at=seq(500,2000,100),
         labels=format(seq(500,2000,100), big.mark=","))
  } else {
    axis(side=2, las=1, mgp=c(3,0.5,0))
  }
  
  mtext(side=2, line=ifelse(n %in% 1:2, 2.75, 2.5), 
        text=c("Global function diversity",
               "Global function diversity",
               "Within-repository diversity",
               "Within-respoitory probability")[n], las=0, cex=0.8)
  
  axis(side=1, at=seq(37,145,12), labels=NA)
  
  if(n %in% c(2,4)){
    par(xpd=NA)
    text(x=seq(49,145,12),
         y=relative.axis.point(-0.04,"y"),
         labels=2014:2022, srt=30, adj=1)
    par(xpd=FALSE)
    
  }
  
  #### Plot per-category data (plots 2 and 4)
  if(n %in% c(2,4)){
    
    predSplit <- split(temp.pred, f=temp.pred$cat)
    
    # if(n==2){
    rawSplit <- split(temp.raw, f=temp.raw$cat)
    #}
    
    sapply(1:3, function(n1){
      
      print(n1)
      subPred <- predSplit[[n1]]
      subRaw <- rawSplit[[n1]]
      
      pointCol <- col2rgb(c("#2369bd", "grey75", "#d12239")[n1])/255
      faintCol <- col2rgb(c("#8eb8ea", "grey90", "#ef9da7")[n1])/255
      linecol <- col2rgb(c("black", "grey75", "#d12239")[n1])/255
      
      if(n==2){
        segments(x0=subRaw[-nrow(subRaw),2], x1=subRaw[-1,2],
                 y0=subRaw[-nrow(subRaw), 1], y1=subRaw[-1,1],
                 pch=16, col=rgb(faintCol[1], faintCol[2], faintCol[3], 1), lwd=1)
        
        points(subRaw[,1] ~ subRaw[,2],
               pch=21, 
               bg=rgb(faintCol[1], faintCol[2], faintCol[3], 1),
               col=rgb(linecol[1], linecol[2], linecol[3], 0.5), cex=0.7)
      }
      
      # log link
      if(n == 2){
        
        polygon(x=c(subPred$month, rev(subPred$month)),
                y=exp(c(subPred$fit + 1.96 * subPred$se.fit,
                        rev(subPred$fit - 1.96 * subPred$se.fit))),
                border=NA, col=rgb(pointCol[1], pointCol[2], pointCol[3], 0.5)) 
        
        if(n1 == 3){pointCol <- linecol}
        lines(exp(subPred$fit) ~ subPred$month,
              col=rgb(pointCol[1], pointCol[2], pointCol[3]), lwd=2,
              lty=c("solid", "solid", "solid")[n1])
        
        propChange <- exp(rev(subPred$fit)[1]) / exp(subPred$fit[1])
        
        text(x=rev(subPred$month)[1],
             y=exp(rev(subPred$fit)[1]),
             pos=4, offset=0.3,
             labels=paste0(c("B", "O", "T")[n1], " (", sprintf("%.1f", round(propChange,1)),"x)"),
             font=2,
             col=rgb(pointCol[1], pointCol[2], pointCol[3]))
      }
      
      # logit link
      if(n == 4){
        
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
             pos=4, offset=0.3,
             labels=paste0(c("B", "O", "T")[n1], " (", sprintf("%.1f", round(propChange,1)),"x)"),
             font=2,
             col=rgb(pointCol[1], pointCol[2], pointCol[3]))
        
      }
      
    }) 
  }
  
  #### Plot global data
  if(n %in% c(1,3)){
    yearVect <- temp.raw[,2] %/% 12
    yearVect <- yearVect - min(yearVect) + 1
    yearBase <- viridis(max(yearVect), option="C", end=0.8)
    yearCol <- yearBase[yearVect]
    
    if(n == 1){targets = seq(-0.5, 0.5, 0.1)}
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
    
  }
  
  text(x=relative.axis.point(0.05, "x"),
       y=relative.axis.point(0.95, "y"),
       labels=paste0("(", LETTERS[c(1,3,2,4)][n],")"), font=2)
  
  box()
  
})

dev.off()

#               H1 plot (Fig S2) ####

pdf(date.wrap("./plots/function_diversityRevisedH1",".pdf"), 
    height=4.5, width=8, useDingbats=FALSE)

par(mfcol=c(2,3), mar=c(0,3,0,1), oma=c(2,1,2,0), ps=10, tcl=-0.25)

ylims <- list(c(60,135), c(5,100), 
              c(11.5,14), c(0,15))

sapply(1:4, function(n){
  
  print(n)
  temp.pred <- list(h1P, h1CatPred, h1PL, h1CatGammaPL)[[n]]
  temp.raw <- list(funHill[,c("h1", "month")], 
                   funHillCat[,c("h1", "month", "functionCat")], 
                   funOccurS[,c("h1", "month")], 
                   funOccurSCatDf[,c("h1", "monthsSinceJan10Created", "functionCat")]
  )[[n]]
  
  if(n == 3){
    temp.raw <- cbind(tapply(temp.raw[,1], temp.raw[,2], mean),
                      sort(unique(temp.raw[,2])))
  }
  
  if(n %in% c(2,4)){
    colnames(temp.pred)[colnames(temp.pred)=="functionCat"] = "cat"
    colnames(temp.raw)[colnames(temp.raw)=="functionCat"] = "cat"
    
    colnames(temp.raw)[2] = "month"
  }
  
  plot(x=NULL, y=NULL, xlim=range(funHill$month) + c(0, 20), ylim=ylims[[n]], axes=FALSE,
       ylab="", xlab="")
  
  if(n == 1){
    axis(side=2, las=1, mgp=c(3,0.5,0),
         at=seq(0,200,20),
         labels=seq(0,200,20))
  } else {
    axis(side=2, las=1, mgp=c(3,0.5,0))
  }
  
  mtext(side=2, line=ifelse(n %in% 1:2, 2.75, 2.5), 
        text=c("Global function diversity",
               "Global function diversity",
               "Within-repository diversity",
               "Within-respoitory diversity")[n], las=0, cex=0.8)
  
  axis(side=1, at=seq(37,145,12), labels=NA)
  
  if(n %in% c(2,4)){
    par(xpd=NA)
    text(x=seq(49,145,12),
         y=relative.axis.point(-0.04,"y"),
         labels=2014:2022, srt=30, adj=1)
    par(xpd=FALSE)
    
  }
  
  #### Plot per-category data (plots 2 and 4)
  if(n %in% c(2,4)){
    
    predSplit <- split(temp.pred, f=temp.pred$cat)
    
    # if(n==2){
    rawSplit <- split(temp.raw, f=temp.raw$cat)
    #}
    
    sapply(1:3, function(n1){
      
      print(n1)
      subPred <- predSplit[[n1]]
      subRaw <- rawSplit[[n1]]
      
      pointCol <- col2rgb(c("#2369bd", "grey75", "#d12239")[n1])/255
      faintCol <- col2rgb(c("#8eb8ea", "grey90", "#ef9da7")[n1])/255
      linecol <- col2rgb(c("black", "grey75", "#d12239")[n1])/255
      
      if(n==2){
        segments(x0=subRaw[-nrow(subRaw),2], x1=subRaw[-1,2],
                 y0=subRaw[-nrow(subRaw), 1], y1=subRaw[-1,1],
                 pch=16, col=rgb(faintCol[1], faintCol[2], faintCol[3], 1), lwd=1)
        
        points(subRaw[,1] ~ subRaw[,2],
               pch=21, 
               bg=rgb(faintCol[1], faintCol[2], faintCol[3], 1),
               col=rgb(linecol[1], linecol[2], linecol[3], 0.5), cex=0.7)
      }
      
      
      polygon(x=c(subPred$month, rev(subPred$month)),
              y=exp(c(subPred$fit + 1.96 * subPred$se.fit,
                      rev(subPred$fit - 1.96 * subPred$se.fit))),
              border=NA, col=rgb(pointCol[1], pointCol[2], pointCol[3], 0.5)) 
      
      if(n1 == 3){pointCol <- linecol}
      lines(exp(subPred$fit) ~ subPred$month,
            col=rgb(pointCol[1], pointCol[2], pointCol[3]), lwd=2,
            lty=c("solid", "solid", "solid")[n1])
      
      propChange <- exp(rev(subPred$fit)[1]) / ifelse(n==2, subRaw[4,1], exp(subPred$fit[1]))
      text(x=rev(subPred$month)[1],
           y=exp(rev(subPred$fit)[1]),
           pos=4, offset=0.3,
           labels=paste0(c("B", "O", "T")[n1], " (", sprintf("%.1f", round(propChange,1)),"x)"),
           font=2,
           col=rgb(pointCol[1], pointCol[2], pointCol[3]))
      
    }) 
  }
  
  #### Plot global data
  if(n %in% c(1,3)){
    yearVect <- temp.raw[,2] %/% 12
    yearVect <- yearVect - min(yearVect) + 1
    yearBase <- viridis(max(yearVect), option="C", end=0.8)
    yearCol <- yearBase[yearVect]
    
    if(n == 1){targets = round(seq(-0.6, 1, 0.2),1)}
    if(n == 2){targets = seq(-0.6, 1, 0.1)}
    if(n %in% c(3:4)){targets = seq(-0.1,0.15,0.02)}
    if(n == 5){targets = seq(-0.8,1,0.05)}
    if(n == 6){targets = seq(-0.8,1,0.1)}
    
    sapply(targets, function(x){
      
      # the first non-strange estimate is #4 - use raw values here
      if(n %in% c(1,3)){
        xpos = temp.raw[,1][4] + x * temp.raw[,1][4]
      } else {
        xpos = temp.raw[,1][4] + x * temp.raw[,1][4]
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
    
  }
  
  text(x=relative.axis.point(0.05, "x"),
       y=relative.axis.point(0.95, "y"),
       labels=paste0("(", LETTERS[c(1,3,2,4)][n],")"), font=2)
  
  box()
  
})

dev.off()


# Compositional change over time (Models for Fig 2, Fig S4 and Fig S6) ####
#               nMDS ####

funOrd <- metaMDS(funMonthProp, distance="bray")
plot(funOrd, display="sites", type = "t")
points(funOrd$species, col="red", cex=0.5)

#               time RDA ####

monthVect <- as.numeric(rownames(funMonthProp))

funMonthSqrt <- sqrt(funMonthProp)
funMonthSqrt <- prop.table(funMonthSqrt, margin=1)

funrda <- dbrda(funMonthSqrt ~ monthVect, distance="bray")

#               Function changes over time ####

funOccur <- table(funTableRepo$monthsSinceJan10Created, funTableRepo$function.)

monthCount <- table(funTableRepo$monthsSinceJan10Created[!duplicated(funTableRepo$urID)])
monthCount <- monthCount[match(rownames(funOccur), names(monthCount))]
monthCount <- as.vector(monthCount)
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

funCount <- as.data.frame(table(funTableSub$function.))
funTop <- merge(funTop, funCount,
                by.x="function.", by.y="Var1",
                all.x=TRUE, all.y=FALSE, sort=FALSE)
funTop$prob2021 <- plogis(funTop$`(Intercept)`)

write.csv(funTop, "./outputs/functionTrendTable.csv")

# now do the same at the package level
packTableSub <- funTableRepo[!duplicated(paste0(funTableRepo$urID,":", funTableRepo$package)),]

packOccur <- as.matrix(table(packTableSub$monthsSinceJan10Created, packTableSub$package))

monthCount <- as.vector(table(packTableSub$monthsSinceJan10Created[!duplicated(packTableSub$urID)]))

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

packTop <- merge(packTrends, funTableRepo[!duplicated(funTableRepo$package), 
                                         c("package", "functionCat")],
                 by.x="package", by.y="package", all.x=TRUE, all.y=FALSE, sort=FALSE)

packTop$cat = packTop$functionCat

packCount <- as.data.frame(table(packTableSub$package))
packTop <- merge(packTop, packCount,
                by.x="package", by.y="Var1",
                all.x=TRUE, all.y=FALSE, sort=FALSE)
packTop$prob2021 <- plogis(packTop$`(Intercept)`)

write.csv(packTop, "./outputs/packageTrendTable.csv")

#               Plot (Fig 2 & Fig S4) ####

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
     axes=FALSE, xlab="", ylab="", xlim=c(-0.08,0.1), ylim=log10(c(8e-7, 1)))
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
      text="Probability of function use in December 2021", las=0)
mtext(side=1, line=0.75, text="Slope of change in function use over time")

rect(xleft=-0.01, xright=0.01,
     ybottom=par('usr')[3], ytop=par('usr')[4],
     col="grey90", border=NA)
abline(v=log10(1), lty="31", col="grey60")

funOut <- funTop[funTop$monthVect >= par("usr")[2],]
funIn <- funTop[funTop$monthVect < par("usr")[2],]


# base package points
baseFun <- funIn[funIn$cat=="base",]
baseHull <- chull(baseFun[,c("(Intercept)", "monthVect")])
basecol <- col2rgb("#2369bd")/255

otherFun <- funIn[funIn$cat == "other",]
otherHull <- chull(otherFun[, c("(Intercept)", "monthVect")])
othercol <- col2rgb("grey75")/255

tidyFun <- funIn[funIn$cat == "tidyverse" | funIn$cat == "tidyextend",]
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

legend(x=relative.axis.point(0.65, "x"), y=relative.axis.point(0.2, "y"),
       pch=c(21,21,16), col=c("black","black", "grey75"), pt.bg=c("#2369bd", "#d12239", NA),
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
     axes=FALSE, xlab="", ylab="", xlim=c(-0.08,0.1), ylim=log10(c(8e-7, 1)))

textFuns <- funTop[funTop$function. %in%
                    c("library", "pipeoperator",
                      "read.csv", "read.table", "read_csv", "read.delim", "read_excel",
                                                 "read_xlsx", "read.xlsx", "read.csv2", "read_tsv", "read_delim",
                                                 "read.xls", "fread",
                      "merge", "inner_join", "left_join", "right_join", 
                        "full_join", "semi_join", "anti_join", "merge.data.table",
                        "join", "join_all", "reshape", "pivot_longer", "pivot_wider",
                        "melt", "cast", "gather" ,"spread", "xtabs",
                        "dcast", "acast", "stack", "unstack"),]
                      

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
     axes=FALSE, xlab="", ylab="", xlim=c(-0.05,0.09), ylim=log10(c(8e-5, 1.1)))
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
      text="Probability of package use in December 2021", las=0)
mtext(side=1, line=1, text="Slope of change in package use over time")

rect(xleft=-0.01, xright=0.01,
     ybottom=par('usr')[3], ytop=par('usr')[4],
     col="grey90", border=NA)
abline(v=log10(1), lty="31", col="grey60")

funOut <- packTop[packTop$monthVect >= par("usr")[2],]
funIn <- packTop[packTop$monthVect < par("usr")[2],]

# base package points
baseFun <- funIn[funIn$cat=="base",]
baseHull <- chull(cbind(log10(plogis(baseFun[,"(Intercept)"])), 
                        baseFun[,"monthVect"]))
basecol <- col2rgb("#2369bd")/255

otherFun <- funIn[funIn$cat=="other",]
otherHull <- chull(cbind(log10(plogis(otherFun[,"(Intercept)"])), 
                         otherFun[,"monthVect"]))
othercol <- col2rgb("grey75")/255

tidyFun <- funIn[funIn$cat == "tidyverse" | funIn$cat == "tidyextend",]
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
     axes=FALSE, xlab="", ylab="", xlim=c(-0.05,0.09), ylim=log10(c(8e-5, 1.1)))

textFuns <- packTop
text(x=textFuns$monthVect,
     y=log10(plogis(textFuns$`(Intercept)`)),
     labels=textFuns$package, cex=0.75,
     col=c("#2369bd", "grey75", "#d12239")[as.factor(textFuns$cat)])
close.screen(5)

close.screen(all.screens=TRUE)

dev.off()

#               Summary stats ####

# how many functions used more than 10% in Dec 2021?
table(plogis(funTop$`(Intercept)`) >= 0.1) / nrow(funTop)
table(funTop$function.[plogis(funTop$`(Intercept)`) >= 0.1])

funTop$cat[funTop$cat == "tidyextend"] = "tidyverse"
funTop <- droplevels(funTop)

# % of functions increasing versus decreasing
funTop$sigIncr <- (funTop$monthVect - 1.96 * funTop$monthSE) > 0
funTop$sigDecr <- (funTop$monthVect + 1.96 * funTop$monthSE) < 0

tapply(funTop$monthVect, funTop$cat, summary)

sigChange <- cbind(table(funTop$sigDecr, funTop$cat)[2,],
                   table(!funTop$sigIncr & !funTop$sigDecr, funTop$cat)[2,],
                   table(funTop$sigIncr, funTop$cat)[2,])
colSums(sigChange) / sum(sigChange)

sigProp <- as.data.frame.matrix(sigChange) / t(t(table(funTop$cat)))

meanChange <- table(funTop$cat,
                    cut(funTop$monthVect, breaks=c(-1, -es, es, 1)))
colSums(meanChange) / sum(meanChange)
meanProp <- as.data.frame.matrix(meanChange) / t(t(table(funTop$cat)))                    


# package stats

packTop$cat[packTop$cat == "tidyextend"] = "tidyverse"
packTop <- droplevels(packTop)

packTop$sigIncr <- (packTop$monthVect - 1.96 * packTop$monthSE) > 0
packTop$sigDecr <- (packTop$monthVect + 1.96 * packTop$monthSE) < 0

pSigChange <- cbind(table(packTop$sigDecr, packTop$cat)[2,],
                    table(!packTop$sigIncr & !packTop$sigDecr, packTop$cat)[2,],
                    table(packTop$sigIncr, packTop$cat)[2,])

pSigProp <- as.data.frame.matrix(pSigChange) / t(t(table(packTop$cat)))

pMeanChange <- table(packTop$cat,
                     cut(packTop$monthVect, breaks=c(-1, -es, es, 1)))

pMeanProp <- as.data.frame.matrix(pMeanChange) / t(t(table(packTop$cat)))                    


# Relationship between occurrence probability and slope (Fig. S6) ####
funTop$cat[funTop$cat == "tidyextend"] = "tidyverse"
funTop <- droplevels(funTop)
funTop$cat <- relevel(funTop$cat, "base")

funTop$absVect <- abs(funTop$monthVect)

intSlopelm <- lm(`(Intercept)` ~ -1 + absVect * cat, data=funTop)
plot(simulateResiduals(intSlopelm))
performance(intSlopelm)
summary(intSlopelm)
plot(intSlopelm)
plot(funTop$`(Intercept)` ~ abs(funTop$monthVect) )

emtrends(intSlopelm, pairwise ~ cat, var = "absVect")
emmeans(intSlopelm, pairwise ~ cat,)

pdf("./plots/commonnessStability.pdf", height=4, width=6)

par(mar=c(3,3.5,1,1), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)
plot(funTop$`(Intercept)` ~ funTop$absVect, pch=16, col=c("#2369bd", "gray50", "#d12239")[funTop$cat], 
     cex=0.4, axes=FALSE, xlab="", ylab="")

axis(side=1, mgp=c(3,0.1,0))
mtext(side=1, line=1.25, text="Change in function use probability over time (abs log odds ratio)")

yPoints <- c(1 * 10^(seq(-5,1,1)), 0.5, 0.8)
axis(side=2, at=log(yPoints / (1-yPoints )), labels=yPoints)
mtext(side=2, line=2.5, las=0, text="Probability of usage in December 2021")

#rect(xleft=par("usr")[1], xright=par('usr')[2],
#ybottom=par('usr')[3], ytop=par('usr')[4],
 #    border=NA, col=rgb(1,1,1,0.7))

predData <- data.frame(absVect = rep(seq(min(funTop$absVect), max(funTop$absVect), len=200), 3),
                       cat = rep(levels(funTop$cat), each=200))
corPreds <- cbind(predData,
                  as.data.frame(predict(intSlopelm, newdata=predData, se.fit=TRUE)))

sapply(1:3, function(n){
  
  x <- split(corPreds, f=corPreds$cat)[[n]]
  x <- x[x$absVect >= min(funTop$absVect[funTop$cat == x$cat[1]]) &
         x$absVect <= max(funTop$absVect[funTop$cat == x$cat[1]]),]

  basecol <- col2rgb("#2369bd")/255
  othercol <- col2rgb("grey50")/255
  tidycol <- col2rgb("#d12239")/255
  
  tempCol <- cbind(basecol, othercol, tidycol)[,n]
  
  polygon(x=c(x$absVect, rev(x$absVect)),
          y=c(x$fit + 1.96 * x$se.fit,
              rev(x$fit - 1.96 * x$se.fit)),
          border=NA, col=rgb(tempCol[1], tempCol[2], tempCol[3], 
                             ifelse(n==2, 0.5, 0.3)))
  
  lines(x$fit ~ x$absVect, col=rgb(tempCol[1], tempCol[2], tempCol[3]), lwd=2)
    
})
box()
dev.off()

# Categorical Compositional change (Models for Fig S5) ####

# funCatProp <- lapply(split(funTableRepo, 
#                            f=funTableRepo$functionCat), function(x){
#   
#   tempMonth <- tapply(x$prop,
#                       list(x$monthsSinceJan10Created,
#                            x$funAsCat),
#                       sum, na.rm=TRUE)
#   tempMonth[is.na(tempMonth)] = 0
#   tempMonth <- tempMonth[ ,colSums(tempMonth)>0]
#   tempMonth <- prop.table(tempMonth, 1)
#   return(tempMonth)
#   
# })

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

#               Plot (Fig S5) ####

pdf(date.wrap("./plots/nmdsCat", ".pdf"), height=10, width=3.65)

mdsLims <- list(c(-0.5,0.5),
                c(-1.5, 1.2),
                c(-1.5,0.5))

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

# Are tidyverse models more rich? (Fig S8) ####

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
mtext(side=2, line=1.5, las=0, text="Function richness")

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

# Vinettes of change (Fig 3) ####
#             Import functions ####

ylims <- log(c(1e-5, 1))
mainYlims=c(0,1)

pdf(date.wrap("./plots/vignetteComposite",".pdf"), 
    height=4, width=12, useDingbats = FALSE)

split.screen(rbind(c(0.1,0.4,0.1,0.99),
                   c(0.4,0.7,0.1,0.99),
                   c(0.7,1,0.1,0.99)))

screen(1)
par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)
importFuns <- functionVignetteRepo1(importFuns = c("read.csv", "read.table", "read_csv", "read.delim", "read_excel",
                                 "read_xlsx", "read.xlsx", "read.csv2", "read_tsv", "read_delim",
                                 "read.xls", "fread"),
                 funSub = droplevels(funTableRepo),
                 plotName = "ImportFunctions2",
                 ylims=ylims,
                 mainYlims=mainYlims)
close.screen(1)

screen(2)
par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)
joinFuns <- functionVignetteRepo1(c("merge", "inner_join", "left_join", "right_join", 
                                   "full_join", "semi_join", "anti_join", "merge.data.table",
                                   "join", "join_all"),
                                 funTableRepo,
                                 "JoinFunctions",
                                 ylims=ylims,
                                 mainYlims=mainYlims)
close.screen(2)

screen(3)
par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)
reshapeFuns <- functionVignetteRepo1(c("reshape", "pivot_longer", "pivot_wider",
                                      "melt", "cast", "gather" ,"spread", "xtabs",
                                      "dcast", "acast", "stack", "unstack"),
                                    funTableRepo,
                                    "reshapeFuns",
                                    ylims=ylims,
                                    mainYlims=mainYlims)
close.screen(3)
close.screen(all.screens = TRUE)
dev.off()

saveRDS(importFuns,"./outputs/importFunsFull.rds")


#             Merge functions ####


#             reshape functions #####



#             plyr vs dplyr (Fig S7) ####

modelFuns <- functionVignetteRepo(c("ddply", "ldply", "mapvalues", "llply", "join_all",
                                    "filter", "select", "group_by", "left_join", "bind_rows", "distinct", "funs", "case_when"),
                                  funTableRepo,
                                  "plyrFunctions",
                                  ylims=log(c(1e-6,0.5)),
                                  mainYlims=c(0,1))
