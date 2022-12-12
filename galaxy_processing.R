rm(list=ls())
gc()

library(parallel)
library(tm)
library(vegan)
library(gh)
library(rvest)

library(ecodist)
library(dbscan)
library(meanShiftR)

library(mgcv)
library(lme4)
library(gamm4)
library(MuMIn)
library(usdm)
library(performance)
library(glmmTMB)
library(boot)

library(brms)
library(loo)

library(igraph)
library(qgraph)

library(vioplot)
library(shape)
library(viridisLite)

# ISSUES ####

# 1.  A small handful of days had > 1000 repos created that day, so we missed those

# SET DETAILS ####

setwd("/Users/uqtstapl/Dropbox/Tim/data/RGalaxy")
sapply(paste0("./functions/", list.files("./functions")), source)

#accessAPI for R scripts 
gitToken = c("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",
             "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",
             "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")
              
# FIND GITHUB R SCRIPTS ####
#         Identify repositories with the "R" language flag ####

# get R-based URLs with "R" language flag
# how to get around the 1000 max limit, we need to redo searches for particular periods of time
# let's try each day (from 1st Jan 2010 to 31st Dec 2020)
lapply(1:1000, function(dayN){
  
  print(dayN)
  start <- as.Date("2020-11-30") + dayN-1
  end <- as.Date("2020-11-30") + dayN
  
  if(start >= as.Date("2022-01-1")){
    return(NULL)
  }
  
  mainCall <- paste0("GET /search/repositories?q=language:R+created:",
                     start,"T00:00:00Z..",end, "T00:00:00Z&order=asc?simple=yes&per_page=100")

  tempAPI <- gh(paste0(mainCall), .token = gitToken)
  
  pageN <- ceiling(tempAPI$total_count / 100)
  
  print(paste0(start, " -> ", end, ": ", tempAPI$total_count, " results"))
  
  # if there's no entries for a day, return nothing and wait for a spell to
  # avoid hitting the 30 requests a minute limit
  if(pageN <= 2){Sys.sleep(3)}
  if(pageN == 0){return(NULL)}
  if(pageN > 10){paste0("WARNING...MORE THAN 10 PAGES")}

  RURLs <- lapply(1:ifelse(pageN > 10, 10, pageN), function(n){

    urlString <- paste0(mainCall, "&page=",n)
    tempAPI <- gh(paste0(urlString),
                  .token = gitToken)
    saveRDS(tempAPI, paste0("./gitCalls/", gsub("-","_",start),".",n, ".rds"))
  })
  
})

#                               Get the user and repo name ####

# now we need to read in all of those files, and extract the login and repo name
# a lot of work for two short strings, but hey, labour of love!
RURLS <- as.data.frame(do.call("rbind", lapply(list.files("./gitCalls"), 
                                        function(path){
  print(path)
  x <- readRDS(paste0("./gitCalls/",path))
  
  cbind(id = sapply(x$items, function(y){y$owner$id}),
        login = sapply(x$items, function(y){y$owner$login}),
        repo = sapply(x$items, function(y){y$name}),
        dateCreated = sapply(x$items, function(y){y$created_at}),
        dateUpdated = sapply(x$items, function(y){y$updated_at}))
  
})),
                      stringsAsFactors=FALSE)

# set unique Repo ID for each user
RURLS <- RURLS[order(as.numeric(RURLS$id)),]

RURLS$repoId <- unlist(sapply(table(as.numeric(RURLS$id)), function(x){1:x}))

# convert dates
RURLS$yearCreated <- as.numeric(substr(RURLS$dateCreated, 1, 4))
RURLS$yearUpdated <- as.numeric(substr(RURLS$dateUpdated, 1, 4))

#              get file path for each .R file in each repo ####

# now get R file paths for each URL. Given the sheer number of repos we have to sift
# through, I've employed multiple github OAuth tokens to parallelize this process.
# If you single-thread it, it will take a long time (like months long)

# split our RURLs up into number of segments equal to our gitTokens 
  
RURLSstore <- RURLS

  # first randomize the order so it's not in date created order anymore
  RURLS$loginFile <- gsub("[[:punct:]]", "_",RURLS$login)
  RURLS$loginFile <- gsub("-", "_", RURLS$loginFile, fixed=TRUE)
  RURLS$loginFile <- gsub(".", "_",RURLS$loginFile, fixed=TRUE)
  
  RURLS$repoFile <- gsub("[[:punct:]]", "_",RURLS$repo)
  RURLS$repoFile <- gsub("-", "_", RURLS$repoFile, fixed=TRUE)
  RURLS$repoFile <- gsub(".", "_",RURLS$repoFile, fixed=TRUE)
  
  # remove repos that have been done
  doneFiles <- c(list.files("./gitURLs"),
                 list.files("./gitURLs1"),
                 list.files("./gitGone"))
  
  # we made an error with swapping the "-" for "_" in saving some files)
  repoDoneOld <- paste0(RURLSstore$login, "-",
                        RURLSstore$repo, ".csv")
  repoDone <- paste0(RURLS$loginFile, "-",
                     RURLS$repoFile, ".csv")
  
  RURLStogo <- RURLS[!repoDone %in% doneFiles &
                     !repoDoneOld %in% doneFiles, ]
  
  # remove duplicates
  RURLStogo <- RURLStogo[!duplicated(paste0(RURLStogo$loginFile, "-",RURLStogo$repoFile)),]
  
  RURLsShuffled <- RURLStogo[sample(1:nrow(RURLStogo), nrow(RURLStogo)),]
  
  if(length(gitToken)>1){
  RURLSlist <- split(RURLsShuffled, f=cut(1:nrow(RURLsShuffled), breaks=length(gitToken)))
  } else {
  RURLSlist <- list(RURLsShuffled)
  }
  
  # now use all tokens simultaneously, and hope we don't get rejected by github
  # (the harshest rejection of all)
  no_cores <- length(gitToken)
  cl <- makeCluster(no_cores)
  setDefaultCluster(cl)
  # export data and load libraries in cluster
  clusterExport(varlist=c("RURLSlist", "gitToken", "doneFiles"), envir=environment())
  clusterCall(cl, "library", "gh", character.only = TRUE)
  options(na.action = "na.fail")
  files <- parLapply(cl, 1:length(gitToken), function(tokenN){
    
    # cycle through
    RURLsub <- RURLSlist[[tokenN]]
    tempToken <- gitToken[tokenN]
    
    sapply(1:nrow(RURLsub), function(n){
      
      temp <- unlist(RURLsub[n,])
      
      if(paste0(temp[2], "-", temp[3], ".csv") %in% doneFiles){
        return(NULL)
      }
      
      # check to see if we have this already (in the case that this apply fails at any point)
      tempFileString <- paste0("./gitURLs1/", temp[2], "-", temp[3], ".csv")
      
      # get repo content tree
      repoURL <- paste0("GET /repos/", temp[2],"/", temp[3], "/git/trees/master?recursive=1")
      repoContents <- try(gh(repoURL, .token = tempToken))
      
      if(class(repoContents)[1]=="try-error"){
        writeChar("A", paste0("./gitGone/",temp[2], "-", temp[3], ".csv"))
        return(NULL)
      }
      
      # get url paths to every file in the repo
      filePaths <- paste0("https://raw.githubusercontent.com/",
                          temp[2],"/",temp[3],"/master/",
                          sapply(repoContents$tree, function(x){x$path}))
      
      # get just .R files (ignore everything else)
      repoExts <- substr(filePaths,
                         sapply(filePaths, function(x){rev(gregexpr("\\.", x)[[1]])[1]}),
                         nchar(filePaths))
      
      repoRs <- filePaths[tolower(repoExts) == ".r"]
      
      print(paste0(temp[2], " has *", sum(tolower(repoExts) == ".r"), "* R scripts"))
      
      write.csv(repoRs, tempFileString, row.names=FALSE)
      
      Sys.sleep(1) # to avoid hitting our 30 API requests a minute limit
      
      return(repoRs)
      
    })
    
  })
  stopCluster(cl=cl)
  
# PROCESS FUNCTION DATA ####  
#             Download .R files ####

# Get all the files in each repo
cl <- makeCluster(6)
setDefaultCluster(cl)
options(na.action = "na.fail")
fileDf <- parLapply(cl, list.files("./gitURLs"), function(path){
print(path)
unlist(read.csv(paste0("./gitURLs/",path), stringsAsFactors = FALSE))
})
stopCluster(cl=cl)

# If we are being lazy and trying to source URLs and functions at the same
# time, this will make sure the names match the URLs in the files
names(fileDf) = list.files("./gitURLs")#[1:length(fileDf)]

# create a dummy .csv file column to match repo data

# we made an error with swapping the "-" for "_" in saving some files)
repoDoneOld <- paste0(RURLSstore$login, "-",
                      RURLSstore$repo, ".csv")
repoDone <- paste0(RURLS$loginFile, "-",
                   RURLS$repoFile, ".csv")

RURLS$fileName = paste0(RURLS$loginFile, "-", RURLS$repoFile, ".csv")
RURLS$oldFileName <- paste0(RURLSstore$login, "-", RURLSstore$repo, ".csv")

# now we remove these from our files to do
fileNames <- gsub("\\.csv", "", names(fileDf))
fileIds <- strsplit(fileNames, "-")
fileIds <- data.frame(fullString = fileNames,
                      userName = sapply(fileIds, function(x){x[[1]]}),
                      repoName = sapply(fileIds, function(x){x[[2]]}))

fileIds$userID <- RURLS$id[match(fileIds$userName, RURLS$login)]

fileIds$repoID <- RURLS$repoId[match(gsub("_", "-", fileIds$repoName), 
                                     gsub("_", "-", RURLS$repo))]

fileComplete <- fileIds[!is.na(fileIds$userID) & !is.na(fileIds$repoID),]
fileRepoOnly <- fileIds[!is.na(fileIds$userID) & is.na(fileIds$repoID),]
fileNoUser <- fileIds[is.na(fileIds$userID),]

# some userIds are not available because of poor formatting on my part in 
# the URLs (user names or repos with "-" in them)
uniqueLogins <- unique(RURLS$login)
potentialCombs <- gregexpr("-", gsub("_", "-", fileNoUser$fullString))
logMatch <- sapply(1:length(potentialCombs), function(n){
  
  print(n)
  tempString <- gsub("_","-",fileNoUser$fullString[n])
  tempCombs <- c(potentialCombs[[n]], nchar(tempString))
  potNames <- sapply(tempCombs, function(x){
    substr(tempString, 1, x-1)
  })

  potMatches <- match(potNames, uniqueLogins)
  
  # return the most complete user name
  potMatches <- potMatches[!is.na(potMatches)]
  rev(potMatches)[1]
  
})

fileNoUser$userName = uniqueLogins[logMatch]
fileNoUser$userID = RURLS$id[match(uniqueLogins[logMatch], RURLS$login)]

# now we add these to the ones missing repos
fileFindRepo <- rbind(fileRepoOnly, fileNoUser)
repoByUser <- split(RURLS$repo, f=RURLS$id)
repoMatch <- sapply(1:nrow(fileFindRepo), function(n){
  
print(n)
  
  tempUsers <- fileFindRepo$userID[n]
  
  # targetRepos
  targRepos <- repoByUser[[tempUsers]]
  tempFile <- fileFindRepo$fullString[n]
  
  # remove userId from string
  tempMatch <-sapply(targRepos, function(x){grepl(gsub("_|-|\\.","",x), 
                                                  gsub("_|-|\\.","",tempFile))})
  tempMatch <- names(tempMatch)[tempMatch]
  
  # if no matching repo, return NA
  if(length(tempMatch) == 0){tempMatch=NA}
  
  # if multiple matches, return the longest matching repo
  if(length(tempMatch) > 1){tempMatch = tempMatch[which.max(nchar(tempMatch))]}
  
  return(tempMatch)
})

fileFindRepo$repoName = repoMatch
fileFindRepo$repoID <- RURLS$repoId[match(paste0(fileFindRepo$userID, "-", fileFindRepo$repoName),
                                          paste0(RURLS$id, "-",RURLS$repo))]

# now we have a record of user and repo Ids we can cross-check against the 
# files we've done
fileComb <- rbind(fileComplete,
                  fileFindRepo)

# reorder fileComb to match fileDf using rowname
fileComb <- fileComb[order(as.numeric(rownames(fileComb))),]

# Download and write scripts as raw text files (subset the script-id off)
doneScripts <- c(list.files("./gitScripts"),
                 list.files("./gitScripts1"))

# convert IDs of doneScripts into repoIds that align with all repos to be done
doneCut <- gregexpr("-", doneScripts)
doneCut <- sapply(doneCut, function(x){x[2]})
doneRepo <- unique(substr(doneScripts, 1, doneCut-1))

# add in repositories that don't have any R files despite the R flag
blankRepo <- list.files("./gitBlanks")
blankRepo <- unique(substr(blankRepo, 1, nchar(blankRepo)-4))

cutRepo <- unique(c(doneRepo, blankRepo))

# remaining files are those with user and repo ids NOT in cutRepo
fullId <- paste0(fileComb$userID,"-",fileComb$repoID)
head(fullId)

fileToDo <- fileDf[!fullId %in% cutRepo]

fileToDo <- fileToDo[sample(1:length(fileToDo))]

cl <- makeCluster(10)
setDefaultCluster(cl)
# export data and load libraries in cluster
clusterExport(varlist=c("RURLS", "repoDoneOld", "fileToDo", "doneScripts"), envir=environment())
clusterCall(cl, "library", "gh", character.only = TRUE)
options(na.action = "na.fail")
parLapply(cl,1:length(fileToDo), function(n){
  
  print(n)
  repoFiles <- fileToDo[[n]]
  
  # get ID number of repo
  repoName = names(fileToDo)[n]
  userID = RURLS$id[match(repoName, RURLS$fileName)]
  repoId = RURLS$repoId[match(repoName, RURLS$fileName)]
  
  if(sum(is.na(c(userID, repoId)))>0){
    userID = RURLS$id[match(repoName, repoDoneOld)]
    repoId = RURLS$repoId[match(repoName, repoDoneOld)]
  }
  
  # for repos with R flag but no .R files.
  if(length(repoFiles)==0){
    try(writeChar("A", paste0("./gitBlanks/",userID,"-",repoId, ".txt")))
    return(NULL)
    
    }
  
  # make sure repo hasn't been done already
if(paste0(userID, "-", repoId, "-1.txt") %in% doneScripts){return(NULL)}
  
  sapply(1:length(repoFiles), function(n1){
  #print(n1)
  x <- repoFiles[n1]
    
  fileName <- substr(x,
                     sapply(x, function(x){rev(gregexpr("/", x, fixed=TRUE)[[1]])[1]})+1,
                     nchar(x)-2)
    
  # swap out spaces for 'html spaces' so we can find the file
  filePath <- gsub(" ", "%20", x)
  
  sc <- try(readChar(filePath, nchar=1e6, useBytes=TRUE))
  
  if(class(sc) == "try-error" | length(nchar(sc, type="bytes")) == 0){
    sc = ""
  }
  
  scUTF <- enc2utf8(sc)
  
  try(writeChar(scUTF, paste0("./gitScripts1/",userID,"-",repoId,"-",n1, ".txt")))
  
  })
  
})
stopCluster(cl=cl)

#             Process scripts to extract function calls ####

doneFuns <- c(list.files("./gitFuns"), list.files("./gitFuns1"))
doneFuns <- substr(doneFuns, 1, nchar(doneFuns)-4)

allScripts <- list.files("./gitScripts")
allScriptsSub <- substr(allScripts, 1, nchar(allScripts)-4)

funsToBeDone <- allScripts[!allScriptsSub %in% doneFuns]

cl <- makeCluster(10)
setDefaultCluster(cl)
# export data and load libraries in cluster
clusterExport(varlist=c("funsToBeDone"), envir=environment())
clusterCall(cl, "library", "gh", character.only = TRUE)
options(na.action = "na.fail")
parLapply(cl, funsToBeDone, function(path){

# read in script, allocating 1MB of space each (should be far more than any individual script needs)
sc <- try(readChar(paste0("./gitScripts/", path), nchar=1e6, useBytes=TRUE))
                             
# Switch %>% from dplyr and others to a 'function term' we can use later on
sc <- gsub("%>%", "pipeoperator(", sc, fixed=TRUE, useBytes = TRUE)
sc <- gsub("%in%", "inoperator(", sc, fixed=TRUE, useBytes = TRUE)

# break script up into separate code lines
scLine <- unlist(strsplit(sc, "\n", fixed=TRUE, useBytes=TRUE))

# remove anything to the right of a comment "#", as this is not "live" code
# and we identify functions via parentheses later
commPos <- regexpr("#", scLine, fixed=TRUE, useBytes=TRUE)
commPos <- ifelse(commPos == -1, nchar(scLine)+1, commPos)

scLive <- substr(scLine, 1, commPos-1)

# remove any empty lines
scLive <- scLive[nchar(scLive) > 0]

# paste back into a single script for function extraction
scLive <- paste0(scLive, collapse=" ")

# split our script into separate lines based on the location of left parentheses
scLive <- unlist(strsplit(scLive, "(", fixed=TRUE, useBytes=TRUE))

# replace all punctuation except periods and underscores with spaces
scLive <- gsub("(?!.)(?!_)[[:punct:]]"," ", scLive, perl=TRUE)
                             
# remove a couple of troublesome punctuation marks and line breaks
scLive <-  gsub("\n", " ", scLive, fixed=TRUE, useBytes = TRUE)
scLive <-  gsub("[", " ", scLive, fixed=TRUE, useBytes = TRUE)
scLive <-  gsub("]", " ", scLive, fixed=TRUE, useBytes = TRUE)
scLive <-  gsub(")", " ", scLive, fixed=TRUE, useBytes = TRUE)
scLive <-  gsub(":", " ", scLive, fixed=TRUE, useBytes = TRUE)
scLive <-  gsub("+", " ", scLive, fixed=TRUE, useBytes = TRUE)
scLive <-  gsub(",", " ", scLive, fixed=TRUE, useBytes = TRUE)
scLive <-  gsub("-", " ", scLive, fixed=TRUE, useBytes = TRUE)
scLive <-  gsub('"', " ", scLive, fixed=TRUE, useBytes = TRUE)
scLive <-  gsub("!", " ", scLive, fixed=TRUE, useBytes = TRUE)
scLive <-  gsub("=", " ", scLive, fixed=TRUE, useBytes = TRUE)
scLive <-  gsub("<", " ", scLive, fixed=TRUE, useBytes = TRUE)
scLive <-  gsub(">", " ", scLive, fixed=TRUE, useBytes = TRUE)
scLive <-  gsub("*", " ", scLive, fixed=TRUE, useBytes = TRUE)
scLive <-  gsub("$", " ", scLive, fixed=TRUE, useBytes = TRUE)
scLive <-  gsub("%", " ", scLive, fixed=TRUE, useBytes = TRUE)
scLive <-  gsub("#", " ", scLive, fixed=TRUE, useBytes = TRUE)
scLive <-  gsub("'", " ", scLive, fixed=TRUE, useBytes = TRUE)
scLive <-  gsub("/", " ", scLive, fixed=TRUE, useBytes = TRUE)
scLive <-  gsub("&", " ", scLive, fixed=TRUE, useBytes = TRUE)
scLive <-  gsub("{", " ", scLive, fixed=TRUE, useBytes = TRUE)
scLive <-  gsub("}", " ", scLive, fixed=TRUE, useBytes = TRUE)
scLive <-  gsub(";", " ", scLive, fixed=TRUE, useBytes = TRUE)
scLive <-  gsub("\t", " ", scLive, fixed=TRUE, useBytes = TRUE)
scLive <-  gsub("\r", " ", scLive, fixed=TRUE, useBytes = TRUE)
scLive <-  gsub("~", " ", scLive, fixed=TRUE, useBytes = TRUE)
scLive <-  gsub("|", " ", scLive, fixed=TRUE, useBytes = TRUE)
scLive <-  gsub("^", " ", scLive, fixed=TRUE, useBytes = TRUE)
scLive <-  gsub("∫", " ", scLive, fixed=TRUE, useBytes = TRUE)
scLive <-  gsub("`", " ", scLive, fixed=TRUE, useBytes = TRUE)
scLive <-  gsub("‘", " ", scLive, fixed=TRUE, useBytes = TRUE)
scLive <-  gsub("?", " ", scLive, fixed=TRUE, useBytes = TRUE)
scLive <-  gsub("\\\\", " ", scLive, fixed=TRUE, useBytes=TRUE)

# make sure there are no leading spaces prior to these identifying punctuation marks
scLive <- gsub(" \\(", "\\(", scLive)
scLive <- gsub(" \\_", "\\_", scLive)
scLive <- gsub(" \\.", "\\.", scLive)

# now the function name should be the final word on each line, so cutting the
# substring between the final space " " and the end of the line should give us
# the function name.

startPos <- sapply(scLive, function(x){rev(gregexpr(" ", x, useBytes = TRUE)[[1]])[1]})
startPos = ifelse(startPos < 0, 0, startPos)

scFun <- try(substr(scLive, startPos, 1000))

# if we have a unicode symbol appear somewhere, the substr function will fail
# and we'll have to remove the offending symbol and try again
if(class(scFun)=="try-error"){
 
  # find which string is the problem 
  subTry <- unlist(sapply(scLive, function(x){
    y = try(substr(x, 1, 1))
    return(class(y)[1])
    }))

  offendingString <- which(subTry == "try-error")
  
  # for now just remove offending strings
  scLive <- scLive[-(offendingString)]
  
  startPos <- sapply(scLive, function(x){rev(gregexpr(" ", x, useBytes = TRUE)[[1]])[1]})
  startPos = ifelse(startPos < 0, 0, startPos)
  scFun <- try(substr(scLive, startPos, 1000))
  

}

# make sure we didn't pick up any extra spaces anywhere
scFun <- gsub(" ", "", scFun)

# exclude empty strings
scFun <- scFun[scFun != ""]

# now paste functions together into a single string so we can use word association
# algorithms to get function correlation scores.

write.csv(scFun, paste0("./gitFuns/", substr(path, 1, nchar(path)-4), ".csv"),
          row.names=FALSE)

})
stopCluster(cl=cl)

# and extract function calls #####

cl <- makeCluster(4)
setDefaultCluster(cl)
options(na.action = "na.fail")
funlist <-  parLapply(cl, list.files("./gitFuns"), function(path){
  x <- as.vector(read.csv(paste0("./gitFuns/",path)))
  if(nrow(x)==0){return(NA)} else {return(x)}
})
stopCluster(cl=cl)
names(funlist) = list.files("./gitFuns")

# Longform script-count data ####

# how many scripts have no functions
funLengths <- sapply(funlist, function(x){
  if(is.na(x)){return(0)}
  nrow(x)
  })
summary(funLengths==0)

# table the functions in each script for a count, bring in repo-level metadata.
funTableList <- mclapply(1:length(funlist), function(n){
  
                      print(n)
                      funtab <- as.data.frame(table(funlist[[n]]))
                      
                      if(nrow(funtab)==0){return(NULL)}
                      colnames(funtab) = c("function", "count")
                      
                      ids <- strsplit(names(funlist)[n], "-")[[1]]
                      
                      funtab$id = ids[1]
                      funtab$repoId = ids[2]
                      funtab$scriptId = substr(ids[3], 1, nchar(ids[3])-4)
                      
                      return(funtab)
}, mc.cores=3)

saveRDS(funTableList, "./outputs/funTableList.rds")
#funTableList <- readRDS("./outputs/funTableList.rds")

# run through and process files in groups of 1000, binding together each little
# script data-frame. We have to do this in bunches because rbind freaks out if
# we give it too many to do at one time.
funCount <- c(seq(1000, length(funTableList), 1000), nrow(funTableList))
lapply(funCount, function(n){

  temp <- do.call("rbind", funTableList[(n-(funCount[1]-1)):n])
  write.csv(temp, paste0("./gitFunTables/",n,".csv"))
  
})

# then we can go back and rbind our bundles of 1000 scripts together.
funTable <- do.call("rbind", 
                    mclapply(list.files("./gitFunTables/"), function(x){
  read.csv(paste0("./gitFunTables/",x))
},  
mc.cores=8))

funTable$id.repo <- paste0(funTable$id, ".", funTable$repoId)
RURLS$id.repo <- paste0(RURLS$id, ".", RURLS$repoId)

# merge in repo-level metadata
funTable <- merge(funTable, RURLS[,c("id.repo", "dateCreated","dateUpdated",
                                     "yearCreated","yearUpdated")],
                  by.x="id.repo", by.y="id.repo", all.x=TRUE, all.y=FALSE,
                  sort=FALSE)

write.csv(funTable, "./outputs/gitFunctionLong.csv")

# FIND FUNCTION PACKAGES ####

#funTable <- read.csv("./outputs/gitFunctionLong.csv")

# what packages am I missing (I can only check packages I have)?
# funSummTemp <- data.frame(Freq = tapply(funTable$count,
#                                         funTable$function.,
#                                         sum, na.rm=TRUE))
funSummTemp <- as.data.frame(table(funTable$function.))
funSummTemp <- funSummTemp[order(funSummTemp$Freq, decreasing=TRUE),]

# remove all currently installed packages
if(FALSE){
  
  # create a list of all installed packages
  ip <- as.data.frame(installed.packages())
  head(ip)
  # if you use MRO, make sure that no packages in this library will be removed
  ip <- subset(ip, !grepl("MRO", ip$LibPath))
  # we don't want to remove base or recommended packages either\
  ip <- ip[!(ip[,"Priority"] %in% c("base", "recommended")),]
  # determine the library where the packages are installed
  path.lib <- unique(ip$LibPath)
  # create a vector with all the names of the packages you want to remove
  pkgs.to.remove <- ip[,1]
  head(pkgs.to.remove)
  # remove the packages
  sapply(pkgs.to.remove, remove.packages, lib = path.lib)
  
}

if(FALSE){
install.packages(c("tidyverse", "testthat", "tinytest", "shiny",
                   "plyr", "raster", "randomForest", "foreach",
                     "gbm", "caret", "here", "doParallel", "sf",
                   "grid", "rgdal", "stopwords", "plotly",
                   "leaflet", "data.table", "shinydashboard",
                   "lme4", "xlsx", "tm", "ROCR", "igraph",
                   "RUnit", "optparse", "sqldf", "DT", "gWidgets2",
                   "cowplot", "car", "phylobase", "pracma",
                   "packrat", "xts", "usethis", "glmnet",
                   "corrplot", "keras", "RODBC", "pacman", "wordcloud",
                   "shinythemes", "doMC", "ggrepel", "pbdTEST",
                   "quantmod", "ggpubr", "janitor","shinyjs", "pheatmap"))

# some packages require specific install pathways
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("IRanges")
BiocManager::install("GenomicRanges")
BiocManager::install("SRAdb")
BiocManager::install("ssize")
BiocManager::install("phyloseq")
}

# compile packages
install.packs <- do.call("rbind", lapply(.libPaths(), function(path){
  x <- installed.packages(path)[,c("Package", "Priority")]
  x[is.na(x[,2]),2] = "z"
  if(length(x)==0){return(NULL)}
  return(x)
  }))

# some functions are listed in multiple packages ("plot") is a classic example
# but hardly the only one. To ensure each function is only listed once per
# package, we will re-order packages to list the "base" packages first, followed
# by optional and recommended packages, and finally by all others.
install.packs <- install.packs[order(install.packs[,2]),]

packsFun <- do.call("rbind", lapply(install.packs[,1], function(x){
  
  print(x)
  
  tempLoad <- try(library(x, character.only=TRUE))
  if(class(tempLoad)=="try-error"){return(NULL)}
  
  tempPacks <- ls(paste0("package:", x))
  if(length(tempPacks)==0){return(NULL)}
    
  data.frame(fun=ls(paste0("package:", x)), package=x)
  
}))

packsFun <- rbind(packsFun,
                  data.frame(fun=c("pipeoperator", "inoperator", "packages", "windows", "filter", "."),
                             package=c("magrittr", "base", "pkgload", "grDevices", "dplyr", "dplyr")))

# now we reduce this to one package per function, assuming that other versions
# of the same-named function are simply iterations of a base function (e.g.,
# "plot", "predict")
packsFun <- packsFun[!duplicated(packsFun$fun),]

funSumm <- base::merge(x=funSummTemp, y=packsFun,
                 by.x="Var1", by.y="fun",
                 all.x=TRUE, all.y=FALSE, sort=FALSE)
missingPacks <- funSumm[is.na(funSumm$package),]
missingPacks <- missingPacks[order(missingPacks$Freq, decreasing=TRUE),]
summary(missingPacks$Freq)
sum(missingPacks$Freq) / sum(funTable$count)
missingPacks[1:50,]

functionTablePack <- merge(x=funTable, y=packsFun,
                           by.x="function.", by.y="fun",
                           all.x=TRUE, all.y=FALSE, sort=FALSE)

functionTablePack <- functionTablePack[,!colnames(functionTablePack) %in% c("X.1", "X", "id")]

write.csv(functionTablePack, "./outputs/gitFunctionLongPackage.csv",
          row.names = FALSE)

# DATA PROCESSING ####

functionTablePack <- read.csv("./outputs/gitFunctionLongPackage.csv")

# reconfigure dates to include months
functionTablePack$monthUpdated <- as.numeric(substr(functionTablePack$dateUpdated, 6, 7))
functionTablePack$monthsSinceJan10 <- (functionTablePack$yearUpdated-2010)*12 + functionTablePack$monthUpdated

functionTablePack$monthCreated <- as.numeric(substr(functionTablePack$dateCreated, 6, 7))
functionTablePack$monthsSinceJan10Created <- (functionTablePack$yearCreated-2010)*12 + functionTablePack$monthCreated

# separate out user Id which got added to a repo id
functionTablePack$userId <- as.numeric(substr(functionTablePack$id.repo, 1, regexpr("\\.", functionTablePack$id.repo)-1))

# add in a script-level ID
functionTablePack$ursID <- paste0(functionTablePack$id, ".", functionTablePack$repoId, ".", functionTablePack$scriptId)
functionTablePack$urID <- paste0(functionTablePack$id, ".", functionTablePack$repoId)

functionTablePack$package[functionTablePack$function.=="filter"] = "dplyr"
functionTablePack$package[functionTablePack$function.=="select"] = "dplyr"
functionTablePack$package[functionTablePack$function.=="."] = "plyr"

# capture tidyverse packages
functionTablePack$tidyverse <- functionTablePack$package %in% c("ggplot2", "dplyr", "tidyr",
                                              "readr", "purrr", "tibble",
                                              "stringr", "forcats", "magrittr", "glue")

functionTablePack$tidyExtend <- functionTablePack$package %in% c("autoplotly", "calendR", "ComplexUpset", "cowplot",
                                               "esquisse", "geomnet", "ggalluvial", "GGally", "ggalt",
                                               "gganimate", "ggasym", "ggbeesswarm", "ggcharts",
                                               "ggdag", "ggdark", "ggdist", "ggedit", "ggExtra",
                                               "ggfittext", "ggforce", "ggformula", "ggfortify",
                                               "gggenes", "ggh4x", "gghalves", "ggip", "ggiraph",
                                               "gglm", "gglorenz", "ggmosaic", "ggmuller", "ggnetwork",
                                               "ggnewscale", "ggpmisc", "ggpointdensity", "ggpol", "ggpubr", "ggQC",
                                               "ggQQunif", "ggquiver", "ggradar", "ggraph", "ggrastr", "ggrepel",
                                               "ggridges", "ggsci", "ggseas", "ggshadow", "ggsom", "ggspectra",
                                               "ggstance", "ggstatsplot", "ggtech", "ggtext", "ggthemes", "ggTimeSeries",
                                               "ggwordcloud", "hrbrthemes", "lemons", "lindia", "plotROC", "qqplotr", "sugrrants",
                                               "survminer", "treemapify", "tvthemes", "xmrr",
                                               
                                               "hms", "lubridate", "feather", "haven", "httr", "jsonlite", "readxl",
                                               "rvest", "xml2", "modelr", "broom",
                                               
                                               "rsample", "recipes", "tune", "parsnip", "yardstick", "patchwork",
                                               "tidygraph")

# capture other internal core base packages
basePacks <- as.data.frame(installed.packages())
basePacks <- basePacks$Package[basePacks$Priority == "base" & !is.na(basePacks$Priority)]

functionTablePack$basePacks <- functionTablePack$package %in% basePacks

#                 Start of 2014 onwards ####

functionTablePack <- functionTablePack[functionTablePack$yearCreated >= 2014,]

# SUMMARY STATISTICS ####

summaryDf <- data.frame(repoCount = length(unique(functionTablePack$id.repo)),
                        scriptCount = length(unique(functionTablePack$ursID)),
                        userCount = length(unique(functionTablePack$userId)),
                        functionCount = sum(functionTablePack$count),
                        uniqueFunCount = length(unique(functionTablePack$function.)),
                        funOnce = sum(table(functionTablePack$function.)==1))

# SUBSET DATA ####
#                   actively maintained repos ####

# exclude repositories that have been updated > 12 months away from their original
# creation date. They are likely updated over time, potentially gitHub repos for
# R packages themselves.
funTableRepo <- functionTablePack[!duplicated(functionTablePack$urID), c("urID", 
                                                               "monthsSinceJan10",
                                                               "monthsSinceJan10Created")]
funTableRepo$updateLag <- funTableRepo$monthsSinceJan10 - funTableRepo$monthsSinceJan10Created

activeRepos <- funTableRepo$urID[funTableRepo$updateLag < 12]

sum(funTableRepo$updateLag >= 12) / nrow(funTableRepo)

functionTablePack <- functionTablePack[functionTablePack$urID %in% activeRepos,]

#                   functions with no identifiable package ####

noPack <- is.na(functionTablePack$package)

noPackFun <- is.na(functionTablePack$package[!duplicated(functionTablePack$function.)])
names(noPackFun) = functionTablePack$function.[!duplicated(functionTablePack$function.)]

sum(noPackFun)
sum(noPack) / nrow(functionTablePack)

noPackTable <- table(functionTablePack$function.[noPack])
sum(table(noPackTable)[1:5]) / sum(noPackFun)


functionTablePack <- droplevels(functionTablePack[!noPack,])

#                   rare functions ####
 
sampleThreshold = 0.001 # remove functions under a certain usage threshold

funList <- with(functionTablePack, table(yearCreated, function.))
yearScriptCount <- as.vector(table(functionTablePack$yearCreated[!duplicated(functionTablePack$ursID)]))

funProp <- funList / yearScriptCount
funBin = ifelse(funProp >= sampleThreshold, 1, 0)
funRep <- colSums(funBin)
funKeep <- names(funRep)[funRep > 0]
funLoss <- names(funRep)[funRep == 0]

functionTablePack$tidyverse <- functionTablePack$tidyverse | functionTablePack$tidyExtend

functionTablePack$other <- !functionTablePack$basePacks &
                           !functionTablePack$tidyverse

catLoss <- t(sapply(c("basePacks", "tidyverse", "other"), function(cat){

  return(cbind(total = sum(!duplicated(functionTablePack$function.[functionTablePack[,cat]])),
                    lose = sum(funLoss %in% functionTablePack$function.[functionTablePack[,cat]]),
                    keep = sum(funKeep %in% functionTablePack$function.[functionTablePack[,cat]])))
    
}))

catLoss[,3] / catLoss[,1]

funBaseLoss <- sum(funLoss %in% functionTablePack$function.[functionTablePack$basePacks])
funBaseKeep <- sum(funKeep %in% functionTablePack$function.[functionTablePack$basePacks])
funBaseLoss / funBaseKeep

funTableSub <- droplevels(functionTablePack[functionTablePack$function. %in% funKeep,])

sum(!duplicated(functionTablePack$function.))
sum(!duplicated(funTableSub$function.))

# remove scripts with only a single function ####
soloFuns <- tapply(funTableSub$count,
                   funTableSub$ursID, sum)
soloFuns <- soloFuns[soloFuns==1]

funTableSolo <- droplevels(funTableSub[funTableSub$ursID %in% names(soloFuns),])
funTableSub <- droplevels(funTableSub[!funTableSub$ursID %in% names(soloFuns),])

write.csv(funTableSub, "./outputs/commonFunctionLong.csv")

# FINAL SUMMARY STATISTICS ####

summaryDf <- rbind(summaryDf,
      data.frame(repoCount = length(unique(funTableSub$id.repo)),
                 scriptCount = length(unique(funTableSub$ursID)),
                 userCount = length(unique(funTableSub$userId)),
                 functionCount = sum(funTableSub$count),
                 uniqueFunCount = length(unique(funTableSub$function.)),
                 funOnce = NA))

# function category richness

sum(funTableSub$basePacks[!duplicated(funTableSub$function.)])
sum(funTableSub$tidyverse[!duplicated(funTableSub$function.)])
sum(funTableSub$tidyExtend[!duplicated(funTableSub$function.)])



summaryDf[3,] <- summaryDf[2,] / summaryDf[1,]

write.csv(summaryDf, "./outputs/summaryStats.csv")

# test sampling filters ####

testThresholds <- c(seq(0.0001, 0.001, 0.0001),
                    seq(0.001, 0.01, 0.001),
                    seq(0.01, 0.1, 0.01))

samplingList <- lapply(testThresholds, function(x){
print(x)
sampleThreshold = x # remove functions under a certain usage threshold

funList <- with(functionTablePack, table(yearCreated, function.))
yearScriptCount <- as.vector(table(functionTablePack$yearCreated[!duplicated(functionTablePack$ursID)]))

funProp <- funList / yearScriptCount
funBin = ifelse(funProp >= sampleThreshold, 1, 0)
funRep <- colSums(funBin)
funKeep <- names(funRep)[funRep > 0]
funLoss <- names(funRep)[funRep == 0]

functionTablePack$tidyverse <- functionTablePack$tidyverse | functionTablePack$tidyExtend

functionTablePack$other <- !functionTablePack$basePacks &
  !functionTablePack$tidyverse

catLoss <- t(sapply(c("basePacks", "tidyverse", "other"), function(cat){
  
  return(cbind(total = sum(!duplicated(functionTablePack$function.[functionTablePack[,cat]])),
               lose = sum(funLoss %in% functionTablePack$function.[functionTablePack[,cat]]),
               keep = sum(funKeep %in% functionTablePack$function.[functionTablePack[,cat]])))
  
}))

return(catLoss)
})

pdf("./plots/samplingFilterTest.pdf", height=4, width=5, useDingbats = FALSE)
par(mar=c(3,3,0.5,0.5), ps=10, mgp=c(3,0.5,0), las=1, tcl=-0.25)

plot(x=NULL, y=NULL, xlim=log(range(testThresholds))+c(0,0.1), ylim=log(c(1,4000)),
     axes=FALSE, xlab="", ylab="")

axis(side=1, at=log(testThresholds), labels=NA, tcl=-0.125)
axis(side=1, at=log(c(1e-4, 1e-3, 0.01, 0.1)), 
     labels=paste0(c(1e-4, 1e-3, 0.01, 0.1)*100,"%"), mgp=c(3,0.2,0))
mtext(side=1, line=2,
      text="Sampling filter\n(occurs in n of repositories in at least one month)")

axis(side=2, at=log(c(seq(1,10,1),
                      seq(10,100,10),
                      seq(100,1000,100),
                      seq(1000,4000,1000))),
     labels=NA, tcl=-0.125)
axis(side=2, at=log(c(0,1, 10, 100, 1000, 4000)+1),
     labels=c(0,1, 10, 100, 1000, 4000), las=1)
mtext(side=2, line=2,
      text="Function count that exceeded filter", las=0)

abline(v=log(0.001), lwd=2, lty="31")

tempKeep <- sapply(samplingList, function(x){
  sum(x[,3])
})
lines(log(tempKeep+1) ~ log(testThresholds), lwd=2)

text(x=rev(log(testThresholds))[10],
     y=rev(log(tempKeep+1))[10],
     pos=4,
     labels=c("All"), font=2)

sapply(1:3, function(n){
  
  tempKeep <- sapply(samplingList, function(x){
    x[n,3]
  })
  
  lines(log(tempKeep+1) ~ log(testThresholds),
        col = c("#2369bd", "#d12239", "grey75")[n], lwd=2)
  
  text(x=rev(log(testThresholds))[1],
       y=rev(log(tempKeep+1))[1],
       pos=4, offset=0.25,
       labels=c("B", "T", "O")[n], font=2,
       col=c("#2369bd", "#d12239", "grey75")[n])
  
  })
box()
dev.off()
