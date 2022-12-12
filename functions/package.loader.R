package.loader <- function(packages){

  # are there any packages that aren't already installed?
  new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  
  # If there's at least 1, install them
  if(length(new.packages) > 0){
    install.packages(new.packages, dependencies=T)
  }
  
  # then load the packages
  sapply(packages, require, character.only=T)
  
}