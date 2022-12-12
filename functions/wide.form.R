wide.form <- function(group.by,
                      spread.by,
                      values,
                      fill = 0 ,
                      prop = FALSE,
                      prune = TRUE){
  
  sum.mat <- tapply(values,
                    list(group.by, spread.by),
                    sum)
  
  sum.mat[is.na(sum.mat)] = fill
  
  if(prune){
    
    sum.mat <- sum.mat[rowSums(sum.mat) > 0, colSums(sum.mat) > 0]
    
  }
  
  if(prop){
    
    sum.mat = prop.table(sum.mat, 1)
    
    }
  
  return(sum.mat)
}
