hillCalc <- function(comm, q = 1, l=NULL){
  
  s = length(comm) # we need number of species
  p = comm / sum(comm) # their relative abundance
  r = 1/p # and the inverse, their relative rareness
  
  if(is.null(l)){l = 1-q} # if we don't have l, calculate it from q
  
  #>>>>>>>>>>>>>>>>>>>
  # I'm gonna break the equation up into steps so you can see each bit in action.
  #>>>>>>>>>>>>>>>>>>>
  
  # Here we take each species one at a time, multiply its relative abundance
  # by its relative rareness to the power of our "l" exponent.
  # This is where we choose how much we want to weight towards common species.
  # if l = 0, we just get our relative abundances back
  theInsideSumBit <- sapply(1:length(p),
                            function(n){
                              p[n] * (r[n]^l)
                            })
  
  # then we sum these together
  theSumBit <- sum(theInsideSumBit)
  
  # then finally we raise the sum to 1/l, which reverses the change in units from the l above,
  # and converts the scores back to "diversity" units.
  theExpBit <- theSumBit^(1/l)
  
  # and that's our hill number!
  return(theExpBit)
}