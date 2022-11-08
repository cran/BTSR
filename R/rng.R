##-------------------------------------------------------------------------
## internal function
## Assigns a initial value to the seed argument passed to fortran,
## based on the type of random number generation chosen by the user
##-------------------------------------------------------------------------
.seed.start <- function(seed, rngtype){
  ##------------------------------------------
  ##  Random number generators
  ##------------------------------------------
  ## 0: Original rng_uniform  (Jason Blevins)
  ## 1: Wichmann-Hill
  ## 2: Mersenne Twister
  ## 3: Marsaglia-MultiCarry (kiss 32)
  ## 4: Marsaglia-MultiCarry (kiss 64)
  ## 5: Knuth (2002)
  ## 6: L'Ecuyer's 1999 (64-bits)
  ##------------------------------------------

  if(rngtype == 0) add <- c(521288629, 362436069, 16163801, 1131199299)
  if(rngtype == 1) add <- c(3026, 3030, 3032)
  if(rngtype  %in% c(2,6)) add <- 521288629
  if(rngtype %in% c(3,4)) add <- c(123456789, 362436069, 521288629, 916191069)
  if(rngtype == 5) add <- c(153587801,-759022222,-759022222,-1718083407,-123456789)
  if(! rngtype %in% 0:6) stop("rngtype must be a number between 0 and 6")

  if(length(seed) != length(add)){
    if(length(seed) > length(add)){
      warning(paste("only the first ", length(add),
                    " values in seed will be used", sep = ""), immediate. = TRUE)
      seed <- seed[1:length(add)]
    }else{
      if(length(seed) > 1){
        warning("only the first value in seed will be used", immediate. = TRUE)
        seed <- seed[1]
      }
    }
  }	
	
  seed <- seed + add
	if(rngtype == 1){
		# These values should be positive integers between 1 and 30,000.
		seed <- sapply(seed, function(x) min(x, 30000))
	}
	
  return(as.integer(seed))
}


