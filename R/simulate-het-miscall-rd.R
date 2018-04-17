

#' Simulate genotypes under read-depth-dependent het miscall
#' 
#' This takes an exising 012 file to match allele freqs and pattern of missing data from it.
#' @param Y an N x L 012,-1 matrix
#' @param R an N x L matrix of read depth categories.  Must be labeled 1, 2, ...
#' @param M a vector of miscall rates.  Will be recycled to the number of categories in R.
simulate_het_miscall_rd <- function(Y, R, M) {
  
  # get allele freqs
  G <- Y
  G[G==-1] <- NA
  ref_f <- 1.0 - (colSums(G, na.rm = TRUE) / (2 * colSums(!is.na(G))))
  
  # now simulate genotypes by independently segregating 0s and 1s
  pp <- rep(ref_f, each = nrow(Y))
  Y_sim <- matrix(
    as.integer(pp < runif(n = length(pp))) + 
    as.integer(pp < runif(n = length(pp))), 
    nrow = nrow(Y))
  
  # no-call all that were no-called in the original data
  # and those that have read depth category less than 1
  Y_sim[Y == -1] <- -1
  Y_sim[R <= 0] <- -1
  
  # now we want to simulate het miscalls
  num_cats <- length(unique(as.vector(R[R>0])))
  mvec <- rep(M, length.out = num_cats)  # recycle M to a vector with entries for each read depth category
  

  hets <- Y_sim == 1  # positions of hets
  Rhets <- R[hets]  # read depth categories of the hets
  mhets <- mvec[Rhets]  # miscall rates of the hets
  
  # now we simulate a 1, 0, or -1 that we add to each het's 1 value to do the miscall
  ones <- as.integer(runif(n = length(mhets)) < mhets)
  sgn <- sample(c(-1, 1), size = length(mhets), replace = TRUE)
  
  newhets <- rep(1, length(mhets)) + ones * sgn
  
  Y_mis <- Y_sim
  Y_mis[hets] <- newhets 
  
  list(Y_no_err = Y_sim, 
       Y_with_err = Y_mis,
       R = R)
}