

#' Mocking up an R-version of estimating m with read depth but no nulls
#' 
#' I had originally thought that I would do this in RCpp, but I think 
#' I will mock it up here first in R and then write in RCpp if it works
#' well but needs to be faster.
#' @param Y the 012,-1 matrix that is N x L giving the observed genotypes of the N individuals
#' at L SNPs.
#' @param R integer matrix that is N x L giving the read depth categories.
#' @param init_m initial value to use for all depth categories of m
estimate_m_with_rd_and_nulls <- function(Y, R, init_m) {
  
  # get data to start developing:
  Y <- read_rds(path = "dev_stuff/g012_mat.rds")
  R <- read_rds(path = "dev_stuff/dp_mat.rds")
  
  # start off with some initialization.  
  # set starting allele frequencies for the reference and alt alleles at each locus
  tmp <- Y
  tmp[tmp==-1] <- NA
  n1 <- colSums(tmp == 1, na.rm = TRUE) + 2 * colSums(tmp == 2, na.rm = TRUE)
  n0 <- colSums(tmp == 1, na.rm = TRUE) + 2 * colSums(tmp == 0, na.rm = TRUE)
  
  # make the null freq about 0.08 to start
  p <- rbind(n0, n1)
  p <- p / rep(colSums(p), each = 2)
  
  # initialize m to be 0.1 across the boards
  m <- rep(0.1, length(table(R)))
  
  # now compute the full conditionals for the X's (considering only 0, 1, and 2)
  # here is a matrix of priors on the genotypes:
  gp <- rbind(p[1,] ^ 2, 2 * p[1,] * p[2,], p[2,] ^ 2)
  
}