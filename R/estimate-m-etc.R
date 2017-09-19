

#' simulate values for the genotypes given the observed genotype, the est alle freq, and the gtyp error rate
#' 
#' This is a helper function for the estimate_m function
#' @param D an 012,-1 matrix of observed genotypes
#' @param p the estimated allele freqs
#' @param m the genotyping error rate (must be a scalar)
simulate_genos_from_posterior <- function(D, p, m) {
  stopifnot(length(m) == 1)
  
  glist <- lapply(1:ncol(D), function(i) {
    obs <- D[, i] # the observed genotypes
    pl <- p[i]  # the alle freq at locus i
    post0 <- c(
      (1 - pl) / (1 - pl + m * pl),  # posterior that observed 0 is truly a 0
      (m * pl) / (1 - pl + m * pl)   # posterior that observed 0 is truly a 1
    )
    post2 <- c(
      (m * (1 - pl)) / (pl + m * (1 - pl)),  # posterior that observed 2 is truly a 1
      pl / (pl + m * (1 - pl))               # posterior that observed 2 is truly a 2
    )
    obs[obs == 0] <- sample(x = c(0, 1), size = sum(obs == 0), replace = TRUE, prob = post0)
    obs[obs == 2] <- sample(x = c(1, 2), size = sum(obs == 2), replace = TRUE, prob = post2)
    obs
  })
  
  # then turn it into a matrix with the same dimensions and dimnames as D
  ret <- matrix(unlist(glist), nrow = nrow(D))
  dimnames(ret) <- dimnames(D)
  ret
}


#' estimate the heterozygote miscall rate from a simple model
#' 
#' @param dat012 An 012 matrix.  Missing data can be -1 or NA
#' @param nreps number of MCMC sweeps to do
#' @param m_init initial starting value for m must be between 0 and 1
#' @param a0 beta parameter for reference alleles
#' @param a1 beta parameter for alternate alleles
#' @param sm standard devation of proposal distribution for m
estimate_m <- function(dat012,
                       nreps = 200,
                       m_init = runif(1),
                       a0 = 0.5,
                       a1 = 0.5,
                       sm = 0.005
) {
  
  stopifnot(m_init > 0 & m_init < 1)
  
  D <- dat012
  D[is.na(D)] <- -1
  
  # get the N variables
  N0 <- colSums(D == 0)
  N1 <- colSums(D == 1)
  N2 <- colSums(D == 2)
  
  # initialize the Zs to the Ns
  Z0 <- N0
  Z1 <- N1
  Z2 <- N2
  
  # make some place to return the m values visited
  m <- rep(NA, nreps)
  m[1] <- m_init
  
  # then do the sweeps
  for (r in 2:nreps) {
    
    # new estimate of frequency of the "1" allele from Gibbs sampling
    p <- rbeta(n = length(Z0), 
               shape1 = a1 + 2 * Z2 + Z1, 
               shape2 = a0 + 2 * Z0 + Z1)
    
    # propose then accept or reject a new value for m
    mprop <- m[r - 1] + rnorm(1, 0, sm)
    reject <- TRUE  # reject it unless we don't
    if (mprop > 0 & mprop < 1) {
      numer <- sum(N0 * log((1 - p)^2 + mprop * p * (1 - p)) +
                     N1 * log((1 - mprop) * 2 * p * (1 - p)) +
                     N2 * log(p ^ 2 + mprop * p * (1 - p)))
      denom <- sum(N0 * log((1 - p)^2 + m[r - 1] * p * (1 - p)) +
                     N1 * log((1 - m[r - 1]) * 2 * p * (1 - p)) +
                     N2 * log(p ^ 2 + m[r - 1] * p * (1 - p)))
      if (log(runif(1)) < numer - denom) {
        reject <- FALSE
      }
    }
    if (reject == FALSE) {
      m[r] <- mprop
    } else {
      m[r] <- m[r - 1]
    }
    
    # new values for Z from Gibbs sampling
    A0 <- rbinom(n = length(N0), size = N0, prob = (m[r] * p) / (1 - p + m[r] * p))
    A2 <- rbinom(n = length(N2), size = N2, prob = (m[r] * (1 - p)) / (p + m[r] * (1 - p)))
    
    Z0 <- N0 - A0
    Z1 <- N1 + A0 + A2
    Z2 <- N2 - A2
    
  }
  # return m, and eventually I need to also return the final Zs and the Ns
  # and I may as well return a new 012 file with "corrected" genotypes, which 
  # I can make by broadcasting the Zs around, for example...
  
  # inferring/realizing/simulating genotypes. I can simulate these from their posterior
  # given the estimated allele freq and the observed genotype.  To do this I will cycle
  # over the columns (the snps) in D, and for each one, I will compute the posterior of the
  # the genotype given the observed genotype (only have to for 0's and 2's) and then I will
  # sample from those posteriors.  We have a separate function that does this
  ret <- list()
  ret$simmed_genos <- simulate_genos_from_posterior(D, p, m[nreps])
  
  # compute an overall genotyping error rate
  diff <- ret$simmed_genos != D
  diff[D == -1] <- NA
  ret$overall_geno_err_est <- mean(diff, na.rm = TRUE)
  
  # return the trace of m values
  ret$m <- m
  
  ret
}