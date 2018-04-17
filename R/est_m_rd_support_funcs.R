


#' bin read depths of SNPs into categories having at least S observations
#' 
#' @param D a matrix of read depths.  Rows are individuals, columns are SNPs.  Cells where are missing
#' in the genotype matrix must be denoted as NA
#' @param S the min number of observations to have in each bin
#' @return This returns a list with two components.  \code{dp_bins} is a matrix of the same
#' shape as D with the bin categories (as 1, 2, ...) and -1 for this cells 
#' corresponding to missing genotypes.  \code{num_cats} is the number of depth bins.
#' \code{tidy_bins} is a long format description of the bins.
#' \code{bin_stats} is a tibble giving summary information about the read depth bins which
#' is useful for plotting things, etc.
bin_depths <- function(D, S) {
  
  # first, count things up and get them into a reasonable format
  dcnts <- as_tibble(as.data.frame(table(D))) %>%
    mutate(D = as.integer(as.character(D))) %>%
    setNames(c("D", "Freq")) %>%
    rename(dp = D,
           n = Freq) %>%
    arrange(dp) # just to make sure it is in increasing order
  
  # now we lump them up.  We just use a for loop for this
  idx <- 1
  sum <- 0
  n <- dcnts$n
  res <- rep(NA, length(n))
  for(i in seq_along(n)) {
    res[i] <- idx
    sum <- sum + n[i]
    if(sum > S) {
      idx <- idx + 1
      sum <- 0
    }
  }
  
  # at the end of that, we take the very last bin and just merge is with the 
  # preceding one, so that it doesn't have just a small number in it.
  res[res == max(res)] <- max(res) - 1
  
  
  # and we add it on there with a mutate which will bark an error is something has gone awry
  tidy_bins <- dcnts %>%
    mutate(bin = as.integer(res))
  
  # now we have to cut the original matrix up into these categories
  left_endpoints <- tidy_bins %>%
    group_by(bin) %>%
    summarise(ends = max(dp) + 0.2)
  
  cut_vec <- c(0, left_endpoints$ends)
  
  
  dp_bins <- D
  dp_bins[] <- cut(D, breaks = cut_vec)
  
  dp_bins[is.na(dp_bins)] <- -1L
  
  # finally summarize the bin stats
  bin_stats <- tidy_bins %>%
    group_by(bin) %>%
    summarise(total_n = sum(n),
              mean_dp = sum(dp * n) / sum(n))
  
  list(dp_bins = dp_bins,
       num_cats = max(tidy_bins$bin),
       tidy_bins = tidy_bins,
       bin_stats = bin_stats
       )
  
}


#' convert a VCF into an 012,-1 matrix and read_depth bin matrix for estimation
#' 
#' @param v a vcfR object into which a VCF file has been read
#' @param DF Field to use for obtaining total read depth.  Choices are DP and AD, but only DP is
#' implemented at this point. And you have to make sure that DP is a field that exists...
#' @param minBin minimum number of observations for each read depth bin
#' @return This sends back a list with \code{mat012}: the 012,-1 matrix of genotypes.
#' \code{dp_bins_list}: the list returned by bin_depths().
prep_vcf_for_est_m_rd <- function(v, DF, minBin) {
  
  # check to make sure the field named in DF exists.
  #vt <- vcfR2tidy(v, info_only = TRUE)
  #if(!(DF %in% vt$meta$ID)) {
  #  stop("The tag ", DF, " does not appear to be in the data...")
  #}
  
  # extract the matrix as an 012 file
  dgt <- extract.gt(v, element = "GT")
  
  d012 <- make_it_012(dgt)
  dimnames(d012) <- dimnames(dgt)
  
  #### Used to do it this way, but it was super slow, so I wrote a CPP function (make_it_012) to do it
  #' reps <- c(
  #'   "0/0" = 0,
  #'   "0|0" = 0,
  #'   "0/1" = 1,
  #'   "0|1" = 1,
  #'   "1/0" = 1,
  #'   "1|0" = 1,
  #'   "1/1" = 2,
  #'   "1|1" = 2
  #' )
  #' 
  #' dgt[] <- reps[dgt]
  #' storage.mode(dgt) <- "integer"
  #' dgt[is.na(dgt)] <- -1L
  #' 
  #' # so, that is our 012 matrix (but will need to be transposed)
  #' d012 <- dgt 
 
  
  # now get the read depth matrix
  dp <- extract.gt(v, element = "DP")
  storage.mode(dp) <- "integer"
  dp[d012 == -1] <- NA
  
  # OK, now dp is a matrix of read depths with NAs where the genotype was unobserved
  
  # now we need to bin those depths up
  bd <- bin_depths(D = dp, S = minBin)
  bd$dp_bins <- t(bd$dp_bins)
  
  
  list(mat012 = t(d012),
       dp_bins_list = bd)
}



#' tidy up the estimate_m_rd output into something you can plot
#' @param E the list returned by the estimation function
#' @param S the bin stats
tidy_m_ests <- function(E, S, burn = 50) {
  
  tibble(bin = 1:nrow(E$Mtrace),
    mean = rowMeans(E$Mtrace[, -(1:burn)]),
    lo95 = apply(E$Mtrace[, -(1:burn)], 1, quantile, probs = 0.05),
    hi95 = apply(E$Mtrace[, -(1:burn)], 1, quantile, probs = 0.95)
  ) %>%
    left_join(S)
}
