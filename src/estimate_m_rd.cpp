#include <Rcpp.h>
using namespace Rcpp;



//' just a quick function for making an 012 matrix
//' 
//' The standard way within R of pulling values out of a named
//' vector really bogs down.  So I will do this instead.
// [[Rcpp::export]]
IntegerMatrix make_it_012(CharacterMatrix M) {
  int N = M.ncol();
  int L = M.nrow();
  int i,l;
  
  IntegerMatrix ret(L, N);
  
  for(l=0;l<L;l++) {
    for(i=0;i<N;i++) {
      if(strcmp(M(l,i), "0/0") == 0) {
        ret(l,i) = 0;
      } else if(strcmp(M(l,i), "0|0") == 0) {
        ret(l,i) = 0;
      } else if(strcmp(M(l,i), "0/1") == 0) {
        ret(l,i) = 1;
      } else if(strcmp(M(l,i), "1/0") == 0) {
        ret(l,i) = 1;
      } else if(strcmp(M(l,i), "0|1") == 0) {
        ret(l,i) = 1;
      } else if(strcmp(M(l,i), "1|0") == 0) {
        ret(l,i) = 1;
      } else if(strcmp(M(l,i), "1/1") == 0) {
        ret(l,i) = 2;
      } else if(strcmp(M(l,i), "1|1") == 0) {
        ret(l,i) = 2;
      } else {
        ret(l,i) = -1;  // this takes care of NAs
      }
    }
  }
  
  return(ret);
}



//' compute full conditional for each X given Y, p, R, and m, and then sample from it
//'
//' This just writes new values into X as if it were an output variable 
// [[Rcpp::export]]
void gibbsX(IntegerMatrix X, IntegerMatrix Y, IntegerMatrix R, NumericVector p, NumericVector M) {
  int i,l,k;
  int N = Y.nrow();
  int L = Y.ncol();
  double probs[3];  // for storing the probs
  double pp, m, normo, rr, cumul;
  int y, r, newX;
  
  // cycle over each individual and locus
  for(i=0;i<N;i++) {
    for(l=0;l<L;l++) {
      pp = p(l);  // get the reference alle freq
      // intialize the probs to the prior from the allele freq
      probs[0] = pp * pp;
      probs[1] = 2.0 * pp * (1 - pp);
      probs[2] = (1 - pp) * (1 - pp);
      
      y = Y(i,l); // temp variable to hold genotype
      r = R(i,l) - 1;   // temp variable to hold read category.  Subtract 1 for base-0 indexing...
      if(y >= 0 && r >= 0) {  // if we have a genotype call and a read depth category then use those to update the priors 
        m = (double)(M(r));
        probs[0] *= (y == 0);
        probs[1] *= ( (m / 2.0) * (y == 0) + (1.0 - m) * (y == 1) + (m / 2.0) * (y == 2));
        probs[2] *= (y == 2);
      }
      
      // now normalize those
      normo = 0.0;
      for(k=0;k<3;k++) normo += probs[k];
      for(k=0;k<3;k++) probs[k] /= normo;
      
      // now sample from them
      rr = R::runif(0,1);;
      cumul = 0.0;
      newX = -2;
      for(k=0;k<3;k++) {
        cumul += probs[k];
        if(cumul > rr) {
          newX = k;
          break;
        }
      }
      
      // and assign that to X(i,l)
      X(i,l) = newX;
        
    }
  }
  
} 



//' simulate new reference allele frequencies from their beta full conditional
//'
//' This just writes new values into P as if it were an output variable 
// [[Rcpp::export]]
void gibbsP(NumericVector p, IntegerMatrix X, NumericVector pri) {
  int i,l;
  int N = X.nrow();
  int L = X.ncol();
  
  double x0;   // to count up the number of zero alleles
  double x1;   // to count up the number of one alleles
  
  for(l=0;l<L;l++) {
    x0 = pri(0);  // initialize to the prior values
    x1 = pri(1);
    for(i=0;i<N;i++) {
      x0 += 2.0 * (X(i,l) == 0) + 1.0 * (X(i,l) == 1);
      x1 += 2.0 * (X(i,l) == 2) + 1.0 * (X(i,l) == 1);
    }
    p(l) = R::rbeta(x0, x1);
  }
}



//' simulate a new miscall rate for each read depth category given X and Y
//'
//' This just writes new values into M as if it were an output variable 
// [[Rcpp::export]]
void gibbsM(NumericVector M, int num_cats, IntegerMatrix X, IntegerMatrix Y, IntegerMatrix R, NumericVector pri) {
  int i,l,k;
  int N = Y.nrow();
  int L = Y.ncol();
  int r, x, y;
  
  NumericVector mc(num_cats, pri(0));  // for counting up the number of hets that were miscalled
  NumericVector cc(num_cats, pri(1));  // for counting up the number of hets that were correctly called
  
  // Note, the above two vectors are initialized to the beta prior for the miscall rate
  
  // cycle over all loci and all individuals
  for(l=0;l<L;l++) {
    for(i=0;i<N;i++) {
      r = R(i, l) - 1;  // the read depth category
      y = Y(i, l);
      x = X(i, l);
      if(r >= 0 && y >= 0) { // as long as genotype is not missing and read depth category is valid
        if(x == 1 && y != 1) mc(r) += 1.0;
        if(x == 1 && y == 1) cc(r) += 1.0;
      }
    }
  }
  
  // now simulate new values from a beta
  for(k=0;k<num_cats;k++) {
    M(k) = R::rbeta(mc(k), cc(k));
  }
  
}
  

//' Estimate heterozygote miscall rate for different read depth categories (no nulls)
//' 
//' This is pretty simple.
//' @param Y the 012,-1 matrix that is N x L giving the observed genotypes of the N individuals
//' at L SNPs.
//' @param R integer matrix that is N x L giving the read depth categories.  These must be indexed from
//' 1 up to num_cats.  Missing data should be -1.
//' @param init_m starting value for m.  Typically you might want to use the m estimated
//' from init_m
//' @param num_cats the number of read depth categories.  T
// [[Rcpp::export]]
List estimate_m_rd(IntegerMatrix Y, 
              IntegerMatrix R,
              double init_m,
              int num_cats,
              NumericVector p_prior,
              NumericVector m_prior,
              int num_reps) {
  
  int i,l, rep;
  int N = Y.nrow();
  int L = Y.ncol();

  // This is to hold the latent genotypes
  IntegerMatrix X(N, L);
  IntegerVector x0(L);   // to count up the number of zero alleles
  IntegerVector x1(L);   // to count up the number of one alleles
  NumericVector p(L); // for reference allele freq
  NumericMatrix Mtrace(num_cats, num_reps + 1);
  
  // this is to hold the het miscall rates
  NumericVector M(num_cats, init_m);
  
  // first compute the frequency of the reference allele
  for(l=0;l<L;l++) {
    x0(l) = 0; // this is like having a beta(1/2, 1/2) prior
    x1(l) = 0;
  }
  
  
  for(l=0;l<L;l++) {
    for(i=0;i<N;i++) {
      x0(l) += 2 * (Y(i,l) == 0);
      x0(l) += 1 * (Y(i,l) == 1);
      x1(l) += 1 * (Y(i,l) == 1);
      x1(l) += 2 * (Y(i,l) == 2);
    }
    p(l) = (double)(x0(l) + 0.5) / (double)(x0(l) + x1(l) + 1.0);
  }

  // then do the reps, and store the M vector each time
  Mtrace(_, 0) = M;
  for(rep=0; rep<num_reps; rep++) {
    gibbsX(X, Y, R, p, M);
    gibbsP(p, X, p_prior);
    gibbsM(M, num_cats, X, Y, R, m_prior);
    Mtrace(_, rep + 1) = M;
  }
  

  return List::create(_["M"] = M, 
                      _["N"] = N, 
                      _["L"] = L,
                      _["p"] = p,
                      _["x0"] = x0,
                      _["x1"] = x1, 
                      _["Y"] = Y, 
                      _["X"] = X,
                      _["Mtrace"] = Mtrace);
}

