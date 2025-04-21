//// make particles for potts examples ////

//// ABC particles ////
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// Generate summary statistics of auxiliary variables for given particles
mat pAuxSamp(Rcpp::RawVector foo, int ncolor, int cycle, Rcpp::NumericVector Designmat, int m, int num){
  
  int thnrow = Designmat.size();             // number of design points   
  //int thnrow = Designmat.n_rows;  
  mat H(thnrow,m);                       // cube summary statistics will be stored. (thnrow by m)
  omp_set_num_threads(num);
  
  
  //////    START OF BIGGEST CHAIN (M)     ////// 
  for (int i = 0; i < thnrow; i++) {
#pragma omp parallel for schedule(static)
    for (int M = 0; M < m; M++) {
      double thetaprop = Designmat[i];
      double sumstat = potts_stat(foo, ncolor, thetaprop, cycle);
      H(i, M) = sumstat;
    }
  }
  return(H);        	
}

