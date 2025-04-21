#include <RcppArmadillo.h>
#include <limits>
#include <omp.h>
#define min(x,y) (((x) < (y)) ? (x) : (y))

using namespace std;
using namespace Rcpp;
using namespace arma;

////// ERGM example //////
double Choose(int x, int y){
  double result1 = 1, result2 = 1, result;
  int iter = y;
  
  if( x< y ){ result = 0;	}else{	
    for(int i = 0; i<iter; i++){
      result1 = result1*x;
      result2 = result2*y;
      y = y-1;
      x = x-1;   	
    }	
    
    result = result1/result2;
  }
  return(result);
}

// [[Rcpp::depends("RcppArmadillo")]]
// Count degree of each node 
vec countDegree(vec rowsum){
  int nrow = rowsum.n_elem, ind;
  vec degree = zeros( nrow );
  for(int i = 0; i < nrow; i++ ){
    ind = rowsum(i); 
    if(ind>0){ degree(ind-1) = degree(ind-1) + 1; }
  }
  return(degree);
}

// [[Rcpp::depends("RcppArmadillo")]]
// Count edgewise shared partners 
vec countShared(vec rowsum, mat X){
  
  int numedges = sum(rowsum)/2, nrow = rowsum.n_elem ,ind;
  vec hist = zeros(numedges);
  int num = 0;
  for(int k = 0; k< nrow; k++){
    for(int j = 0; j< k+1; j++){
      if( X(k,j) == 1){  for(int i = 0; i<nrow; i++){ hist(num) =  hist(num) + X(i,k)*X(k,j)*X(j,i); }
      num = num + 1; }
    }
  }
  vec shared = zeros(nrow);
  for(int i = 0; i < numedges; i++ ){
    ind = hist(i);	
    if(ind>0){ shared(ind-1) = shared(ind-1) + 1; }
  }
  return(shared);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// Calculate summary statistics
vec Summary(mat X, vec grade, vec sex){
  int nrow = X.n_rows;
  double decay = 0.25;
  vec rowsum = sum(X,1), result = zeros(9);  // 크기를 9로 조정 (0~8)
  
  // count degree of each node
  vec degree = countDegree(rowsum);
  
  // count edgewise shared partners 
  vec shared = countShared(rowsum,X);
  
  // Calculate summary statistics
  for(int i = 0; i< nrow; i++){
    result(0) = result(0) + Choose(rowsum(i),1);
    result(7) = result(7) + ( 1 - pow( (1-exp(-decay)), i+1) ) * degree(i);
    result(8) = result(8) + ( 1 - pow( (1-exp(-decay)), i+1) ) * shared(i);
    
    for(int j = 0; j< i; j++){
      result(1) = result(1) + X(i,j)*( (grade[i]==7)*(grade[j]==7)   );
      result(2) = result(2) + X(i,j)*( (grade[i]==8)*(grade[j]==8)   );
      result(3) = result(3) + X(i,j)*( (grade[i]==9)*(grade[j]==9)   );	
      result(4) = result(4) + X(i,j)*( (grade[i]==10)*(grade[j]==10) );
      result(5) = result(5) + X(i,j)*( (grade[i]==11)*(grade[j]==11) );		
      result(6) = result(6) + X(i,j)*( (grade[i]==12)*(grade[j]==12) );		   
    }
  }
  result(0) = result(0) / 2;
  result(7) = exp(decay) * result(7);
  result(8) = exp(decay) * result(8);   
  
  return result;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// random scan gibbs update
vec ergm_gibbs(mat X, vec grade, vec sex, vec coef, int cycle){
  int nrow = X.n_rows, indi, indj, prei, prej, res;
  double decay = 0.25;
  vec star = sum(X,1), changestat = zeros(9), Sumstat = Summary(X,grade,sex);
  rowvec ivec = trans(zeros(nrow)), jvec = trans(zeros(nrow));
  vec geoweight = zeros(nrow);
  for(int i = 0; i < nrow; i++){ 
    geoweight(i) = 1 - pow( (1 - exp(-decay)), i + 1); 
  }
  
  for(int l = 0; l < cycle; l++){
    for(int i = 1; i < nrow; i++){
      for(int j = 0; j < i; j++){
        ivec(i) = 1; jvec(j) = 1;
        res = 0; 
        
        // When Xij is 0
        if(X(i,j) == 0){  
          indi = star(i) + 1, indj = star(j) + 1; 
          // change statistics of edge
          changestat(0) = ( Choose(indi,1) + Choose(indj,1) - Choose(star(i),1) - Choose(star(j),1) ) / 2;
          
          // change statistics of nodematch
          changestat(1) = ( (grade[i] == 7) * (grade[j] == 7) );
          changestat(2) = ( (grade[i] == 8) * (grade[j] == 8) );
          changestat(3) = ( (grade[i] == 9) * (grade[j] == 9) );	       
          changestat(4) = ( (grade[i] == 10) * (grade[j] == 10) );
          changestat(5) = ( (grade[i] == 11) * (grade[j] == 11) );            
          changestat(6) = ( (grade[i] == 12) * (grade[j] == 12) );       		               
          
          // change statistics of gwd
          changestat(7) = exp(decay) * ( geoweight(indi - 1) + geoweight(indj - 1) );
          if(indi - 1 > 0){ changestat(7) = changestat(7) - exp(decay) * geoweight(indi - 1 - 1); }
          if(indj - 1 > 0){ changestat(7) = changestat(7) - exp(decay) * geoweight(indj - 1 - 1); }
          
          // change statistics of gwes
          changestat(8) = 0;     
          for(int k = 0; k < nrow; k++){
            if( ( X(i,k) == 1 ) && ( X(j,k) == 1  ) ){
              prei = sum( X.row(i) % X.row(k) );
              prej = sum( X.row(j) % X.row(k) );
              
              changestat(8) = changestat(8) + exp(decay) * geoweight(prei + 1 - 1); 
              if(prei > 0){ changestat(8) = changestat(8) - exp(decay) * geoweight(prei - 1);} 
              changestat(8) = changestat(8) + exp(decay) * geoweight(prej + 1 - 1); 
              if(prej > 0){ changestat(8) = changestat(8) - exp(decay) * geoweight(prej - 1);} 
              res = res + 1; // X(k,i) multiply X(i,j) +1 multiply X(j,k)
            }
          }
          if(res > 0){ changestat(8) = changestat(8) + exp(decay) * geoweight(res - 1); }
          
          // accept reject step  	        
          // probability of Xij = 1 for given others are fixed 
          vec r = exp(trans(coef) * changestat);         
          double p = r(0) / (1 + r(0));   		   
          if( randu() < p  ){
            X(i,j) = X(j,i) = 1; 
            star(i) = indi; star(j) = indj;		 
            Sumstat = Sumstat + changestat;
          }
          
          // When Xij is 1  
        } else {
          indi = star(i) - 1, indj = star(j) - 1;	
          // change statistics of edge
          changestat(0) = ( Choose(star(i),1) + Choose(star(j),1) - Choose(indi,1) - Choose(indj,1) ) / 2;
          
          // change statistics of nodematch
          changestat(1) = ( (grade[i] == 7) * (grade[j] == 7) );
          changestat(2) = ( (grade[i] == 8) * (grade[j] == 8) );
          changestat(3) = ( (grade[i] == 9) * (grade[j] == 9) );	       
          changestat(4) = ( (grade[i] == 10) * (grade[j] == 10) );
          changestat(5) = ( (grade[i] == 11) * (grade[j] == 11) );            
          changestat(6) = ( (grade[i] == 12) * (grade[j] == 12) );       		               
          
          // change statistics of gwd
          changestat(7) = exp(decay) * ( geoweight(indi - 1 + 1) + geoweight(indj - 1 + 1) );
          if(indi - 1 >= 0){ changestat(7) = changestat(7) - exp(decay) * geoweight(indi - 1); }       
          if(indj - 1 >= 0){ changestat(7) = changestat(7) - exp(decay) * geoweight(indj - 1); }        		   
          
          // change statistics of gwesp 
          changestat(8) = 0;
          for(int k = 0; k < nrow; k++){
            if( ( X(i,k) == 1 ) && ( X(j,k) == 1  ) ){
              prei = sum( X.row(i) % X.row(k) );
              prej = sum( X.row(j) % X.row(k) );
              if(prei - 1 > 0){ changestat(8) = changestat(8) - exp(decay) * geoweight(prei - 1 - 1); } 
              changestat(8) = changestat(8) + exp(decay) * geoweight(prei - 1); 
              if(prej - 1 > 0){ changestat(8) = changestat(8) - exp(decay) * geoweight(prej - 1 - 1); } 
              changestat(8) = changestat(8) + exp(decay) * geoweight(prej - 1); 
              res = res + 1; // X(k,i) multiply X(i,j)  multiply X(j,k)
            }
          }
          if(res > 0){ changestat(8) = changestat(8) + exp(decay) * geoweight(res - 1); }
          
          // accept reject step  	        
          // probability of Xij = 1 for given others are fixed 
          vec r = exp(trans(coef) * changestat);         
          double p = r(0) / (1 + r(0));   		   
          if( randu() > p  ){
            X(i,j) = X(j,i) = 0; 
            star(i) = indi; star(j) = indj;  
            Sumstat = Sumstat - changestat;
          }             
        }
        // End of a single update  
        ivec(i) = 0; jvec(j) = 0;
      }
    }   
  }
  
  return(Sumstat); 
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// Generate summary statistics of auxiliary variables for given particles
List pAuxSamp_ergm(mat X, vec grade, vec sex, int cycle, int num, int numcore, mat Designmat){
  
  omp_set_num_threads(numcore);
  int numpoint = Designmat.n_rows;       // number of designpoints
  int ndim = Designmat.n_cols;           // parameter dimension
  
  mat Mu(numpoint, ndim);
  mat Cov(numpoint*ndim, ndim);
  //parallel computing
  int i;
#pragma omp parallel shared(Designmat) private(i)
{	
#pragma omp for schedule(static)  
  for(i = 0; i < numpoint; i++){
    mat pseudo(num,ndim);
    for(int j = 0; j < num; j++){ pseudo.row(j) =  trans( ergm_gibbs(X, grade, sex, trans( Designmat.row( i ) ), cycle) ); }
    Mu.row(i) =  mean(pseudo,0) ;    
    Cov.submat( span(ndim*i, ndim*i+ndim-1), span(0, ndim-1) ) = cov(pseudo);
  }
}

return List::create(Named("Cov") = Cov, Named("Mu") = Mu);	
}

// [[Rcpp::export]]
// Generate summary statistics of auxiliary variables for given particles to construct IS estimate
mat pResponse_ergm(mat X, int cycle, vec hatparameter, int m, const vec& Grade, const vec& Sex, int num) {
  int n_stats = hatparameter.n_elem;
  mat H(m, n_stats, fill::zeros); 
  
  omp_set_num_threads(num);  // parallel computing
  
#pragma omp parallel shared(H) 
{
#pragma omp for schedule(static)
  for (int M = 0; M < m; M++) {
    vec sumstat = ergm_gibbs(X, Grade, Sex, hatparameter, cycle); // Gibbs update
    for (int stat = 0; stat < n_stats; stat++) {
      H(M, stat) = sumstat[stat];
    }
  }
}

return H;
}

// [[Rcpp::export]]
mat computeSigma(const mat& Designmat, const vec& phi, double sigmasq) {
  int n = Designmat.n_rows; 
  int d = Designmat.n_cols; 
  
  mat Sigma(n, n, fill::zeros);
  
  for (int i = 0; i < n; i++) {
    for (int j = 0; j <= i; j++) {
      double distance = 0.0;
      for (int k = 0; k < d; k++) {
        double diff = fabs(Designmat(i, k) - Designmat(j, k));
        distance += pow(diff / phi[k], 2.0); // normalize with phi
      }
      distance = sqrt(distance); // Euclidean distance
      
      double matern_kernel = sigmasq * (1 + sqrt(3) * distance) * exp(-sqrt(3) * distance);
      Sigma(i, j) = Sigma(j, i) = matern_kernel;
    }
  }
  return Sigma;
}


// [[Rcpp::export]]
// GP emulator
List ergm_GPEm(int Niter, mat theta, mat COV, double lhXZ, vec betahat, vec phi, double sigmasq, mat Designmat, vec y) {
  int thnrow = Designmat.n_rows;
  int nCOVcols = COV.n_cols;
  vec thetaprev(nCOVcols, fill::zeros); 
  double lhXZp, logprob, u;
  // normal prior
  double sigma_prior = 10.0;
  double prior_coeff = -1.0 / (2.0 * sigma_prior);
  
  mat Sigma = computeSigma(Designmat, phi, sigmasq);
  mat InvSigma = inv(Sigma);
  
  // Designmat 준비
  mat Xth = ones(thnrow, 1);
  Xth.insert_cols(1, Designmat);
  vec logprobs(Niter);
  
  for (int k = 0; k < Niter; k++) {
    if (k > 1000 && k % 1000 == 0) {
      COV = cov(theta);
    }
    
    for(int i = 0; i< nCOVcols; i++){
      thetaprev[i] = theta(k,i);
    }
    
    
    vec Znormal = randn(nCOVcols);
    vec thetaprop = trans(  trans(thetaprev) + trans(Znormal)*chol(COV)  );
    vec h1dcross(thnrow, fill::zeros);
    for (int i = 0; i < thnrow; i++) {
      h1dcross[i] = 0.0;
      for (int j = 0; j < nCOVcols; j++) {
        double diff = fabs(thetaprop[j] - Designmat(i, j));
        h1dcross[i] += pow(diff / phi[j], 2.0); // normalize with phi
      }
      h1dcross[i] = sqrt(h1dcross[i]); // Euclidean distance
    }
    
    mat Sigmacross(thnrow, 1, fill::zeros);
    for (int i = 0; i < thnrow; i++) {
      Sigmacross(i, 0) = sigmasq * (1 + sqrt(3) * h1dcross[i]) * exp(-sqrt(3) * h1dcross[i]);
    }
    
    vec xpoint = ones(1);
    xpoint.insert_rows(1, thetaprop);
    
    double prior_prev = prior_coeff * dot(thetaprev, thetaprev);
    double prior_prop = prior_coeff * dot(thetaprop, thetaprop);
    
    lhXZp = as_scalar(trans(xpoint) * betahat + trans(Sigmacross) * InvSigma * (y - Xth * betahat));
    logprob = lhXZp - lhXZ + (prior_prop - prior_prev); 
    logprobs[k] = logprob;
    
    u = log(randu());
    if (u < logprob) {
      theta.insert_rows(k + 1, thetaprop.t());
      lhXZ = lhXZp;
    } else {
      theta.insert_rows(k + 1, thetaprev.t());
    }
  }
  
  return List::create(
    Named("theta") = theta,
    Named("logprob") = logprob,
    Named("logprobs") = logprobs
  );
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// Baseline AVM
List ergm_AVM(mat X, vec grade, vec sex, mat COV, mat theta, int outer, int cycle, int adaptInterval, double adaptFactorExponent, int adapIter, int thin){
  
  // Initializing part
  double logprob,u; 
  double rhat = 0, gamma1 = 0, gamma2 = 0;
  vec accprob = zeros(outer);
  double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234;
  int nCOVcols = COV.n_cols;
  vec thetaprev(nCOVcols);          
  vec stat = Summary(X, grade,sex), statprop(nCOVcols);
  double sigma2 = (5.76)/nCOVcols;
  mat cholCOV=COV;
  int mj = nCOVcols;
  cholCOV = trans( chol( sigma2 * (COV + 0.001 * diagmat(ones(mj)) ) )) ;
  
  //// Start of OUTER MCMC Chain 
  for(int l = 0; l< outer; l++){
    
    // adaptively update COV
    if( (l+1 >= adaptInterval) && (l+1 - (adaptInterval * trunc((l+1) / adaptInterval)) == 0) && (l < outer*0.2) ){
      double dummyacc=0;
      for(int acc=l+1-adaptInterval; acc<l; acc++){
        dummyacc = dummyacc + accprob(acc);
      }
      rhat =  dummyacc / (adaptInterval-1);
      gamma1 = 1 / pow(adapIter, c1);
      gamma2 = c0 * gamma1;
      sigma2 = exp( log(sigma2) + gamma2 * (rhat - ropt) );
      
      COV = COV + gamma1 * ( cov( theta.rows(l+1-adaptInterval, l-1) ) - COV );
      cholCOV = trans( chol( sigma2 * (COV + 0.001 * diagmat(ones(mj)) ) )) ;
      
      adapIter = adapIter + 1;
    }
    
    
    for(int i = 0; i< nCOVcols; i++){
      thetaprev[i] = theta(l,i);
    }
    
    vec Znormal = randn(nCOVcols);                                           // multivariate proposal by using Cholesky factorization
    vec thetaprop = trans(  trans(thetaprev) + trans(Znormal)*cholCOV  );  // proposed parameter
    
    // proposed auxiliary variable		
    statprop = ergm_gibbs(X, grade, sex, thetaprop, cycle);
    
    // log probability ratio to determine acceptance of Outer MCMC (normal prior)  
    vec dummy = ( -0.05*trans(thetaprop)*thetaprop + 0.05*trans(thetaprev)*thetaprev + trans(thetaprev - thetaprop)*(statprop - stat) );
    logprob = dummy[0];
    u = log( randu() );
    if( u< logprob ){
      theta.insert_rows(l+1,trans(thetaprop));
      accprob(l)=1;
    }else{
      theta.insert_rows(l+1,trans(thetaprev));
    }
    
    // Iteration
    if ((l + 1) % 100 == 0) {
      Rcpp::Rcout << "Iteration: " << (l + 1) << " completed." << std::endl;
    }
  }
  
  double acceptanceRate = mean(accprob);
  Rcpp::Rcout << "Final Acceptance Rate: " << acceptanceRate << std::endl;
  
  // Return theta and acceptance rate as a list
  return List::create(Named("theta") = theta,
                      Named("Accentance_rate") = acceptanceRate);
}

// [[Rcpp::export]]
// DA-AVM with function emulation
List ergm_DAAVM_GP(mat X, mat COV, mat theta, int outer, int cycle, 
                   const vec& grade, const vec& sex, double lhXZ, vec betahat, vec phihat, double sigmasq,
                   mat Designmat, vec y, int adaptInterval, double adaptFactorExponent, int adapIter, int thin) {
  // Initializing part
  double logprob, u;
  double rhat = 0, gamma1 = 0, gamma2 = 0;
  vec accprob = zeros(outer);
  double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234;
  int nCOVcols = COV.n_cols;
  vec thetaprev(nCOVcols);
  vec stat = Summary(X, grade, sex), statprop(nCOVcols);
  double sigma2 = (5.76) / nCOVcols;
  mat cholCOV = trans(chol(sigma2 * (COV + 0.001 * diagmat(ones(nCOVcols)))));
  int gp_accepted_count = 0;
  int final_accepted_count = 0;
  int thnrow = Designmat.n_rows;
  double sigma_prior = 10.0;
  double prior_coeff = -1.0 / (2.0 * sigma_prior);
  mat Sigma = computeSigma(Designmat, phihat, sigmasq);
  mat InvSigma = inv(Sigma);
  mat Xth = ones(thnrow, 1);
  Xth.insert_cols(1, Designmat);
  
  for (int l = 0; l < outer; l++) { // adaptively update covariance
    if ( (l+1 >= adaptInterval) && (l+1 - (adaptInterval * trunc((l+1) / adaptInterval)) == 0) && (l < outer*0.2) ){
      double dummyacc = 0;
      for (int acc = l+1-adaptInterval; acc < l; acc++) {
        dummyacc += accprob(acc);
      }
      rhat = dummyacc / (adaptInterval-1);
      gamma1 = 1 / pow(adapIter, c1);
      gamma2 = c0 * gamma1;
      sigma2 = exp(log(sigma2) + gamma2 * (rhat - ropt));
      COV = COV + gamma1 * (cov(theta.rows(l+1-adaptInterval, l-1)) - COV);
      cholCOV = trans(chol(sigma2 * (COV + 0.001 * diagmat(ones(nCOVcols)))));
      adapIter++;
    }
    
    for (int i = 0; i < nCOVcols; i++) {
      thetaprev[i] = theta(l, i);
    }
    
    vec Znormal = randn(nCOVcols);
    vec thetaprop = trans(trans(thetaprev) + trans(Znormal) * cholCOV);
    
    vec h1dcross(thnrow, fill::zeros);
    for (int i = 0; i < thnrow; i++) {
      h1dcross[i] = 0.0;
      for (int j = 0; j < nCOVcols; j++) {
        double diff = fabs(thetaprop[j] - Designmat(i, j));
        h1dcross[i] += pow(diff / phihat[j], 2.0);
      }
      h1dcross[i] = sqrt(h1dcross[i]);
    }
    
    mat Sigmacross(thnrow, 1, fill::zeros);
    for (int i = 0; i < thnrow; i++) {
      Sigmacross(i, 0) = sigmasq * (1 + sqrt(3) * h1dcross[i]) * exp(-sqrt(3) * h1dcross[i]);
    }
    
    vec xpoint = ones(1);
    xpoint.insert_rows(1, thetaprop);
    double prior_prev = prior_coeff * dot(thetaprev, thetaprev);
    double prior_prop = prior_coeff * dot(thetaprop, thetaprop);
    
    double lhXZp = as_scalar(trans(xpoint) * betahat + trans(Sigmacross) * InvSigma * (y - Xth * betahat));
    double GP_log_likelihood = lhXZp - lhXZ + (prior_prop - prior_prev);
    
    u = log(randu());
    if (u < GP_log_likelihood) {
      gp_accepted_count++;
      
      statprop = ergm_gibbs(X, grade, sex, thetaprop, cycle);
      vec dummy = (-0.05 * trans(thetaprop) * thetaprop +
        0.05 * trans(thetaprev) * thetaprev +
        trans(thetaprev - thetaprop) * (statprop - stat))-GP_log_likelihood;
      logprob = dummy[0];
      
      u = log(randu());
      if (u < logprob) {
        theta.insert_rows(l + 1, trans(thetaprop));
        final_accepted_count++;
        lhXZ = lhXZp;
      } else {
        theta.insert_rows(l + 1, trans(thetaprev));
      }
    } else {
      theta.insert_rows(l + 1, trans(thetaprev));
    }
    
    if ((l + 1) % 1000 == 0) {
      Rcpp::Rcout << "Iteration: " << (l + 1) << " completed." << std::endl;
    }
  }
  
  double gp_acceptance_rate = (double)gp_accepted_count / outer;
  double final_acceptance_rate = (double)final_accepted_count / outer;
  
  return List::create(
    Named("theta") = theta,
    Named("First_acceptance_count") = gp_accepted_count,
    Named("Final_acceptance_count") = final_accepted_count,
    Named("First_acceptance_rate") = gp_acceptance_rate,
    Named("Final_acceptance_rate") = final_acceptance_rate
  );
}

// Using normal surrogate model
// [[Rcpp::export]]
List ergm_DAAVM_F(mat X, mat COV, mat theta, int outer, int cycle, vec mu, mat covariance,
                  const vec& grade, const vec& sex, int adaptInterval, 
                  double adaptFactorExponent, int adapIter, int thin) {
  
  // Initializing part
  double logprob, u;
  double rhat = 0, gamma1 = 0, gamma2 = 0;
  vec accprob = zeros(outer);
  double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234;
  int nCOVcols = COV.n_cols;
  vec thetaprev(nCOVcols);
  vec stat = Summary(X, grade, sex), statprop(nCOVcols);
  double sigma2 = (5.76) / nCOVcols;
  
  int mj = nCOVcols;
  mat cholCOV = trans(chol(sigma2 * (COV + 0.001 * diagmat(ones(nCOVcols)))));
  int mvn_accepted_count = 0;   
  int final_accepted_count = 0; 
  
  for (int l = 0; l < outer; l++) {
    if( (l+1 >= adaptInterval) && (l+1 - (adaptInterval * trunc((l+1) / adaptInterval)) == 0) && (l < outer*0.2) ){
      double dummyacc=0;
      for(int acc=l+1-adaptInterval; acc<l; acc++){
        dummyacc = dummyacc + accprob(acc);
      }
      rhat =  dummyacc / (adaptInterval-1);
      gamma1 = 1 / pow(adapIter, c1);
      gamma2 = c0 * gamma1;
      sigma2 = exp( log(sigma2) + gamma2 * (rhat - ropt) );
      
      COV = COV + gamma1 * ( cov( theta.rows(l+1-adaptInterval, l-1) ) - COV );
      cholCOV = trans( chol( sigma2 * (COV + 0.001 * diagmat(ones(mj)) ) )) ;
      
      adapIter = adapIter + 1;
    }
    
    
    for(int i = 0; i< nCOVcols; i++){
      thetaprev[i] = theta(l,i);
    }
    
    vec Znormal = randn(nCOVcols);                                           // multivariate proposal by using Cholesky factorization
    vec theta_prop = trans(  trans(thetaprev) + trans(Znormal)*cholCOV  );  // proposed parameter
    
    // surrogate density
    double log_density_ratio = 
      -0.5 * as_scalar((theta_prop - mu).t() * inv(covariance) * (theta_prop - mu)) +
      0.5 * as_scalar((thetaprev - mu).t() * inv(covariance) * (thetaprev - mu));
    
    u = log(randu());
    
    if (u < log_density_ratio) {
      mvn_accepted_count++;  
      
      statprop = ergm_gibbs(X, grade, sex, theta_prop, cycle);
      vec dummy = (-0.05 * trans(theta_prop) * theta_prop +
        0.05 * trans(thetaprev) * thetaprev +
        trans(thetaprev - theta_prop) * (statprop - stat))-log_density_ratio;
      logprob = dummy[0];
      
      u = log(randu());
      if (u < logprob) {
        theta.insert_rows(l + 1, trans(theta_prop));
        final_accepted_count++;
      } else {
        theta.insert_rows(l + 1, trans(thetaprev));
      }
    } else {
      theta.insert_rows(l + 1, trans(thetaprev));
    }
    
    if ((l + 1) % 1000 == 0) {
      Rcpp::Rcout << "Iteration: " << (l + 1) << " completed." << std::endl;
    }
  }
  
  // Accept rates
  double mvn_acceptance_rate = (double)mvn_accepted_count / outer;
  double final_acceptance_rate = (double)final_accepted_count / outer;
  
  return List::create(
    Named("theta") = theta,
    Named("First_acceptance_count") = mvn_accepted_count,
    Named("Final_acceptance_count") = final_accepted_count,
    Named("First_acceptance_rate") = mvn_acceptance_rate,
    Named("Final_acceptance_rate") = final_acceptance_rate
  );
}
