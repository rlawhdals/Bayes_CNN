// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
// [[Rcpp::depends("RcppArmadillo")]]

#include <RcppArmadillo.h>
#include <limits>
#include <omp.h>
#include <cmath>
#include <vector>
#define min(x,y) (((x) < (y)) ? (x) : (y))

using namespace std;
using namespace Rcpp;
using namespace arma;


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// calculate r1 and r2
vec coeff(vec parameter){
	
vec r(2);
double lambda = parameter[0], t1 = parameter[1], t2 = parameter[2], t3 = parameter[3], R = 0; 
double a = (t2-R)*(t2-R)/(t1*t3*t3); 
double b = t3*t3*(t1 - 1);
double d = 27*a*b*b;
double c = pow(d + sqrt((d+2)*(d+2)-4) + 2,1.0/3.0); 
double deltar = 1/(3*b) * ( c/pow(2,1.0/3.0) + pow(2,1.0/3.0)/c + 1 );
r[0] = a / pow(deltar,1.5)+ t2;
r[1] = r[0] - sqrt(deltar);

return(r);	
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// evaluate unnormalized likelihood at proposed parameters for data.
double evalunNormX(vec thetaprop, mat initialdist){
    int nrow = initialdist.n_rows; 
	double R = 0;
	

	vec rp = coeff(thetaprop);
    double r1p = rp[0], r2p = rp[1], lambdap= thetaprop[0], t1p = thetaprop[1], t2p = thetaprop[2], t3p = thetaprop[3];
   
     
    mat resultXp(nrow,nrow);    
	for(int i = 0; i < nrow; i++){
    for(int j = 0; j <= i; j++){
    	
        if(initialdist(i,j)>100){
        resultXp(i,j) = resultXp(j,i) = log(1);
        }else if ( (initialdist(i,j)>r1p) && (initialdist(i,j)<=100) ){
        resultXp(i,j) = resultXp(j,i) = log( 1 + 1/(t3p*t3p*(initialdist(i,j)-r2p)*(initialdist(i,j)-r2p)) );
        }else if ( (initialdist(i,j)>R) && (initialdist(i,j)<=r1p) ){
        resultXp(i,j) = resultXp(j,i) = log( t1p - (sqrt(t1p)*(initialdist(i,j)-t2p)/(t2p-R))*(sqrt(t1p)*(initialdist(i,j)-t2p)/(t2p-R)) );
        }else {
        resultXp(i,j) = resultXp(j,i) = 0;
        }
        }
    }  
    
    rowvec rhoXpsum = sum(resultXp);                    
    double lhXp = 0;
	for(int i = 0; i<nrow; i++){
    	lhXp = min(rhoXpsum[i],1.2) + lhXp;
        }
    lhXp = lhXp + (nrow)*log(lambdap);              

return(lhXp);    
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// evaluate unnormalized likelihood at both current and proposed parameters for auxiliary variable. 
vec evalunNormY(vec thetaprev, vec thetaprop, mat Y){
	int count = Y.n_rows;
    double R = 0;
		
	// calculate coefficients for rho function	
	vec rp = coeff(thetaprop), r = coeff(thetaprev);
	double r1p = rp[0], r2p = rp[1], lambdap= thetaprop[0], t1p = thetaprop[1], t2p = thetaprop[2], t3p = thetaprop[3];
    double r1 = r[0], r2 = r[1], lambda= thetaprev[0], t1 = thetaprev[1], t2 = thetaprev[2], t3 = thetaprev[3];
	
	 // lhY and lhYp
    double distanceY;
    mat resultY(count,count), resultYp(count,count);
    for(int i = 0; i < count; i++){
    for(int j = 0; j <= i; j++){
    distanceY = sqrt( (Y(i,0) - Y(j,0))*(Y(i,0) - Y(j,0)) + (Y(i,1) - Y(j,1))*(Y(i,1) - Y(j,1)) );	
    
        if(distanceY>100){
        resultY(i,j) = resultY(j,i) = log(1);
        }else if ( (distanceY>r1) && (distanceY<=100) ){
        resultY(i,j) = resultY(j,i) = log( 1 + 1/(t3*t3*(distanceY-r2)*(distanceY-r2)) );
        }else if ( (distanceY>R) && (distanceY<=r1) ){
        resultY(i,j) = resultY(j,i) = log( t1 - (sqrt(t1)*(distanceY-t2)/(t2-R))*(sqrt(t1)*(distanceY-t2)/(t2-R)) );
        }else {
        resultY(i,j) = resultY(j,i) = 0;
        }
        
        if(distanceY>100){
        resultYp(i,j) = resultYp(j,i) = log(1);
        }else if ( (distanceY>r1p) && (distanceY<=100) ){
        resultYp(i,j) = resultYp(j,i) = log( 1 + 1/(t3p*t3p*(distanceY-r2p)*(distanceY-r2p)) );
        }else if ( (distanceY>R) && (distanceY<=r1p) ){
        resultYp(i,j) = resultYp(j,i) = log( t1p - (sqrt(t1p)*(distanceY-t2p)/(t2p-R))*(sqrt(t1p)*(distanceY-t2p)/(t2p-R)) );
        }else {
        resultYp(i,j) = resultYp(j,i) = 0;
        }
        
        }
    }  
    
    rowvec rhoYsum = sum(resultY), rhoYpsum = sum(resultYp);  
    double lhY = 0, lhYp = 0; 
    for(int i = 0; i<count; i++){
    	lhY = min(rhoYsum[i],1.2) + lhY;
    	lhYp = min(rhoYpsum[i],1.2) + lhYp;
        }
    lhY = lhY + (count)*log(lambda);   
    lhYp = lhYp + (count)*log(lambdap);  

    vec result(2);
    result[0] = lhY, result[1] = lhYp;

return(result);	
}


// [[Rcpp::export]]
// Birth-Death MCMC
mat BDmcmc(vec center, double range, mat initial, mat initialdist, vec parameter, int n){
	

vec r = coeff(parameter);
double r1 = r[0], r2 = r[1], lambda= parameter[0], t1 = parameter[1], t2 = parameter[2], t3 = parameter[3], R = 0;
              
int nrow = initialdist.n_rows;       
mat result(nrow,nrow);               
rowvec rho0sum(nrow);                
rowvec xprop(2);                     
double logratio,logv;               

for(int i = 0; i < nrow; i++){
    for(int j = 0; j <= i; j++){
    	
        if(initialdist(i,j)>100){
        result(i,j) = result(j,i) = log(1);
        }else if ( (initialdist(i,j)>r1) && (initialdist(i,j)<=100) ){
        result(i,j) = result(j,i) = log( 1 + 1/(t3*t3*(initialdist(i,j)-r2)*(initialdist(i,j)-r2)) );
        }else if ( (initialdist(i,j)>R) && (initialdist(i,j)<=r1) ){
        result(i,j) = result(j,i) = log( t1 - (sqrt(t1)*(initialdist(i,j)-t2)/(t2-R))*(sqrt(t1)*(initialdist(i,j)-t2)/(t2-R)) );
        }
        
        }
    } 
        	
rho0sum = sum(result);                       
double distance;                                
double likelihood,likelihoodprev;           
mat loc = initial;                          

		
double p1 = 0.5, p2 = 0.5;    
for(int k = 0; k < n; k++){

    int count = loc.n_rows;
    likelihoodprev = 0; 
    for(int i = 0; i<count; i++){
  	    likelihoodprev = min(rho0sum[i],1.2) + likelihoodprev;   
   	    }  	
   	    
	double u = randu();
	if( u < p1 ){  

	double radius = sqrt(   pow((range),2)*randu()   );
	double t = 2*(M_PI)*randu();
	double xx = center[0]+radius*cos(t), yy = center[1]+radius*sin(t);
	xprop[0] = xx, xprop[1] = yy;
	
	rowvec logrb = zeros<rowvec>(count);
	for(int i = 0; i< count; i++){
	   	distance = sqrt( (xprop[0] - loc(i,0))*(xprop[0] - loc(i,0)) + (xprop[1] - loc(i,1))*(xprop[1] - loc(i,1)) );
	   	
        if(distance>100){
        logrb[i] = log(1);
        }else if ( (distance>r1) && (distance<=100) ){
        logrb[i] = log( 1 + 1/(t3*t3*(distance-r2)*(distance-r2)) );
        }else if ( (distance>R) && (distance<=r1) ){
        logrb[i] = log( t1 - (sqrt(t1)*(distance-t2)/(t2-R))*(sqrt(t1)*(distance-t2)/(t2-R)) );
        }	  	
	}
	
	rowvec rhosum = rho0sum + logrb;
    rhosum.insert_cols(count,1);
    rhosum[count] = sum(trans(logrb));
    
    likelihood = 0; 
	for(int i = 0; i<count+1; i++){
    	likelihood = min(rhosum[i],1.2) + likelihood;
    }
    
    logratio = ( log(0.5) + log(lambda) + likelihood + log( M_PI*pow((range),2) ) ) - 
               ( log(0.5) + likelihoodprev + log(count+1) ) ;
    
    logv = log( randu() );
    if( logv < logratio ){
	loc.insert_rows(count,xprop);
    rho0sum = zeros<rowvec>(count+1);
	rho0sum = rhosum;  			
	}
	
	}else{         

	int Death = count*randu();
	xprop[0] = loc(Death,0), xprop[1] = loc(Death,1);
	
	rowvec logrd = zeros<rowvec>(count);
	for(int i = 0; i< count; i++){
	   	distance = sqrt( (xprop[0] - loc(i,0))*(xprop[0] - loc(i,0)) + (xprop[1] - loc(i,1))*(xprop[1] - loc(i,1)) );
	   	
        if(distance>100){
        logrd[i] = log(1);
        }else if ( (distance>r1) && (distance<=100) ){
        logrd[i] = log( 1 + 1/(t3*t3*(distance-r2)*(distance-r2)) );
        }else if ( (distance>R) && (distance<=r1) ){
        logrd[i] = log( t1 - (sqrt(t1)*(distance-t2)/(t2-R))*(sqrt(t1)*(distance-t2)/(t2-R)) );
        }  	
	}
	
	rowvec rhosum = rho0sum - logrd;
	rhosum.shed_col(Death);
	
	likelihood = 0; 
	for(int i = 0; i<count-1; i++){
    	likelihood = min(rhosum[i],1.2) + likelihood;
    }
   
    logratio = ( log(0.5) + likelihood + log(count) ) - 
               ( log(0.5) + log(lambda) + likelihoodprev + log( M_PI*pow((range),2) ) );
			   	
	logv = log( randu() );
    if( logv < logratio ){
    loc.shed_row(Death);
	rho0sum = zeros<rowvec>(count-1);		
   	rho0sum = rhosum;  			
	}			
	}
	
  }
  
return(loc);        	       	
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// Birth Death MCMC generating auxiliary variable for given parameter and 
// returns unnormalized likelihood at current and proposed parmeters  
vec BDmcmcLike(vec center, double range, mat initial, mat initialdist, vec thetaprev, vec thetaprop, int n, double numshard, double shard){
	
    // Initializing part 
    vec r = coeff(thetaprop);
    double r1 = r[0], r2 = r[1], lambda= thetaprop[0], t1 = thetaprop[1], t2 = thetaprop[2], t3 = thetaprop[3], R = 0;
                  
    int nrow = initialdist.n_rows;       // initial Data's distance matrix row numbers
    mat result(nrow,nrow);               // rho matrix result
    rowvec rho0sum(nrow);                // row sum of rho matrix
    rowvec xprop(2);                     // proposed points in Birth Death MCMC
    double logratio,logv;                // log acceptance ratio used in Birth Death MCMC
    
    for(int i = 0; i < nrow; i++){
        for(int j = 0; j <= i; j++){
            
            if(initialdist(i,j)>100){
            result(i,j) = result(j,i) = log(1);
            }else if ( (initialdist(i,j)>r1) && (initialdist(i,j)<=100) ){
            result(i,j) = result(j,i) = log( 1 + 1/(t3*t3*(initialdist(i,j)-r2)*(initialdist(i,j)-r2)) );
            }else if ( (initialdist(i,j)>R) && (initialdist(i,j)<=r1) ){
            result(i,j) = result(j,i) = log( t1 - (sqrt(t1)*(initialdist(i,j)-t2)/(t2-R))*(sqrt(t1)*(initialdist(i,j)-t2)/(t2-R)) );
            }
            
            }
        } 
                
    rho0sum = sum(result);                       // previous rho sum 
    double distance;                                
    double likelihood,likelihoodprev;            // likelihood (exponential part)
    mat loc = initial;                           // will be updated point location
    
    // Birth Death MCMC part 			
    double p1 = 0.5, p2 = 0.5;    // Birth and Death probability 
    for(int k = 0; k < n; k++){
    
        int count = loc.n_rows;
        likelihoodprev = 0; 
        for(int i = 0; i<count; i++){
              likelihoodprev = min(rho0sum[i],1.2) + likelihoodprev;   
               }  	
               
        double u = randu();
        if( u < p1 ){  // Birt step 
        // uniformly proposed points around circle
        double radius = sqrt(   pow((range),2)*randu()   );
        double t = (1/numshard)*2*(M_PI)*randu();
        double xx = center[0]+radius*cos(t + shard*2*(M_PI)/numshard), yy = center[1]+radius*sin(t + shard*2*(M_PI)/numshard);
        xprop[0] = xx, xprop[1] = yy;
        
        rowvec logrb = zeros<rowvec>(count);
        for(int i = 0; i< count; i++){
               distance = sqrt( (xprop[0] - loc(i,0))*(xprop[0] - loc(i,0)) + (xprop[1] - loc(i,1))*(xprop[1] - loc(i,1)) );
               
            if(distance>100){
            logrb[i] = log(1);
            }else if ( (distance>r1) && (distance<=100) ){
            logrb[i] = log( 1 + 1/(t3*t3*(distance-r2)*(distance-r2)) );
            }else if ( (distance>R) && (distance<=r1) ){
            logrb[i] = log( t1 - (sqrt(t1)*(distance-t2)/(t2-R))*(sqrt(t1)*(distance-t2)/(t2-R)) );
            }	  	
        }
        
        rowvec rhosum = rho0sum + logrb;
        rhosum.insert_cols(count,1);
        rhosum[count] = sum(trans(logrb));
        
        likelihood = 0; 
        for(int i = 0; i<count+1; i++){
            likelihood = min(rhosum[i],1.2) + likelihood;
        }
        
        logratio = ( log(0.5) + log(lambda) + likelihood + log( (1/numshard)*M_PI*pow((range),2) ) ) - 
                   ( log(0.5) + likelihoodprev + log(count+1) ) ;               
                   
        logv = log( randu() );
        if( logv < logratio ){
        loc.insert_rows(count,xprop);
        rho0sum = zeros<rowvec>(count+1);
        rho0sum = rhosum;  			
        }
        
        }else{         // Death step
        // uniformly proposed death
        int Death = count*randu();
        xprop[0] = loc(Death,0), xprop[1] = loc(Death,1);
        
        rowvec logrd = zeros<rowvec>(count);
        for(int i = 0; i< count; i++){
               distance = sqrt( (xprop[0] - loc(i,0))*(xprop[0] - loc(i,0)) + (xprop[1] - loc(i,1))*(xprop[1] - loc(i,1)) );
               
            if(distance>100){
            logrd[i] = log(1);
            }else if ( (distance>r1) && (distance<=100) ){
            logrd[i] = log( 1 + 1/(t3*t3*(distance-r2)*(distance-r2)) );
            }else if ( (distance>R) && (distance<=r1) ){
            logrd[i] = log( t1 - (sqrt(t1)*(distance-t2)/(t2-R))*(sqrt(t1)*(distance-t2)/(t2-R)) );
            }  	
        }
        
        rowvec rhosum = rho0sum - logrd;
        rhosum.shed_col(Death);
        
        likelihood = 0; 
        for(int i = 0; i<count-1; i++){
            likelihood = min(rhosum[i],1.2) + likelihood;
        }
       
        logratio = ( log(0.5) + likelihood + log(count) ) - 
                   ( log(0.5) + log(lambda) + likelihoodprev + log( (1/numshard)*M_PI*pow((range),2) ) );
                       
        logv = log( randu() );
        if( logv < logratio ){
        loc.shed_row(Death);
        rho0sum = zeros<rowvec>(count-1);		
           rho0sum = rhosum;  			
        }			
        }
        
      }
      
    vec likeval = evalunNormY(thetaprev, thetaprop, loc); 
      
    return(likeval);        	       	
    }

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// Birth Death MCMC generating auxiliary variable for given parameter and 
// returns unnormalized likelihood at current and proposed parmeters  
vec BDmcmcLikeIdeal(vec center, double range, mat initial, mat initialdist, vec thetaprev, vec thetaprop, int n, double numshard, double shard,double angle_number){
	
    // Initializing part 
    vec r = coeff(thetaprop);
    double r1 = r[0], r2 = r[1], lambda= thetaprop[0], t1 = thetaprop[1], t2 = thetaprop[2], t3 = thetaprop[3], R = 0;
                  
    int nrow = initialdist.n_rows;       // initial Data's distance matrix row numbers
    mat result(nrow,nrow);               // rho matrix result
    rowvec rho0sum(nrow);                // row sum of rho matrix
    rowvec xprop(2);                     // proposed points in Birth Death MCMC
    double logratio,logv;                // log acceptance ratio used in Birth Death MCMC
    
    for(int i = 0; i < nrow; i++){
        for(int j = 0; j <= i; j++){
            
            if(initialdist(i,j)>100){
            result(i,j) = result(j,i) = log(1);
            }else if ( (initialdist(i,j)>r1) && (initialdist(i,j)<=100) ){
            result(i,j) = result(j,i) = log( 1 + 1/(t3*t3*(initialdist(i,j)-r2)*(initialdist(i,j)-r2)) );
            }else if ( (initialdist(i,j)>R) && (initialdist(i,j)<=r1) ){
            result(i,j) = result(j,i) = log( t1 - (sqrt(t1)*(initialdist(i,j)-t2)/(t2-R))*(sqrt(t1)*(initialdist(i,j)-t2)/(t2-R)) );
            }
            
            }
        } 
                
    rho0sum = sum(result);                       // previous rho sum 
    double distance;                                
    double likelihood,likelihoodprev;            // likelihood (exponential part)
    mat loc = initial;                           // will be updated point location
    
    // Birth Death MCMC part 			
    double p1 = 0.5, p2 = 0.5;    // Birth and Death probability 
    for(int k = 0; k < n; k++){
    
        int count = loc.n_rows;
        likelihoodprev = 0; 
        for(int i = 0; i<count; i++){
              likelihoodprev = min(rho0sum[i],1.2) + likelihoodprev;   
               }  	
               
        double u = randu();
        if( u < p1 ){  // Birth step 
        // uniformly proposed points around circle
        double angle_width = M_PI / angle_number;


        double radius = sqrt(pow(range, 2) * randu());


        double center_angle = shard * 2 * M_PI / numshard;


        double t = center_angle - angle_width / 2 + angle_width * randu();


        double xx = center[0] + radius * cos(t);
        double yy = center[1] + radius * sin(t);
        xprop[0] = xx;
        xprop[1] = yy;
        rowvec logrb = zeros<rowvec>(count);
        for(int i = 0; i< count; i++){
               distance = sqrt( (xprop[0] - loc(i,0))*(xprop[0] - loc(i,0)) + (xprop[1] - loc(i,1))*(xprop[1] - loc(i,1)) );
               
            if(distance>100){
            logrb[i] = log(1);
            }else if ( (distance>r1) && (distance<=100) ){
            logrb[i] = log( 1 + 1/(t3*t3*(distance-r2)*(distance-r2)) );
            }else if ( (distance>R) && (distance<=r1) ){
            logrb[i] = log( t1 - (sqrt(t1)*(distance-t2)/(t2-R))*(sqrt(t1)*(distance-t2)/(t2-R)) );
            }	  	
        }
        
        rowvec rhosum = rho0sum + logrb;
        rhosum.insert_cols(count,1);
        rhosum[count] = sum(trans(logrb));
        
        likelihood = 0; 
        for(int i = 0; i<count+1; i++){
            likelihood = min(rhosum[i],1.2) + likelihood;
        }
        
        logratio = ( log(0.5) + log(lambda) + likelihood + log( (4/numshard)*M_PI*pow((range),2) ) ) - 
                   ( log(0.5) + likelihoodprev + log(count+1) ) ;               
                   
        logv = log( randu() );
        if( logv < logratio ){
        loc.insert_rows(count,xprop);
        rho0sum = zeros<rowvec>(count+1);
        rho0sum = rhosum;  			
        }
        
        }else{         // Death step
        // uniformly proposed death
        int Death = count*randu();
        xprop[0] = loc(Death,0), xprop[1] = loc(Death,1);
        
        rowvec logrd = zeros<rowvec>(count);
        for(int i = 0; i< count; i++){
               distance = sqrt( (xprop[0] - loc(i,0))*(xprop[0] - loc(i,0)) + (xprop[1] - loc(i,1))*(xprop[1] - loc(i,1)) );
               
            if(distance>100){
            logrd[i] = log(1);
            }else if ( (distance>r1) && (distance<=100) ){
            logrd[i] = log( 1 + 1/(t3*t3*(distance-r2)*(distance-r2)) );
            }else if ( (distance>R) && (distance<=r1) ){
            logrd[i] = log( t1 - (sqrt(t1)*(distance-t2)/(t2-R))*(sqrt(t1)*(distance-t2)/(t2-R)) );
            }  	
        }
        
        rowvec rhosum = rho0sum - logrd;
        rhosum.shed_col(Death);
        
        likelihood = 0; 
        for(int i = 0; i<count-1; i++){
            likelihood = min(rhosum[i],1.2) + likelihood;
        }
       
        logratio = ( log(0.5) + likelihood + log(count) ) - 
                   ( log(0.5) + log(lambda) + likelihoodprev + log( (4/numshard)*M_PI*pow((range),2) ) );
                       
        logv = log( randu() );
        if( logv < logratio ){
        loc.shed_row(Death);
        rho0sum = zeros<rowvec>(count-1);		
           rho0sum = rhosum;  			
        }			
        }
        
      }
      
    vec likeval = evalunNormY(thetaprev, thetaprop, loc); 
      
    return(likeval);        	       	
    }

// [[Rcpp::export]]
// code for standard AVM
mat DMH(vec center, double range, mat initial, mat initialdist, double lhX, int Niter, mat theta, mat COV, int n){


double R = 0, negativeInf = -std::numeric_limits<float>::infinity();;	
double logprob,u,lhXp, lhY, lhYp;                  
int nCOVcols = COV.n_cols;                         
vec thetaprev(nCOVcols);                           
int nrow = initialdist.n_rows;                       



for(int l = 0; l< Niter; l++){

	if( (l > 1000) && (l <= 10000) ){ 
	COV = cov(theta)+1e-20;
    }	

    for(int i = 0; i< nCOVcols; i++){
    	thetaprev[i] = theta(l,i);
    }
    
    vec Znormal = randn(nCOVcols);                                          
    vec thetaprop = trans(  trans(thetaprev) + trans(Znormal)*chol(COV)  );  


	if( thetaprop[0] > 6e-4 || thetaprop[0] < 2e-4 || thetaprop[1] > 2 || thetaprop[1] < 1 || thetaprop[2] > 20 || thetaprop[2] < 0  || thetaprop[3] > 1 || thetaprop[3] < 0 ){
	logprob = negativeInf;	
	}else{
		

	mat Y = BDmcmc(center, range, initial, initialdist, thetaprop,  n);
	int count = Y.n_rows;
	

	vec rp = coeff(thetaprop), r = coeff(thetaprev);
    double r1p = rp[0], r2p = rp[1], lambdap= thetaprop[0], t1p = thetaprop[1], t2p = thetaprop[2], t3p = thetaprop[3];
    double r1 = r[0], r2 = r[1], lambda= thetaprev[0], t1 = thetaprev[1], t2 = thetaprev[2], t3 = thetaprev[3];


 
    mat resultXp(nrow,nrow);    
	for(int i = 0; i < nrow; i++){
    for(int j = 0; j <= i; j++){
    	
        if(initialdist(i,j)>100){
        resultXp(i,j) = resultXp(j,i) = log(1);
        }else if ( (initialdist(i,j)>r1p) && (initialdist(i,j)<=100) ){
        resultXp(i,j) = resultXp(j,i) = log( 1 + 1/(t3p*t3p*(initialdist(i,j)-r2p)*(initialdist(i,j)-r2p)) );
        }else if ( (initialdist(i,j)>R) && (initialdist(i,j)<=r1p) ){
        resultXp(i,j) = resultXp(j,i) = log( t1p - (sqrt(t1p)*(initialdist(i,j)-t2p)/(t2p-R))*(sqrt(t1p)*(initialdist(i,j)-t2p)/(t2p-R)) );
        }else {
        resultXp(i,j) = resultXp(j,i) = 0;
        }
        }
    }  
    
    rowvec rhoXpsum = sum(resultXp);                    
    lhXp = 0;
	for(int i = 0; i<nrow; i++){
    	lhXp = min(rhoXpsum[i],1.2) + lhXp;
        }
    lhXp = lhXp + (nrow)*log(lambdap);              


    double distanceY;
    mat resultY(count,count), resultYp(count,count);
    for(int i = 0; i < count; i++){
    for(int j = 0; j <= i; j++){
    distanceY = sqrt( (Y(i,0) - Y(j,0))*(Y(i,0) - Y(j,0)) + (Y(i,1) - Y(j,1))*(Y(i,1) - Y(j,1)) );	
    
        if(distanceY>100){
        resultY(i,j) = resultY(j,i) = log(1);
        }else if ( (distanceY>r1) && (distanceY<=100) ){
        resultY(i,j) = resultY(j,i) = log( 1 + 1/(t3*t3*(distanceY-r2)*(distanceY-r2)) );
        }else if ( (distanceY>R) && (distanceY<=r1) ){
        resultY(i,j) = resultY(j,i) = log( t1 - (sqrt(t1)*(distanceY-t2)/(t2-R))*(sqrt(t1)*(distanceY-t2)/(t2-R)) );
        }else {
        resultY(i,j) = resultY(j,i) = 0;
        }
        
        if(distanceY>100){
        resultYp(i,j) = resultYp(j,i) = log(1);
        }else if ( (distanceY>r1p) && (distanceY<=100) ){
        resultYp(i,j) = resultYp(j,i) = log( 1 + 1/(t3p*t3p*(distanceY-r2p)*(distanceY-r2p)) );
        }else if ( (distanceY>R) && (distanceY<=r1p) ){
        resultYp(i,j) = resultYp(j,i) = log( t1p - (sqrt(t1p)*(distanceY-t2p)/(t2p-R))*(sqrt(t1p)*(distanceY-t2p)/(t2p-R)) );
        }else {
        resultYp(i,j) = resultYp(j,i) = 0;
        }
        
        }
    }  
    
    rowvec rhoYsum = sum(resultY), rhoYpsum = sum(resultYp);  
    lhY = 0, lhYp = 0; 
    for(int i = 0; i<count; i++){
    	lhY = min(rhoYsum[i],1.2) + lhY;
    	lhYp = min(rhoYpsum[i],1.2) + lhYp;
        }
    lhY = lhY + (count)*log(lambda);   
    lhYp = lhYp + (count)*log(lambdap);  


    logprob = ( lhXp + lhY ) - ( lhX + lhYp ) ;  
	}

    u = log( randu() );
    if( u< logprob ){
    theta.insert_rows(l+1,trans(thetaprop));
    lhX = lhXp;	
	}else{
	theta.insert_rows(l+1,trans(thetaprev));
	}
		
}

return theta;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// Delayed acceptance MCMC for attraction repulsion example
// subsample DMH is used in first stage
List DADMHIdealRandom(vec center, double range, mat initial, mat initialdist, double lhX, List initialsublist, List initialdistsublist, NumericVector lhXsublist, 
    int Niter, mat theta, mat COV, int n, double numshard,double angle_number){
    
    // Initializing part	
    double logprob1,u1,logprob2,u2,lhXp,lhY,lhYp,lhXpsub,lhYsub,lhYpsub;    // used in Outer MCMC
    int nCOVcols = COV.n_cols;                                              // number of parameters	
    int smallN = floor(4*n/numshard);
    
    vec thetaprev(nCOVcols);                                                // before propose in Outer MCMC, previous parameters
    int nrow = initialdist.n_rows;                                          // initial Data's number of points           
    int numDelay = 0;                                                       // count number of delay rejection
    double eps=1e-20;
    bool is_lhxsub_initialized = false;
    double lhXsub;
    //// Start of OUTER MCMC Chain 
    for(int l = 0; l< Niter-1; l++){
    
        // adaptive udpate of mean and cov for heavy tail approximations 
        if( (l > 1000) && (l <= 10000) ){ 
            COV =  cov(theta)+eps*eye(nCOVcols,nCOVcols);
        }
        
        for(int i = 0; i< nCOVcols; i++){
            thetaprev[i] = theta(l,i);
            
        }
        int shard_idx = floor(numshard * randu());
        mat initialsub = initialsublist[shard_idx];
        mat initialdistsub = initialdistsublist[shard_idx];
        if (is_lhxsub_initialized) {
            lhXsub = evalunNormX(thetaprev, initialdistsub);
        } else {
            lhXsub = lhXsublist[shard_idx];  // 아직 초기화 안됐으면 리스트에서 가져옴
        }
    
        
        vec Znormal = randn(nCOVcols);                                           // multivariate proposal by using Cholesky factorization
        vec thetaprop = trans(  trans(thetaprev) + trans(Znormal)*chol(COV)  );  // proposed parameter
        
        // constraints on prior space
        if( thetaprop[0] > 6e-4 || thetaprop[0] < 2e-4 || thetaprop[1] > 2 || thetaprop[1] < 1 || thetaprop[2] > 20 || thetaprop[2] < 0  || thetaprop[3] > 1 || thetaprop[3] < 0){
        theta.insert_rows(l+1,trans(thetaprev));	
        }else{
        
        // delayed rejection step 
        // subsample DMH is used as first stage 
        // generate auxiliary variable at proposed parameter and evaluate lhY and lhYp	
        vec genauxsub = BDmcmcLikeIdeal(center, range, initialsub, initialdistsub, thetaprev, thetaprop, smallN, numshard, shard_idx,angle_number);	
        lhYsub = genauxsub[0], lhYpsub = genauxsub[1];	
        
        // evaluate unnormalized likelihood for data at proposed parameter
        lhXpsub = evalunNormX(thetaprop, initialdistsub);
        
        // log probability ratio in the first stage
        logprob1 = ( lhXpsub + lhYsub ) - ( lhXsub + lhYpsub );   	    
        u1 = log( randu() );
        if( u1 > logprob1 ){
        theta.insert_rows(l+1,trans(thetaprev));
        numDelay = numDelay + 1;
        Rcpp::Rcout << "numDelay: " << numDelay << std::endl;
        }else{
    
    
        // generate auxiliary variable at proposed parameter and evaluate lhY and lhYp	
        vec genaux = BDmcmcLike(center, range, initial, initialdist, thetaprev, thetaprop, n, 1, 0);	
        lhY = genaux[0], lhYp = genaux[1];	
        
        // evaluate unnormalized likelihood for data at proposed parameter
        lhXp = evalunNormX(thetaprop, initialdist);
            
        // log probability ratio to determine acceptance of Outer MCMC
        logprob2 = -logprob1 +  ( lhXp + lhY ) - ( lhX + lhYp );   	
        u2 = log( randu() );
        if( u2 > logprob2 ){
        theta.insert_rows(l+1,trans(thetaprev));	
        }else{
        theta.insert_rows(l+1,trans(thetaprop));
        lhX = lhXp;	
        lhXsub = lhXpsub;
        if (!is_lhxsub_initialized) {
            is_lhxsub_initialized = true;
        } 
        }	
        
        }
        }
        Rcpp::Rcout << "Niter: " << l << std::endl;
    }
    //// End of Outer MCMC ////
    return List::create(
        Named("theta") = theta,
        Named("GP_reject_count") = numDelay
      );
    }

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
//  importance sampling
mat pResponse(vec center, double range, mat initial, mat initialdist, vec hatparameter, int inner, mat Designmat, int m, int num){

int thnrow = Designmat.n_rows;             
mat H(m,thnrow);                          
omp_set_num_threads(num);


vec rhat = coeff(hatparameter);
double r1 = rhat[0], r2 = rhat[1], lambda= hatparameter[0], t1 = hatparameter[1], t2 = hatparameter[2], t3 = hatparameter[3], R = 0;



int M;
#pragma omp parallel shared(Designmat) private(M)
{	
#pragma omp for schedule(static)  
for(M = 0; M < m; M++){                 

    mat loc = BDmcmc(center, range, initial, initialdist, hatparameter, inner);    
    int count = loc.n_rows;                                                      
    double locdistance;                                                     
    	   
    for(int k = 0; k < thnrow; k++){
        

        vec parameter = trans( Designmat.row( k ) );
        vec r = coeff(parameter);
        double R1 = r[0], R2 = r[1], Lambda= parameter[0], T1 = parameter[1], T2 = parameter[2], T3 = parameter[3];	
   
        double likelihood = 0;     
        double sumResult = 0;     
        double result = 0;         
        
        double Likelihood = 0;     
        double SumResult = 0;      
        double Result = 0;         
        	    	
        for(int i = 0; i < count; i++){
        	    sumResult = 0;
	            SumResult = 0;
	        
                for(int j = 0; j < count; j++){
	                result = 0;
			        Result = 0; 

		
                    locdistance = sqrt( (loc(i,0)-loc(j,0))*(loc(i,0)-loc(j,0)) + (loc(i,1)-loc(j,1))*(loc(i,1)-loc(j,1)) );
            
                
                    if(locdistance>100){
                        result = log(1);
                    }else if ( (locdistance>r1) && (locdistance<=100) ){
                        result = log( 1 + 1/(t3*t3*(locdistance-r2)*(locdistance-r2)) );
                    }else if ( (locdistance>R) && (locdistance<=r1) ){
                        result = log( t1 - (sqrt(t1)*(locdistance-t2)/(t2-R))*(sqrt(t1)*(locdistance-t2)/(t2-R)) );
                    }
                
                  
                    if(locdistance>100){
                        Result  = log(1);
                    }else if ( (locdistance>R1) && (locdistance<=100) ){
                        Result  = log( 1 + 1/(T3*T3*(locdistance-R2)*(locdistance-R2)) );
                    }else if ( (locdistance>R) && (locdistance<=R1) ){
                        Result  = log( T1 - (sqrt(T1)*(locdistance-T2)/(T2-R))*(sqrt(T1)*(locdistance-T2)/(T2-R)) );
                    }
                                         
                    sumResult = sumResult + result ;
                    SumResult = SumResult + Result ;
                    }   
				    likelihood  = likelihood  + min(sumResult,1.2);    
                    Likelihood = Likelihood + min(SumResult,1.2); 
    	        }


        H(M,k) = ( Likelihood + count*log(Lambda) ) - ( likelihood + count*log(lambda) );                                    
        }

   }
}
return(H);        	
}




// [[Rcpp::export]]
// function emulator
List GPmcmcL(int Niter, mat theta, mat COV, double lhXZ, vec betahat, vec phihat, mat Designmat, vec y, mat Distmat){
	int thnrow = Designmat.n_rows;   
                                                
	int nCOVcols = COV.n_cols;                                                    
	int nrow = Distmat.n_rows; 
                                                   
    mat result(nrow,nrow);                                                       
    rowvec rhosum(nrow);  
                                                          
	vec thetaprev(nCOVcols);                                                      
	
	double lhXZp,logprob,u;                                   
	double negativeInf = -std::numeric_limits<float>::infinity();;	               
	double phi1hat = phihat[0], phi2hat = phihat[1], phi3hat = phihat[2], phi4hat = phihat[3], sigmasqhat = phihat[4]; 
	

	mat h1(thnrow,thnrow), h2(thnrow,thnrow), h3(thnrow,thnrow), h4(thnrow,thnrow); 
    vec h1dcross(thnrow), h2dcross(thnrow), h3dcross(thnrow), h4dcross(thnrow);    
    
    for(int i = 0; i < thnrow; i++){
        for(int j = 0; j <= i; j++){
	        h1(i,j) = h1(j,i) = fabs(Designmat(i,0)-Designmat(j,0));
	        h2(i,j) = h2(j,i) = fabs(Designmat(i,1)-Designmat(j,1));
	        h3(i,j) = h3(j,i) = fabs(Designmat(i,2)-Designmat(j,2));
	        h4(i,j) = h4(j,i) = fabs(Designmat(i,3)-Designmat(j,3));
	    }
	}
	mat Sigma = sigmasqhat*(1+sqrt(3)*h1/phi1hat)%exp(-sqrt(3)*h1/phi1hat)%(1+sqrt(3)*h2/phi2hat)%exp(-sqrt(3)*h2/phi2hat)%
	                       (1+sqrt(3)*h3/phi3hat)%exp(-sqrt(3)*h3/phi3hat)%(1+sqrt(3)*h4/phi4hat)%exp(-sqrt(3)*h4/phi4hat);                      
    mat InvSigma = inv(Sigma);	    
	mat Xth = ones(thnrow,1);
    Xth.insert_cols(1,Designmat);	



	for(int k = 0; k< Niter; k++){

        vec Znormal = randn(nCOVcols); 
    	for(int i = 0; i< nCOVcols; i++){
    		thetaprev[i] = theta(k,i);
		}
	
    	vec thetaprop = trans(  trans(thetaprev) + trans(Znormal)*chol(COV)  );
     	vec r = coeff(thetaprop);
        double r1 = r[0], r2 = r[1],  lambda = thetaprop[0], t1 = thetaprop[1], t2 = thetaprop[2], t3 = thetaprop[3], R=0; 

        if( lambda > 6e-4 || lambda < 2e-4 || t1 > 2 || t1 < 1 || t2 > 20 || t2 < 0  || t3 > 1 || t3 < 0 ){
		logprob = negativeInf;	
		}else{			
		for(int i = 0; i< thnrow; i++){  
    		h1dcross[i] =  fabs(lambda-Designmat(i,0));
			h2dcross[i] =  fabs(t1-Designmat(i,1));
			h3dcross[i] =  fabs(t2-Designmat(i,2));
			h4dcross[i] =  fabs(t3-Designmat(i,3));		
		}
       
    	mat Sigmacross = sigmasqhat*(1+sqrt(3)*h1dcross/phi1hat)%exp(-sqrt(3)*h1dcross/phi1hat)%(1+sqrt(3)*h2dcross/phi2hat)%exp(-sqrt(3)*h2dcross/phi2hat)%
	    	                        (1+sqrt(3)*h3dcross/phi3hat)%exp(-sqrt(3)*h3dcross/phi3hat)%(1+sqrt(3)*h4dcross/phi4hat)%exp(-sqrt(3)*h4dcross/phi4hat);
        vec xpoint = ones(1);
    	xpoint.insert_rows(1,thetaprop);
    	lhXZp = (trans(xpoint)*betahat + trans(Sigmacross)* InvSigma*(y-Xth*betahat))[0]; 
   
    	logprob = lhXZp - lhXZ; 
        } 

        u = log( randu() );
    	if( u< logprob ){
    	theta.insert_rows(k+1,trans(thetaprop));
          lhXZ = lhXZp;		
		}else{
	    theta.insert_rows(k+1,trans(thetaprev));
		}

   }	
return List::create(
        Named("theta") = theta,
        Named("logprob") = logprob,
	    Named("L1_prop") = lhXZ
    );
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
//DA-AVM with GP emulator
List DelayedAcceptanceMCMC(
    int Niter,
    mat theta,
    mat COV,
    double L1_current,
    double lhX,
    vec betahat,
    vec phihat,
    mat Designmat,
    vec y,
    mat Distmat,
    vec center,
    double range,
    mat initial,
    mat initialdist,
    int n,
    int inner
)
{
    double negativeInf = -numeric_limits<float>::infinity();
    int gp_accepted_count = 0;
    int thnrow = Designmat.n_rows;
    for(int k=0; k<Niter; k++){
        if((k>1000) && (k<=10000)){
            COV = cov(theta);
        }
        vec thetaprev = theta.row(k).t();
        List gpResult = GPmcmcL(inner, thetaprev.t(), COV, L1_current, betahat, phihat, Designmat, y, Distmat);
        arma::mat GP_result = gpResult["theta"];
        vec thetaprop = GP_result.row(1).t();
        double logprob = gpResult["logprob"];
        double L1_prop = gpResult["L1_prop"];
        bool GP_rejected = all(thetaprop == thetaprev);
        if(!GP_rejected){
        gp_accepted_count++;
        mat Y = BDmcmc(center, range, initial, initialdist, thetaprop, n);
        int count = Y.n_rows;
        int nrow = initialdist.n_rows;
        double R = 0.0;
        vec rp = coeff(thetaprop), r = coeff(thetaprev);
        double r1p = rp[0], r2p = rp[1], lambdap= thetaprop[0], t1p = thetaprop[1], t2p = thetaprop[2], t3p = thetaprop[3];
        double r1 = r[0], r2 = r[1], lambda= thetaprev[0], t1 = thetaprev[1], t2 = thetaprev[2], t3 = thetaprev[3];
        mat resultXp(nrow,nrow);
        for(int i = 0; i < nrow; i++){
        for(int j = 0; j <= i; j++){
    	
            if(initialdist(i,j)>100){
            resultXp(i,j) = resultXp(j,i) = log(1);
            }else if ( (initialdist(i,j)>r1p) && (initialdist(i,j)<=100) ){
            resultXp(i,j) = resultXp(j,i) = log( 1 + 1/(t3p*t3p*(initialdist(i,j)-r2p)*(initialdist(i,j)-r2p)) );
            }else if ( (initialdist(i,j)>R) && (initialdist(i,j)<=r1p) ){
            resultXp(i,j) = resultXp(j,i) = log( t1p - (sqrt(t1p)*(initialdist(i,j)-t2p)/(t2p-R))*(sqrt(t1p)*(initialdist(i,j)-t2p)/(t2p-R)) );
            }else {
            resultXp(i,j) = resultXp(j,i) = 0;
            }
            }
        }  
        rowvec rhoXpsum = sum(resultXp);
        double lhXp = 0.0;
        for(int i=0; i<nrow; i++){
            lhXp += min(rhoXpsum[i], 1.2);
        }
        lhXp += nrow * log(thetaprop[0]);
        mat resultY(count, count, fill::zeros);
        mat resultYp(count, count, fill::zeros);
        double distanceY;
        double lambdaOld = thetaprev[0];
        double t1Old = thetaprev[1];
        double t2Old = thetaprev[2];
        double t3Old = thetaprev[3];
        double lambdaNew = thetaprop[0];
        double t1New = thetaprop[1];
        double t2New = thetaprop[2];
        double t3New = thetaprop[3];
        for(int i = 0; i < count; i++){
        for(int j = 0; j <= i; j++){
        distanceY = sqrt( (Y(i,0) - Y(j,0))*(Y(i,0) - Y(j,0)) + (Y(i,1) - Y(j,1))*(Y(i,1) - Y(j,1)) );	
    
            if(distanceY>100){
            resultY(i,j) = resultY(j,i) = log(1);
            }else if ( (distanceY>r1) && (distanceY<=100) ){
            resultY(i,j) = resultY(j,i) = log( 1 + 1/(t3*t3*(distanceY-r2)*(distanceY-r2)) );
            }else if ( (distanceY>R) && (distanceY<=r1) ){
            resultY(i,j) = resultY(j,i) = log( t1 - (sqrt(t1)*(distanceY-t2)/(t2-R))*(sqrt(t1)*(distanceY-t2)/(t2-R)) );
            }else {
            resultY(i,j) = resultY(j,i) = 0;
            }
        
            if(distanceY>100){
            resultYp(i,j) = resultYp(j,i) = log(1);
            }else if ( (distanceY>r1p) && (distanceY<=100) ){
            resultYp(i,j) = resultYp(j,i) = log( 1 + 1/(t3p*t3p*(distanceY-r2p)*(distanceY-r2p)) );
            }else if ( (distanceY>R) && (distanceY<=r1p) ){
            resultYp(i,j) = resultYp(j,i) = log( t1p - (sqrt(t1p)*(distanceY-t2p)/(t2p-R))*(sqrt(t1p)*(distanceY-t2p)/(t2p-R)) );
            }else {
            resultYp(i,j) = resultYp(j,i) = 0;
            }
        
            }
        }  
        rowvec rhoYsum = sum(resultY);
        rowvec rhoYpsum = sum(resultYp);
        double lhY = 0.0, lhYp_val = 0.0;
        for(int i=0; i<count; i++){
            lhY += min(rhoYsum[i], 1.2);
            lhYp_val += min(rhoYpsum[i],1.2);
        }
        lhY += count * log(lambdaOld);
        lhYp_val += count * log(lambdaNew);
        double logprob2 = ((lhXp + lhY) - (lhX + lhYp_val))-logprob;
        double u2 = log(randu());
        if(u2 < logprob2){
            theta.insert_rows(k+1, thetaprop.t());
            L1_current = L1_prop;
            lhX = lhXp;
        } else {
            theta.insert_rows(k+1, thetaprev.t());
        }
    }else{
        theta.insert_rows(k+1,thetaprev.t());
    }
    if ((k + 1) % 100 == 0) {
      Rcpp::Rcout << "Iteration: " << (k + 1) << " completed." << std::endl;
    }
  }
  return List::create(
    Named("theta") = theta,
    Named("GP_acceptance_count") = gp_accepted_count
  );
}

