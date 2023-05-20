#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <stdlib.h>
#include <math.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

double factorial(int x );
double power(double x,int n);
  
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec Exp_dpois(int n,  arma::vec lambdai) {
  
  arma::sp_mat lambda(lambdai);
  //x está okay
  arma::sp_mat x(lambda.n_rows,lambda.n_cols);
  
  for (arma::uword i = 0; i < lambda.n_rows; i++) {
        // x(i) = log(power((double)(lambda(i)),n));
        x(i) = power((double)(lambda(i)),n);
  }
  
    arma::sp_mat y(lambda.n_rows,lambda.n_cols)  ;
                      
  //y está okay                 
  for (arma::uword j = 0; j < lambda.n_rows; j++) {
        y(j) = exp(-lambda(j));
  }
  
  double z = factorial(n);
  
  
  arma::vec matriz_final(x%y);
  
  if (z > 0){
    
    for (arma::uword j = 0; j < lambda.n_rows; j++) {
      
          matriz_final(j) = matriz_final(j)/z;
      } 
    
  }
  return matriz_final;
  
}

  
double power(double x,int n){
  
  double y = 0;
  
  for(int i=0;i < n;i++){
    
    y = y + log(x);
    
  }
  
  return exp(y) ;
    
}

double factorial(int x ){
  
  if(x<2)
    return 1;
  
  return factorial(x-1) * x ;
  
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec Ext_dnorm(NumericVector Tbi, double mean, double sd ) {
  
  // NumericVector x(Tbi.size());
  // 
  // for(int i=0;i< Tbi.size(); i++){
  //   x(i) = dnorm(Tbi(i) , mean , sd);
  // }
  // 
  return dnorm(Tbi , mean , sd); ;
  
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec Ext_Integrando(NumericVector Tbi, double Tbeta0, double Tsigma2, double Ty ) {
  
  
  arma::vec out =  Exp_dpois(Ty,exp(Tbeta0+Tbi)) % Ext_dnorm(Tbi, 0, sqrt(Tsigma2));
  
  return(out);
  
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec Ext_Integrando_Completo(arma::vec Tbi, double Tbeta0, double Tsigma2, double Ty){
  
  arma::sp_mat lambda(exp(Tbeta0+Tbi));
  arma::sp_mat x(lambda.n_rows,lambda.n_cols);
  
  for (arma::uword i = 0; i < lambda.n_rows; i++) {
    // x(i) = log(power((double)(lambda(i)),n));
    x(i) = power((double)(lambda(i)),Ty);
  }
  
  double z = factorial(Ty);
  
  arma::vec matriz_final(x% arma::sp_mat(exp(-arma::vec(lambda))));
  
  if (z > 0){
    
    for (arma::uword j = 0; j < lambda.n_rows; j++) {
      
      matriz_final(j) = matriz_final(j)/z;
    } 
    
  }
  
  return matriz_final % arma::vec(dnorm((as<NumericVector>(wrap(Tbi))), 0, sqrt(Tsigma2)));
  
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
NumericVector Ext_Integrando_NumericVector(NumericVector Tbi, double Tbeta0, double Tsigma2, double Ty){
  
  // arma::sp_mat lambda(exp(Tbeta0+Tbi));
  // arma::sp_mat x(lambda.n_rows,lambda.n_cols);
  // 
  // for (arma::uword i = 0; i < lambda.n_rows; i++) {
  //   // x(i) = log(power((double)(lambda(i)),n));
  //   x(i) = power((double)(lambda(i)),Ty);
  // }
  // 
  // double z = factorial(Ty);
  // 
  // arma::vec matriz_final(x% arma::sp_mat(exp(-arma::vec(lambda))));
  // 
  // if (z > 0){
  //   
  //   for (arma::uword j = 0; j < lambda.n_rows; j++) {
  //     
  //     matriz_final(j) = matriz_final(j)/z;
  //   } 
  //   
  // }
  //return matriz_final % arma::vec(dnorm((as<NumericVector>(wrap(Tbi))), 0, sqrt(Tsigma2)));
  
  
  double z = factorial(Ty);
  
  NumericVector x(Tbi.size());
  NumericVector lambda(exp(Tbeta0+Tbi));
 
 if (z > 0){
   for (int i = 0;i< Tbi.size() ;i++){
    
     x[i] = (power((double)(lambda[i]),Ty) * exp(-lambda[i]))/z;
    
    
  }
    
 }
   
   return x;
  
}




// ((pow(lambda,x))*(exp(-lambda)))/factorial(x);
