// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp; using namespace arma;



// [[Rcpp::export]]
List config_mod(arma::mat org_corr, double r, double tol, bool verbose) {
  //Rcout << "boop" << std::endl;
  vec s = arma::sum(org_corr, 0).t();
  //Rcout << "boop" << std::endl;
  mat K = arma::inv(org_corr);
  //Rcout << "boop" << std::endl;
  int dim = K.n_rows;
  vec alpha = arma::diagvec(K);
  vec beta = arma::zeros<vec>(dim);
  //Rcout << "boop" << std::endl;
  double error = 1000;
  mat temp = arma::zeros<mat>(dim,dim);
  
  mat K_est = K;
  mat C_con = arma::inv(K_est);
  vec s_con = s;
  int it = 0;
  while(error > tol){
    
    temp.each_col() = beta;
    
    K_est = temp + temp.t() + arma::diagmat(alpha);
    
    C_con = arma::inv(K_est);
    
    s_con = arma::sum(C_con, 0).t();
    
    alpha = alpha + r*(arma::diagvec(C_con) - arma::diagvec(org_corr));
    beta = beta + r* (1.0/double(dim))*(s_con - s);
    
    if(it % 1000 == 0){
      error = (arma::as_scalar(arma::accu(arma::abs((arma::diagvec(org_corr) - arma::diagvec(C_con)) / arma::diagvec(org_corr)))) + arma::as_scalar(arma::accu(arma::abs((s - s_con)/s))))/ double(dim);
      if(verbose){
      Rcout << it << std::endl;
      Rcout << error << std::endl;
      }
    }
    
    it++;
  }
  
  List toReturn = List::create(C_con, alpha, beta, it);
  

  return(toReturn);
}



