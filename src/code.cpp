// [[Rcpp::depends(RcppEigen)]]

// #include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>
#include "penppml_types.h"

using namespace Rcpp;
using namespace Eigen;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
//   more tips: https://stackoverflow.com/questions/13956292/rcpp-inline-creating-and-calling-additional-
//
//   sparse matrices? https://eigen.tuxfamily.org/dox/group__SparseQuickRefPage.html

// [[Rcpp::export]]
SEXP eigenMatMult(Eigen::MatrixXd A, Eigen::MatrixXd B){
  Eigen::MatrixXd C;
  C = A * B;

  return wrap(C);
}


// [[Rcpp::export]]
SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd C;
  C = A * B;

  return Rcpp::wrap(C);
}

Eigen::MatrixXd AtA(const MapMatd& A) {
  int  n(A.cols());
  return  Eigen::MatrixXd(n,n).setZero().selfadjointView<Eigen::Lower>().rankUpdate(A.adjoint());
}

// [[Rcpp::export]]
NumericVector fastolsCpp(MapMatd X, MapVecd y){
  Eigen::MatrixXd A=AtA(X);
  const Eigen::LLT<Eigen::MatrixXd> llt(A);
  const Eigen::VectorXd betahat(llt.solve(X.adjoint() * y));
  return wrap(betahat);
}

// [[Rcpp::export]]
NumericVector fastridgeCpp(MapMatd X, MapVecd y, double lambda){
  Eigen::MatrixXd A=AtA(X);
  int k = A.cols();
  Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(k, k);
  A=A+lambda*Id;
  const Eigen::LLT<Eigen::MatrixXd> llt(A);
  const Eigen::VectorXd betahat(llt.solve(X.adjoint() * y));
  return wrap(betahat);
}

// [[Rcpp::export]]
NumericVector faststddev(Eigen::ArrayXXd X, Eigen::ArrayXd w){
  w = w/(w.sum());
  Eigen::ArrayXXd Xw      = X.colwise()*w;
  Array<double, Dynamic, Dynamic> step1 = (((X.rowwise() - Xw.colwise().sum()).square()).colwise()*w);
  Array<double, Dynamic, 1> step2 = step1.colwise().sum();
  Array<double, Dynamic, 1> step3 = step2.sqrt();
  //Array<double, 1, Dynamic> std_dev = (w*((X.rowwise() - Xw.colwise().sum()).square().colwise()).colwise().sum()).sqrt();
  return wrap(step3);
}

// [[Rcpp::export]]
NumericVector fastwmean(Eigen::ArrayXXd X, Eigen::ArrayXd w){
  //Eigen::ArrayXXd Xw = X*w;
  // convert weights to sum to 1?
  int n = X.rows();
  int k = X.cols();
  w = w/(w.sum());
  Eigen::ArrayXXd Xw = X.colwise()*w;
  Array<double, Dynamic, 1> wmean = Xw.colwise().sum();
  return wrap(wmean);
}


// h/t: https://stackoverflow.com/questions/49206780/column-wise-initialization-and-calculation-of-standard-deviation-in-eigen-librar

// many outer products
// [[Rcpp::export]]
SEXP manyouter(const Eigen::Map<Eigen::VectorXd> A, Eigen::Map<Eigen::VectorXd> B, int c){
  int r = A.rows();
  Eigen::MatrixXd M(r,c);
  int g = r/c;
  for (int n = 0; n <= g-1; n++) {
    M.block(n*c,0,c,c) =  A.segment(n*c,c)*B.segment(n*c,c).transpose();
  }
  return Rcpp::wrap(M);
}

// Block of size (p,q), starting at (i,j): matrix.block(i,j,p,q);
// Vector segment containing n elements, starting at position i: vector.segment(i,n);

// [[Rcpp::export]]
SEXP xeex(const Eigen::Map<Eigen::MatrixXd> X, Eigen::Map<Eigen::MatrixXd> ee){
  int K  = X.cols();  //  A is x's arranged as g (TxK) matrices vertically stacked
  int gT = ee.rows();  //  B is g TxT blocks arranged vertically
  int T  = ee.cols();
  int g = gT/T;
  Eigen::MatrixXd M(K,K);
  M = Eigen::MatrixXd::Zero(K, K);
  for (int n = 0; n <= g-1; n++) {
    M = M + X.block(n*T,0,T,K).transpose()*ee.block(n*T,0,T,T)*X.block(n*T,0,T,K);
  }
  return Rcpp::wrap(M);
}
