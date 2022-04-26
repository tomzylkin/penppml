// [[Rcpp::depends(RcppEigen)]]

// #include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>
#include "penppml_types.h"

using namespace Rcpp;
using namespace Eigen;

//' Faster Matrix Multiplication
//'
//' Faster matrix multiplication using C++.
//'
//' @param A,B Matrices.
// [[Rcpp::export]]
SEXP eigenMatMult(Eigen::MatrixXd A, Eigen::MatrixXd B){
  Eigen::MatrixXd C;
  C = A * B;

  return wrap(C);
}

//' @rdname eigenMatMult
// [[Rcpp::export]]
SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd C;
  C = A * B;

  return Rcpp::wrap(C);
}

//' Computing A'A
//'
//' Computes A'A using C++.
//'
//' @param A A matrix.
// [[Rcpp::export]]
Eigen::MatrixXd AtA(const MapMatd& A) {
  int  n(A.cols());
  return  Eigen::MatrixXd(n,n).setZero().selfadjointView<Eigen::Lower>().rankUpdate(A.adjoint());
}


//' Faster Least Squares Estimation
//'
//' Finds Least Squares solutions using C++.
//'
//' @param X Regressor matrix.
//' @param y Dependent variable (a vector).
//'
//' @return The vector of parameter (beta) estimates.
// [[Rcpp::export]]
NumericVector fastolsCpp(MapMatd X, MapVecd y){
  Eigen::MatrixXd A=AtA(X);
  const Eigen::LLT<Eigen::MatrixXd> llt(A);
  const Eigen::VectorXd betahat(llt.solve(X.adjoint() * y));
  return wrap(betahat);
}

//' Faster Ridge Regression
//'
//' Finds Ridge solutions using C++.
//'
//' @param X Regressor matrix.
//' @param y Dependent variable (a vector).
//' @param lambda Penalty parameter (a number).
//'
//' @return The vector of parameter (beta) estimates.
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

//' Faster Standard Deviation
//'
//' Computes standard deviation using C++.
//'
//' @param X Regressor matrix.
//' @param w Weights.
//'
//' @return Vector of standard deviations of the parameter estimates.
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

//' Faster Weighted Mean
//'
//' Computes weighted mean using C++.
//'
//' @param X Regressor matrix.
//' @param w Weights.
//'
//' @return Weighted mean.
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

//' Many Outer Products
//'
//' Compute a large number of outer products (useful for clustered SEs) using C++.
//'
//' @param A,B Numeric vectors.
//' @param c Integer.
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

MatrixXd selfouter(Eigen::VectorXd A){
  int r = A.rows();
  MatrixXd M = A*A.transpose();
  return M;
}

//' XeeX Matrix Computation
//'
//' Given matrix ee' and matrix X, compute X(k)'ee'X(k) for each regressor X.
//'
//' @param X Regressor matrix.
//' @param e Residuals.
//' @param S Cluster sizes.
//'
//' @return The matrix product X(k)'ee'X(k).
// [[Rcpp::export]]
SEXP xeex(const Eigen::MatrixXd X,  const Eigen::VectorXd e, const Eigen::VectorXd S){
  int K  = X.cols();   //  A is x's arranged as g (TxK) matrices vertically stacked
  int g = S.size();    // number of clusters
  int s = 0;           // size of current cluster
  int p = 0;           // placekeeper for where to start
  MatrixXd M(K,K);
  M = Eigen::MatrixXd::Zero(K, K);
  for (int n = 0; n <= g-1; n++) {
    s = S(n);
    M = M + selfouter(X.block(p,0,s,K).transpose()*e.segment(p,s)); // https://eigen.tuxfamily.org/dox/group__TutorialMatrixClass.html
    p = p + s;
  }
  return Rcpp::wrap(M);
}
