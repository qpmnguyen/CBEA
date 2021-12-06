#include <Rcpp.h>
#include <Rmath.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::plugins("cpp11")]]

//' Geometric mean of a vector
//' @param vec A vector of values with length \code{n}
//' @description Compute geometric mean of a vector
//'     using \code{exp(mean(log(.x)))} format
//' @return A numeric value of the
//'     geometric mean of the vector \code{vec}
//' @examples
//' ex <- abs(rnorm(10))
//' gmean(ex)
//' @export
// [[Rcpp::export]]
double gmean(NumericVector vec){
    double output = exp(mean(log(vec)));
    return output;
}

//' Geometric mean of rows of a matrix
//' @param X A numeric matrix with \code{n} rows and \code{p} columns
//' @description This function computes the geometric
//'     mean by row of a numeric matrix
//' @return A numeric vector of the geometric
//'     mean of the matrix \code{X} with length \code{n}
//' @examples
//' ex <- matrix(rnorm(100), nrow = 10, ncol = 10)
//' ex <- abs(ex)
//' gmeanRow(ex)
//' @export
// [[Rcpp::export]]
DoubleVector gmeanRow(NumericMatrix X){
    if (std::any_of(X.begin(), X.end(), [](double i){return i <= 0;})){
        stop("X has 0 or negative values");
    }
    R_xlen_t nrow = X.nrow();
    DoubleVector output (nrow);
    for (R_xlen_t i = 0; i < nrow; i++){
        output[i] = gmean(X(i,_));
    }
    return output;
}


