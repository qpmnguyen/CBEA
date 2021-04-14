#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins("cpp11")]]

//' Get geometric mean of a data matrix with an index column
//' @param data A numeric matrix of data counts
//' @param index An integer index of set membership
//' @export
// [[Rcpp::export]]
DoubleVector getScore(NumericMatrix data, IntegerVector index) {
    int p = data.ncol();
    int size = index.size();
    double scale = sqrt(size * (p - size)/p);
    Rcpp::Function geometricmeanRow("geometricmeanRow", Rcpp::Environment::namespace_env("compositions"));
    NumericMatrix values;
    for (int i = 0; i < size; i++){
        int col_id = index[i];
        values(_,i) = data(_, col_id);
    }
    DoubleVector out = geometricmeanRow(values);
    return out;
}
