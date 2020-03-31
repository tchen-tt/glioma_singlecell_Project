#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]

mat filtercnv(mat cnvdata, double meancnvnorm) {
  mat results(cnvdata.n_cols, 2, fill::zeros);
  cnvdata = cnvdata - meancnvnorm;
  mat sums = cnvdata.t() * cnvdata / cnvdata.n_rows;
  results.col(0) = sums.diag();
  for(int i = 0; i < cnvdata.n_cols; i++) {
    mat colsvecs = colvec(cnvdata.n_cols, fill::ones);
    colsvecs.at(i, 0) = 0;
    mat rmove = cnvdata * colsvecs / (cnvdata.n_cols - 1);
    mat t = cor(cnvdata.col(i), rmove);
    results.at(i, 1) = t.at(0, 0);
  }
  return results;
}
