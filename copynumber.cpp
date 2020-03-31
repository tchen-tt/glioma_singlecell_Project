#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat cnvs(mat& expression, int refer = 2, int bin = 4)
{
  // define characters
  int nrows, ncols;
  mat results, R, p, base, basem, diagss;
  /*
   * R = expression * p;
   * base = R.head_rows(refer);
   */
  nrows = expression.n_rows;
  ncols = expression.n_cols;
  
  p.zeros(ncols, ncols - bin);
  for(int j=0; j<ncols-bin; j++)
  {
    for(int i=j; i<j+bin+1; i++)
    {
      p(i,j) = 1;
    }
  }
  
  R = expression * p /(bin + 1);
  
  diagss.eye(R.n_rows, R.n_rows);
  diagss.diag() = stddev(R, 0, 1);
  
  //R= inv(diagss) * (R - mean(R, 1) * ones(1, R.n_cols));
  R = R - mean(R, 1) * ones(1, R.n_cols);
  
  
  
  base = R.head_rows(refer);
  
  basem.zeros(2, base.n_cols);
  for(int j=0; j<ncols-bin; j++)
  {
    basem(0, j) = base.col(j).max();
    basem(1, j) = base.col(j).min();
  }
  
  
  
  //results.zeros(R.n_elem);
  results.zeros(R.n_rows, R.n_cols);
  for(int i=0; i<nrows; i++)
  {
    for(int j=0; j<ncols-bin; j++)
    {
      if(R(i,j) > basem(0,j) + 0.2)
      {
        results(i,j) = R(i,j) - basem(0,j);
      } else if(R(i,j) < basem(1, j) - 0.2)
      {
        results(i,j) = R(i,j) - basem(1,j);
      } else {results(i,j) = 0.0;}
    }
  }
  return results;
}



