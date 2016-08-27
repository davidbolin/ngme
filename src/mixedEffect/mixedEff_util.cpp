#include <Rcpp.h>
#include "MixedEffect.h"

Eigen::VectorXd sample_Nc(const Eigen::VectorXd & b,const  Eigen::MatrixXd & Q) {
  Eigen::VectorXd Z = Rcpp::as< Eigen::VectorXd >(Rcpp::rnorm( b.size()) );
  Eigen::LLT<Eigen::MatrixXd> lltOfQ(Q);
  Eigen::VectorXd U = lltOfQ.matrixL().solve(b);
  return lltOfQ.matrixU().solve(U + Z);
}
