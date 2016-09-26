#include <Rcpp.h>
#include "MixedEffect.h"

Eigen::VectorXd sample_Nc(const Eigen::VectorXd & b,const  Eigen::MatrixXd & Q) {
  Eigen::VectorXd Z = Rcpp::as< Eigen::VectorXd >(Rcpp::rnorm( b.size()) );
  Eigen::LLT<Eigen::MatrixXd> lltOfQ(Q);
  Eigen::VectorXd U = lltOfQ.matrixL().solve(b);
  return lltOfQ.matrixU().solve(U + Z);
}

Eigen::VectorXd sample_Nc_par(const Eigen::VectorXd & b,const  Eigen::MatrixXd & Q, std::mt19937 & random_engine)
{
  std::normal_distribution<double> normal;
  Eigen::VectorXd Z;
  Z.setZero(b.size());
  for(int j =0; j < b.size(); j++)
    Z[j] =  normal(random_engine);

  Eigen::LLT<Eigen::MatrixXd> lltOfQ(Q);
  Eigen::VectorXd U = lltOfQ.matrixL().solve(b);
  return lltOfQ.matrixU().solve(U + Z);
}
