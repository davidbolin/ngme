#include <Rcpp.h>
#include "MixedEffect.h"
#include "GIG.h"


double EiV_NGIG(const Eigen::VectorXd & U,
				const Eigen::MatrixXd & Sigma,
        const Eigen::VectorXd & delta,
				const Eigen::VectorXd & mu,
				const double p_GIG,
				const double a_GIG,
				const double b_GIG)
{
  double p = p_GIG - 0.5 * U.size();
  Eigen::VectorXd X_delta = U - delta;
  Eigen::LLT<Eigen::MatrixXd> llt_Sigma(Sigma);
  Eigen::VectorXd temp = llt_Sigma.solve(X_delta);
  double b = X_delta.dot( temp) + b_GIG;
  temp = llt_Sigma.solve(mu);
  double a = mu.dot(temp) + a_GIG;

  return EiV_GIG(p, a, b);
}

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
