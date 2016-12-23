#include <Rcpp.h>
#include "MixedEffect.h"
#include "GIG.h"



void dU_ddU_NIG(
              Eigen::VectorXd & dU,
              Eigen::MatrixXd & ddU,
        const Eigen::VectorXd & U,
        const Eigen::MatrixXd & Sigma,
        const Eigen::VectorXd & delta,
        const Eigen::VectorXd & mu,
        const double p_GIG,
        const double a_GIG,
        const double b_GIG,
        const Eigen::VectorXd & res,
        const Eigen::MatrixXd & Q_noise,
        const Eigen::MatrixXd & B)
{
  //computing EiV 
  double p = p_GIG - 0.5 * U.size();
  Eigen::VectorXd X_delta = U - delta;
  Eigen::MatrixXd iSigma = Sigma.inverse();
  Eigen::VectorXd temp = iSigma * X_delta;
  double b = X_delta.dot( temp) + b_GIG;
  temp =  iSigma  * mu;
  double a = mu.dot(temp) + a_GIG;
  double EiV = EiV_GIG(p, a, b);

  dU.setZero(U.size());
  dU = - B.transpose() * (Q_noise * (res - B * U) );
  dU += iSigma * ( EiV * (U - delta) - mu);

  ddU.setZero(U.size(), U.size());
  ddU += B.transpose() * Q_noise * B; 
  ddU += EiV * iSigma;

  double db_EiV = db_EiV_GIG(p, a, b);
  Eigen::MatrixXd D =   2 * (iSigma * (U - delta)) * (iSigma * (U - delta)).transpose();
  D *= db_EiV;
  ddU += D;
}

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

Eigen::VectorXd dU_EiV_NGIG(const Eigen::VectorXd & U,
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
  Eigen::VectorXd Qx = llt_Sigma.solve(X_delta);
  double b = X_delta.dot( Qx) + b_GIG;
  Eigen::VectorXd temp2 = llt_Sigma.solve(mu);
  double a = mu.dot(temp2) + a_GIG;
  double db_EiV = db_EiV_GIG(p, a, b);
  Qx *= 2 * db_EiV;
  return Qx;
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
