#include "NGIG.h"
#include <Rcpp.h>


double logNIG(const Eigen::VectorXd & U,
              const Eigen::VectorXd & mu,
              const Eigen::VectorXd & delta,
              const Eigen::MatrixXd & iSigma,
              const double nu )
{
	double p = -0.5 * ( 1 + U.size());
	Eigen::VectorXd U_delta = U - delta;
	double b = U_delta.dot( iSigma * U_delta) + nu;
	double a = mu.dot( iSigma * mu) + nu;
	double logf = U_delta.dot( iSigma * mu);
	logf -= 0.75 * log(b);
	double sqrt_ab  = sqrt(a * b);

	 double K1 = R::bessel_k(sqrt_ab, p, 2);
	 logf += log(K1) - sqrt_ab;
	 return(logf); 
}