#include "GHmisc.h"
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
	
	
	double sqrt_ab  = sqrt(a * b);
	logf += p * log(sqrt_ab);

	 double K1 = R::bessel_k(sqrt_ab, p, 2);
	 logf += log(K1) - sqrt_ab;
	 return(logf); 
}

double logGH(const Eigen::VectorXd & U,
              const Eigen::VectorXd & mu,
              const Eigen::VectorXd & delta,
              const Eigen::MatrixXd & iSigma,
              const double p,
              const double a,
              const double b )
{
	double d = (double) U.size();
	Eigen::VectorXd U_delta        =  (U - delta);
	Eigen::VectorXd iSigma_U_delta = iSigma * U_delta;
	double c1  = b + U_delta.dot( iSigma_U_delta);
	c1        *= a + mu.dot( iSigma * mu);
	c1         = sqrt(c1);

	double logf = mu.dot( iSigma_U_delta);
	logf += (p - 0.5 * d ) * log(c1);
	 double K1 = R::bessel_k(c1, p - 0.5 * d, 2);
	 logf += log(K1) - c1;

	 return(logf);
}