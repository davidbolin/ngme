#ifndef __NGIG__H
#define __NGIG__H
#include <Eigen/Dense>

// loGNIG - log denisty for multivariate NIG distribution
/*
    U      - (d x 1) the values
    mu     - (d x 1)  shift parameter
    delta  - (d x 1) location parameter
    iSigma - (d x d) inverse scale parameter
    nu     - ( double ) shape parameter 
    

*/
double logNIG(const Eigen::VectorXd & U,
              const Eigen::VectorXd & mu,
              const Eigen::VectorXd & delta,
              const Eigen::MatrixXd & iSigma,
              const double nu );
#endif 