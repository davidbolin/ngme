#ifndef __GHMISC__H
#define __GHMISC__H
#include <Eigen/Dense>

// logNIG - propotonal to log denisty for multivariate NIG distribution
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


// logGNH- propotonal to log denisty for multivariate GH distribution
/*
    U      - (d x 1)    the values
    mu     - (d x 1)    shift parameter
    delta  - (d x 1)    location parameter
    iSigma - (d x d)    inverse scale parameter
    p      - ( double ) GIG distribution param
    a      - ( double ) GIG distribution param
    b      - ( double ) GIG distribution param
    

*/
double logGH(const Eigen::VectorXd & U,
              const Eigen::VectorXd & mu,
              const Eigen::VectorXd & delta,
              const Eigen::MatrixXd & iSigma,
              const double p,
              const double a,
              const double b );
#endif 


