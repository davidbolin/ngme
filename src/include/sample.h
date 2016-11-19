#ifndef LDMOD_SAMPLE_H_
#define LDMOD_SAMPLE_H_

#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector>
// sampling without repleacement
// nans - number of samples
// p    - (n x 1) probabilility of sampling
Eigen::VectorXi ProbSampleNoReplace( int nans, Eigen::VectorXd p);

#endif 