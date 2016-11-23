#ifndef LDMOD_SSD_H_
#define LDMOD_SSD_H_

#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Cholesky>

/*
	prints the total inner and outer variance due to the subsampling
	Vgrad_inner - matrix inner covariance
	                  cov[x_{.j}]
	grad_outer  - the expectation of for each inner
						E[x_{.j}]
	n           - number of samples
 	ratio       - precentage of observations sampled
*/
void subSampleDiag( Eigen::MatrixXd ,
							   Eigen::MatrixXd ,
							   Eigen::VectorXd ,
							   int ,
							   int 
							   );

/*
	Computing w_i = \grad_i^T Cov^{-1} \grad_i
	grad_outer  - the expectation of for each inner
						E[x_{.j}]

*/
Eigen::VectorXd gradientWeight(Eigen::MatrixXd );

#endif 