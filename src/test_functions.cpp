#include <Rcpp.h>
#include <RcppEigen.h>
#include "sample.h"
#include "solver.h"

// [[Rcpp::export]]
double test_PreDiagsolver(Rcpp::List in_list)
{
	Eigen::SparseMatrix<double, 0, int> Q = Rcpp::as<Eigen::SparseMatrix<double,0,int> >(in_list["Q"]);
	Eigen::VectorXd z = Rcpp::as<Eigen::VectorXd > (in_list["z"]);
	Eigen::VectorXd b = Rcpp::as<Eigen::VectorXd > (in_list["b"]);
	cholesky_solver Solver;
	Solver.init(Q.rows(), 0, 0, 0);
	Solver.analyze(Q);
  	Solver.compute(Q);
  	
    Eigen::VectorXd res = Solver.rMVN(b, z);
    
    Eigen::SparseMatrix<double, 0, int> Q2(Q);
    Eigen::VectorXd D_12 = Q.diagonal().cwiseInverse().cwiseSqrt();
    Q2 = D_12.asDiagonal() * Q * D_12.asDiagonal();
    Solver.analyze(Q2);
    Solver.compute(Q2);
    Eigen::VectorXd b_ = b.cwiseProduct(D_12);
    Eigen::VectorXd res2 = Solver.rMVN(b_, z);
    res2 = res2.cwiseProduct(D_12);
    return((res - res2).sum());
}

// [[Rcpp::export]]
Eigen::VectorXi   sampleR(int n, Eigen::VectorXd w_in)
{
	Eigen::VectorXd w;
  	w = w_in;
  	w /= w.sum();
  	return(ProbSampleNoReplace(n, w));
}