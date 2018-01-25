#include "operatorMatrix.h"
#include "error_check.h"

using namespace std;


double fd2Operator::trace_variance( const Eigen::SparseMatrix<double,0,int>& A, int i){

	if( nop == 1)
		i = 0;

	if(1){
	  Eigen::VectorXd  obs_loc = A * loc[i].tail(loc[i].size()-1); //Adjust for dirichlet BC
	  obs_loc.array() -= m_loc[i];
	  // covariance / h
	  // covariance = t^3/3 (integrated brownian moition
	  return(obs_loc.array().pow(3).sum() / (3 * h_average[i] * tau));
	} else {
	  Eigen::SparseMatrix<double,0,int> Qt = Q[i].transpose();
	  Eigen::SparseMatrix<double,0,int> At = A.transpose();
	  Eigen::SparseLU< Eigen::SparseMatrix<double,0,int> > chol(Qt);
	  Eigen::MatrixXd QtAt = chol.solve(At);
	  Eigen::MatrixXd QQ = QtAt*h[i].asDiagonal()*QtAt.transpose();
	  return(QQ.trace());
	}


}
