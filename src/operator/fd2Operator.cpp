#include "operatorMatrix.h"
#include "error_check.h"

using namespace std;

double fd2Operator::trace_variance( const Eigen::SparseMatrix<double,0,int>& A, int i){

	if( nop == 1)
		i = 0;
  Eigen::VectorXd  obs_loc = A * loc[i];
  obs_loc.array() -= m_loc[i];
  // covariance / h
  // covariance = t^3/3 (integrated brownian moition
  return(obs_loc.array().pow(3).sum() / (3 * h_average[i] * tau));

}

Rcpp::List fd2Operator::output_list()
{
  Rcpp::List  List;
  List["fd2"] = "fd2";
  List["tau"] = tauVec(tauVec.size() -1 );
  List["tauVec"] = tauVec;
  //List["Q"] = Q;
  List["loc"] = loc;
  List["nIter"] = tauVec.size();
  List["h"] = h;
  List["Cov_theta"]   = Cov_theta;
  return(List);
}
