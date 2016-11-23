#include "subSampleDiagnostic.h"

Eigen::VectorXd gradientWeight(Eigen::MatrixXd grad_outer)
{
	Eigen::VectorXd w;
	
	Eigen::MatrixXd centered = grad_outer.rowwise() - grad_outer.colwise().mean();
	Eigen::MatrixXd cov = (centered.transpose() * centered) / double(grad_outer.rows() - 1);
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(cov, Eigen::EigenvaluesOnly);
	double pos_def = eig.eigenvalues().minCoeff();
  if(pos_def <= 1e-16){
    	Eigen::VectorXd D_inv_cov  = cov.diagonal().cwiseInverse();
      w  =(grad_outer.cwiseProduct(D_inv_cov.asDiagonal() * grad_outer)).colwise().sum();
  }else{
      // w =  grad_outer.cwiseProduct(cov.ldlt().solve(  grad_outer )).colwise().sum();
  }
  Eigen::MatrixXd gT = grad_outer.transpose();
  Eigen::VectorXd D_inv_cov  = cov.diagonal().cwiseInverse();
  w  = (gT.cwiseProduct(D_inv_cov.asDiagonal() * gT)).colwise().sum();
  Eigen::MatrixXd ww = (gT.cwiseProduct(D_inv_cov.asDiagonal() * gT)).colwise().sum();
  Eigen::MatrixXd w2 =  (D_inv_cov.asDiagonal() * gT).transpose();
	return w;
}


void subSampleDiag( Eigen::MatrixXd Vgrad_inner,
							   Eigen::MatrixXd grad_outer,
							   Eigen::VectorXd Ebias_inner,
							   int n,
							   int ratio
							   ){

	//Compute variance due to subsampling:
    Vgrad_inner.array() /=  double(grad_outer.rows() - 1);
    Eigen::VectorXd std_grad = Vgrad_inner.diagonal();
    std_grad.array() /= n;

    //Compute Bias due to not converged
    Ebias_inner.array() /=  double(grad_outer.rows() );
    Eigen::MatrixXd centered = grad_outer.rowwise() - grad_outer.colwise().mean();


	  //Compute MC variance
	  Eigen::MatrixXd cov = (centered.transpose() * centered) / double(grad_outer.rows() - 1);
	  Eigen::VectorXd std_grad_outer = cov.diagonal();

	  //std_grad_outer.array() /= nSubsample;
	  std_grad_outer.array() -= std_grad.array();

	  std_grad.array()       /= grad_outer.rows(); // we want the variance of the mean
	  std_grad_outer.array() /=  grad_outer.rows(); // we want the variance of the mean
	  std_grad_outer.array() *= (1-ratio); // still variance form not std

	  //compute total variance
	  Eigen::VectorXd grad_var = std_grad;
	  grad_var.array() += std_grad_outer.array();


	  std_grad_outer.array()  = std_grad_outer.array().sqrt();
	  std_grad.array()  = std_grad.array().sqrt();
	  grad_var.array()  = grad_var.array().sqrt();
	  if(0){
	      Rcpp::Rcout << "Gibbs std = " << std_grad.transpose() << "\n";
	      Rcpp::Rcout << "Outer std = " << std_grad_outer.transpose() << "\n";
	      Rcpp::Rcout << "Total std = " << grad_var.transpose() << "\n";
	      Rcpp::Rcout << "E[bias]    = " << Ebias_inner.transpose()   << "\n";
	    } else {
	      Rcpp::Rcout << "MC std = " << grad_var.sum()
	      			  << " (Gibbs = "<< std_grad.sum() 
	      			  << ", Outer = " << std_grad_outer.sum() 
                << "), bias = "  << Ebias_inner.array().abs().sum()/ grad_outer.rows()
	      			  << "\n";
	    }


}
