#include "measError.h"
#include "error_check.h"



void nsGaussianMeasurementError::printIter()
{
  Rcpp::Rcout << "theta = " << theta;
  
}
void nsGaussianMeasurementError::setupStoreTracj(const int Niter) // setups to store the tracjetory
{
  theta_vec.resize(Niter, B[0].cols());
  vec_counter = 0;
  store_param = 1;
}


nsGaussianMeasurementError::nsGaussianMeasurementError(){
  counter = 0;
  EV  = 1.;  // if there the random variance in the Noise E[V]
  EiV = 1.;
  noise = "nsNormal";
  npars = 1;
  store_param = 0;
}

Rcpp::List nsGaussianMeasurementError::toList()
{
  out_list["theta"]  = theta;
  out_list["gsd"]     = sigma;
  out_list["noise"]       = noise;
  if(store_param){
    out_list["theta_vec"] = theta_vec;
  }
  
  return(out_list);
}

void nsGaussianMeasurementError::initFromList(Rcpp::List const &init_list)
{
  
  
  Rcpp::List B_list = init_list["B"];
  nrep = B_list.length();
  B.resize(nrep);
  int count =0;
  for( Rcpp::List::iterator it = B_list.begin(); it != B_list.end(); ++it ) {
    B[count++] = Rcpp::as < Eigen::MatrixXd >( it[0]);
  }
  npars = B[0].cols();
  if(init_list.containsElementNamed("theta"))
    theta = Rcpp::as < Eigen::VectorXd >( init_list["theta"]);
  else
    theta.resize(npars);
  
  sigma.resize(nrep);
  for(int i=0;i<nrep;i++){
    sigma[i] = B[i]*theta;
    sigma[i] = sigma[i].array().exp();
  }
  
}

void nsGaussianMeasurementError::gradient(const int i,
                                        const Eigen::VectorXd& res,
                                        const double weight)
{
  counter++;
  VectorXd sigma_res, sigma_res2,sigma2;
  for(int j =0; j < npars; j++){
    sigma2 = sigma[i].array().pow(2);
    sigma_res = res.cwiseQuotient(sigma2); 
    sigma_res = sigma_res.cwiseProduct(res);
    sigma_res = sigma_res.cwiseProduct(B[i].col(j));
    dtheta[j] += 0.5*weight *(-B[i].col(j).sum() + sigma_res.sum());  
    for(int k =j; k < npars; k++){
      sigma_res2 = sigma_res.cwiseProduct(B[i].col(k));
      ddtheta(j,k) +=  -weight * sigma_res2.sum();  
      if(k>j){
        ddtheta(k,j) = ddtheta(j,k);
      }
    }
  }
}

void nsGaussianMeasurementError::step_theta(const double stepsize,
                                          const double learning_rate,
                                          const double polyak_rate,
                                          const int burnin)
{
  
  
  step  = ddtheta.ldlt().solve(dtheta);
  dtheta_old.array() *= learning_rate;
  dtheta_old += step;
  theta  += stepsize * dtheta_old;
  
  for(int i=0;i<nrep;i++){
    sigma[i] = B[i]*theta;
    sigma[i] = sigma[i].array().exp();
  }
  
  clear_gradient();
  counter = 0;
  ddtheta.setZero(npars,npars);
  if(store_param){
    if(vec_counter == 0 || polyak_rate == -1){
      theta_vec.col(vec_counter) = theta;
    } else{
      theta_vec.col(vec_counter) = polyak_rate * theta;
      step = theta_vec.col(vec_counter-1);
      step *= (1 - polyak_rate);
      theta_vec.col(vec_counter) += step;
    }
    vec_counter++;
  }
}


std::vector< Eigen::VectorXd > nsGaussianMeasurementError::simulate(std::vector< Eigen::VectorXd > Y)
{
  std::vector< Eigen::VectorXd > residual( Y.size());
  for(int i = 0; i < Y.size(); i++){
    residual[i] =  sigma[i].cwiseProduct(Rcpp::as< Eigen::VectorXd >(Rcpp::rnorm( Y[i].size()) ));
  }
  return(residual);
}


Eigen::VectorXd  nsGaussianMeasurementError::simulate(const Eigen::VectorXd & Y)
{
  Rcpp::Rcout << "Warning, must specify patient number to know what sigma is\n";
  Eigen::VectorXd residual =  (Rcpp::as< Eigen::VectorXd >(Rcpp::rnorm( Y.size()) ));
  return(residual);
}


Eigen::VectorXd  nsGaussianMeasurementError::simulate_par(const Eigen::VectorXd & Y,std::mt19937 & random_engine)
{
  Rcpp::Rcout << "Warning, must specify patient number to know what sigma is\n";
  std::normal_distribution<double> normal;
  Eigen::VectorXd residual;
  residual.setZero(Y.size());
  return(residual);
}

void nsGaussianMeasurementError::clear_gradient()
{
  dtheta.setZero(npars);
}

Eigen::VectorXd nsGaussianMeasurementError::get_gradient()
{
  Eigen::VectorXd g(npars);
  g = dtheta;
  return(g);
}