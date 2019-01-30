#include "measError.h"
#include "error_check.h"



void GaussianMeasurementError::printIter()
{
	Rcpp::Rcout << "sigma = " << sigma;

}
void GaussianMeasurementError::setupStoreTracj(const int Niter) // setups to store the tracjetory
{
	sigma_vec.resize(Niter);
	vec_counter = 0;
	store_param = 1;
}


GaussianMeasurementError::GaussianMeasurementError(): MeasurementError(){
  nsSigma = 0;
  counter = 0;
  sigma   = 1;
  dsigma  = 0;
  dsigma_old = 0;
  ddsigma = 0;
  EV  = 1.;  // if there the random variance in the Noise E[V]
  EiV = 1.;
  noise = "Normal";
  npars = 1;
  store_param = 0;
}

Rcpp::List GaussianMeasurementError::toList()
{
  Rcpp::List out;
  out["sigma"]  = sigma;
  out["sd"]     = sigma;
  out["noise"]       = noise;
  out["Cov_theta"]   = Cov_theta;
  if(store_param){
  	out["sigma_vec"] = sigma_vec;
  	out["sigma"]  = sigma_vec[sigma_vec.size() - 1];
  	}

  return(out);
}

void GaussianMeasurementError::initFromList(Rcpp::List const &init_list)
{
  if(init_list.containsElementNamed("sigma"))
    sigma = Rcpp::as < double >( init_list["sigma"]);
  else
    sigma = 1.;
}

void GaussianMeasurementError::gradient(const int i,
                                 const Eigen::VectorXd& res,
                                 const double weight)
{
    counter++;
    dsigma += weight *(- res.size()/sigma + res.array().square().sum() / pow(sigma, 3) );
    // Expected fisher infromation
    // res.size()/pow(sigma, 2) - 3 * E[res.array().square().sum()] /pow(sigma, 4);
    ddsigma +=  weight * (- 2 * res.size()/pow(sigma, 2));
}
void GaussianMeasurementError::step_theta(const double stepsize,
										  const double learning_rate,
										  const double polyak_rate,
										  const int burnin)
{
  double sigma_temp = -1;
  dsigma /= ddsigma;
  dsigma_old = dsigma_old * learning_rate + dsigma;
  double stepsize_temp  = stepsize;
  while(sigma_temp < 0)
  {
    sigma_temp = sigma - stepsize_temp * dsigma_old;
    stepsize_temp *= 0.5;
    if(stepsize_temp <= 1e-16)
        throw("in GaussianMeasurementError:: can't make sigma it positive \n");
  }
  sigma = sigma_temp;
  clear_gradient();
  counter = 0;
  ddsigma = 0;
  if(store_param){
  	if(vec_counter == 0 || polyak_rate == -1)
  		sigma_vec[vec_counter] = sigma;
  	else
  		sigma_vec[vec_counter] = polyak_rate * sigma  + (1 - polyak_rate) * sigma_vec[vec_counter - 1];
  	vec_counter++;
  	}
}


std::vector< Eigen::VectorXd > GaussianMeasurementError::simulate(std::vector< Eigen::VectorXd > Y)
{
	std::vector< Eigen::VectorXd > residual( Y.size());
	for(int i = 0; i < Y.size(); i++)
		residual[i] =  sigma * (Rcpp::as< Eigen::VectorXd >(Rcpp::rnorm( Y[i].size()) ));

	return(residual);
}


Eigen::VectorXd  GaussianMeasurementError::simulate(const Eigen::VectorXd & Y)
{
	Eigen::VectorXd residual =  sigma * (Rcpp::as< Eigen::VectorXd >(Rcpp::rnorm( Y.size()) ));
	return(residual);
}


Eigen::VectorXd  GaussianMeasurementError::simulate_par(const Eigen::VectorXd & Y,std::mt19937 & random_engine)
{
  std::normal_distribution<double> normal;
  Eigen::VectorXd residual;
  residual.setZero(Y.size());
  for(int j =0; j < Y.size(); j++)
    residual[j] =  sigma*normal(random_engine);

  //Eigen::VectorXd residual =  sigma * (Rcpp::as< Eigen::VectorXd >(Rcpp::rnorm( Y.size()) ));
  return(residual);
}

Eigen::VectorXd  GaussianMeasurementError::simulate_par(const int i,std::mt19937 & random_engine, int nsim)
{
  std::normal_distribution<double> normal;
  Eigen::VectorXd residual;
  residual.setZero(nsim);
  for(int j =0; j < nsim; j++)
    residual[j] =  sigma*normal(random_engine);
  
  //Eigen::VectorXd residual =  sigma * (Rcpp::as< Eigen::VectorXd >(Rcpp::rnorm( Y.size()) ));
  return(residual);
}


void GaussianMeasurementError::clear_gradient()
{
	dsigma = 0;
}

Eigen::VectorXd GaussianMeasurementError::get_gradient()
{
	Eigen::VectorXd g(npars);
	g[0] = dsigma;
	return(g);
}
