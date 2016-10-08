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


GaussianMeasurementError::GaussianMeasurementError(){
  counter = 0;
  sigma   = 0;
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
  out["noise"]  = noise;
  out["Cov_theta"]   = Cov_theta;
  if(store_param)
  	out["sigma_vec"] = sigma_vec;

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
                                 const Eigen::VectorXd& res)
{
    counter++;
    dsigma += - res.size()/sigma + res.array().square().sum() / pow(sigma, 3);
    // Expected fisher infromation
    // res.size()/pow(sigma, 2) - 3 * E[res.array().square().sum()] /pow(sigma, 4);
    ddsigma += - 2 * res.size()/pow(sigma, 2);
}
void GaussianMeasurementError::step_theta(const double stepsize, const double learning_rate)
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
  if(store_param)
  	sigma_vec[vec_counter++] = sigma;
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
