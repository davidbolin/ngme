#include "measError.h"
#include "error_check.h"


void NormalVarianceMixtureBaseError::printIter() 
{
	Rcpp::Rcout << "sigma = " << sigma;

}
void NormalVarianceMixtureBaseError::setupStoreTracj(const int Niter) // setups to store the tracjetory
{
	sigma_vec.resize(Niter);
	vec_counter = 0;
	store_param = 1;
}




NormalVarianceMixtureBaseError::NormalVarianceMixtureBaseError(){

  store_param = 0;
  counter   = 0;
  sigma     = 1;
  dsigma    = 0;
  ddsigma   = 0;
  common_V  = 0;
  npars     = 0;
  rgig.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());
}
Rcpp::List NormalVarianceMixtureBaseError::toList()
{
  Rcpp::List out;
  out["noise"]       = noise;
  out["sigma"]       = sigma;
  out["nu"]          = nu;
  out["Cov_theta"]   = Cov_theta;
  out["noise"]       = noise;
  
  if(store_param)
  	out["sigma_vec"] = sigma_vec;

  return(out);
}
void NormalVarianceMixtureBaseError::initFromList(Rcpp::List const &init_list)
{
  if(init_list.containsElementNamed("sigma"))
    sigma = Rcpp::as < double >( init_list["sigma"]);
  else
  	sigma  = 1.;

   npars += 1;
  if(init_list.containsElementNamed("common_V"))
    common_V = Rcpp::as < int >( init_list["common_V"]);
  else
    common_V  = 0;


 int i = 0;

 if( init_list.containsElementNamed("Vs" )){
 	Rcpp::List Vs_list = init_list["Vs"];
 	Vs.resize(Vs_list.length());
    for( Rcpp::List::iterator it = Vs_list.begin(); it != Vs_list.end(); ++it ) {
      Vs[i++] = Rcpp::as < Eigen::VectorXd >( it[0]);
    }
 }else
 	  throw("in NormalVarianceMixtureBaseError::initFromList Vs must be set! \n");


}

double NormalVarianceMixtureBaseError::simulate_V()
{
	return -1;
}
Eigen::VectorXd  NormalVarianceMixtureBaseError::simulate(const Eigen::VectorXd & Y)
{
	Eigen::VectorXd residual =  sigma * (Rcpp::as< Eigen::VectorXd >(Rcpp::rnorm( Y.size()) ));
	if(common_V == 0){
		for(int ii = 0; ii < residual.size(); ii++)
	    {
	      double V = simulate_V();
	      residual[ii] *=  sqrt(V);
	    }
	
	}else{
		  double V = simulate_V();
	      residual.array() *=  sqrt(V);
	}
	
	return(residual);
}

std::vector< Eigen::VectorXd > NormalVarianceMixtureBaseError::simulate(std::vector< Eigen::VectorXd > Y)
{

	std::vector< Eigen::VectorXd > residual( Y.size());
	for(int i = 0; i < Y.size(); i++){
		residual[i] =   sigma * (Rcpp::as< Eigen::VectorXd >(Rcpp::rnorm( Y[i].size()) ));
	  if(common_V == 0){
	    for(int ii = 0; ii < residual[i].size(); ii++)
	    {
	      double V = simulate_V();
	      Vs[i][ii] = V;
	      residual[i][ii] *=  sqrt(V);
	    }
	  } else {
	    double V = simulate_V();
	    for(int ii = 0; ii < residual[i].size(); ii++)
	    {
	      Vs[i][ii] = V;
	      residual[i][ii] *=  sqrt(V);
	    }
	  }


	}
	return(residual);
}

void NormalVarianceMixtureBaseError::sampleV(const int i, const Eigen::VectorXd& res, int n_s )
{
	if(n_s == -1)
		n_s = Vs[i].size();
	if(common_V == 0){
    	for(int j = 0; j < n_s; j++){
    		double res2 = pow(res[j]/sigma, 2);
	    	Vs[i][j] = sample_V(res2, -1);
		}
	    
	} else {
	  double tmp = res.array().square().sum()/pow(sigma, 2);
	  double cv = sample_V(tmp, n_s);
	  for(int j = 0; j < n_s; j++)
	    Vs[i][j] = cv;
  }
}

void NormalVarianceMixtureBaseError::gradient(const int i,
                                 const Eigen::VectorXd& res)
{
    counter++;
    Eigen::VectorXd res_ = res;
    Eigen::VectorXd iV = Vs[i].cwiseInverse();
    //res_.array() *= iV.array();
    dsigma += - res.size()/sigma + (res_.array().square()*iV.array()).sum() / pow(sigma, 3);
    // Expected fisher infromation
    // res.size()/pow(sigma, 2) - 3 * E[res.array().square().sum()] /pow(sigma, 4);
    ddsigma += - 2 * res.size()/pow(sigma, 2);
}

void NormalVarianceMixtureBaseError::step_theta(double stepsize)
{
  step_sigma(stepsize);
  NormalVarianceMixtureBaseError::clear_gradient();
  
  counter = 0;
  
  sigma_vec[vec_counter++] = sigma;
  
}

void NormalVarianceMixtureBaseError::step_sigma(double stepsize)
{
double sigma_temp = -1;
  dsigma /= ddsigma;
  while(sigma_temp < 0)
  {
    sigma_temp = sigma - stepsize * dsigma;
    stepsize *= 0.5;
    if(stepsize <= 1e-16)
        throw("in NormalVarianceMixtureBaseError:: can't make sigma it positive \n");
  }
  sigma = sigma_temp;
  ddsigma = 0;
}



void NormalVarianceMixtureBaseError::clear_gradient()
{
	dsigma = 0;
}

Eigen::VectorXd NormalVarianceMixtureBaseError::get_gradient()
{
	Eigen::VectorXd g(npars);
	g[0] = dsigma;
	return(g);
}
