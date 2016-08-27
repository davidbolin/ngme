#include "measError.h"
#include "error_check.h"


void NIGMeasurementError::printIter() 
{
	Rcpp::Rcout << "(sigma , nu)= " << sigma << " ," << nu;

}
void NIGMeasurementError::setupStoreTracj(const int Niter) // setups to store the tracjetory
{
	sigma_vec.resize(Niter);
	nu_vec.resize(Niter);
	vec_counter = 0;
	store_param = 1;
}




NIGMeasurementError::NIGMeasurementError(){

  store_param = 0;
  counter   = 0;
  sigma     = 1;
  nu        = 1;
  dsigma    = 0;
  ddsigma   = 0;
  dnu       = 0;
  ddnu      = 0;
  common_V  = 0;
  npars     = 0;
  noise = "NIG";
  rgig.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());
}
Rcpp::List NIGMeasurementError::toList()
{
  Rcpp::List out;
  out["noise"]       = noise;
  out["sigma"]       = sigma;
  out["nu"]          = nu;
  out["Vs"]          = Vs;
  out["Cov_theta"]   = Cov_theta;
  out["noise"]       = noise;
  
  if(store_param){
  	out["sigma_vec"] = sigma_vec;
  	out["nu_vec"]    = nu_vec;
  }
  return(out);
}
void NIGMeasurementError::initFromList(Rcpp::List const &init_list)
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

  if(init_list.containsElementNamed("nu"))
    nu = Rcpp::as < double >( init_list["nu"]);
  else
    nu = 1.;

    EV  = 1.;
    EiV = 1. + 1./nu;

   npars += 1;
 int i = 0;

 if( init_list.containsElementNamed("Vs" )){
 	Rcpp::List Vs_list = init_list["Vs"];
 	Vs.resize(Vs_list.length());
    for( Rcpp::List::iterator it = Vs_list.begin(); it != Vs_list.end(); ++it ) {
      Vs[i++] = Rcpp::as < Eigen::VectorXd >( it[0]);
    }
 }else
 	  throw("in NigMeasurementError::initFromList Vs must be set! \n");


}

std::vector< Eigen::VectorXd > NIGMeasurementError::simulate(std::vector< Eigen::VectorXd > Y)
{

	std::vector< Eigen::VectorXd > residual( Y.size());
	for(int i = 0; i < Y.size(); i++){
		residual[i] =   sigma * (Rcpp::as< Eigen::VectorXd >(Rcpp::rnorm( Y[i].size()) ));
	  if(common_V == 0){
	    for(int ii = 0; ii < residual[i].size(); ii++)
	    {
	      double V = rgig.sample(-0.5, nu, nu);
	      Vs[i][ii] = V;
	      residual[i][ii] *=  sqrt(V);
	    }
	  } else {
	    double V = rgig.sample(-0.5, nu, nu);
	    for(int ii = 0; ii < residual[i].size(); ii++)
	    {
	      Vs[i][ii] = V;
	      residual[i][ii] *=  sqrt(V);
	    }
	  }


	}
	return(residual);
}

Eigen::VectorXd  NIGMeasurementError::simulate(const Eigen::VectorXd & Y)
{
	Eigen::VectorXd residual =  sigma * (Rcpp::as< Eigen::VectorXd >(Rcpp::rnorm( Y.size()) ));
	if(common_V == 0){
		for(int ii = 0; ii < residual.size(); ii++)
	    {
	      double V = rgig.sample(-0.5, nu, nu);
	      residual[ii] *=  sqrt(V);
	    }
	
	}else{
		  double V = rgig.sample(-0.5, nu, nu);
	      residual.array() *=  sqrt(V);
	}
	
	return(residual);
}


void NIGMeasurementError::sampleV(const int i, const Eigen::VectorXd& res, int n_s )
{
	if(n_s == -1)
		n_s = Vs[i].size();
	if(common_V == 0){
    for(int j = 0; j < n_s; j++)
	    Vs[i][j] = rgig.sample(-1., nu, pow(res[j]/sigma, 2) + nu);
	} else {
	  double tmp = res.array().square().sum()/pow(sigma, 2);
	  double cv = rgig.sample(-0.5*(n_s+1), nu, tmp + nu);
	  for(int j = 0; j < n_s; j++)
	    Vs[i][j] = cv;
  }
}

void NIGMeasurementError::gradient(const int i,
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

    dnu  += 0.5 * ( res.size() / nu + 2 * res.size() -   (Vs[i].array().sum() + iV.array().sum()) );
    ddnu += - 0.5*res.size()/( nu * nu);
}
void NIGMeasurementError::step_theta(double stepsize)
{
  step_sigma(stepsize);
  step_nu(stepsize);
  clear_gradient();
  counter = 0;
  
if(store_param){
  	sigma_vec[vec_counter] = sigma;
  	nu_vec[vec_counter++] = nu;
  }
}

void NIGMeasurementError::step_sigma(double stepsize)
{
double sigma_temp = -1;
  dsigma /= ddsigma;
  while(sigma_temp < 0)
  {
    sigma_temp = sigma - stepsize * dsigma;
    stepsize *= 0.5;
    if(stepsize <= 1e-16)
        throw("in NIGMeasurementError:: can't make sigma it positive \n");
  }
  sigma = sigma_temp;
  ddsigma = 0;
}

void NIGMeasurementError::step_nu(double stepsize)
{
double nu_temp = -1;
  dnu /= ddnu;
  while(nu_temp < 0)
  {
    nu_temp = nu - stepsize * dnu;
    stepsize *= 0.5;
    if(stepsize <= 1e-16)
        throw("in NIGMeasurementError:: can't make nu it positive \n");
  }
  nu = nu_temp;
  EV  = 1.;
  EiV = 1. + 1./nu;
  ddnu = 0;

}

void NIGMeasurementError::clear_gradient()
{
	dsigma = 0;
	dnu    = 0;
}

Eigen::VectorXd NIGMeasurementError::get_gradient()
{
	Eigen::VectorXd g(npars);
	g[0] = dsigma;
	g[1] = dnu;
	return(g);
}
