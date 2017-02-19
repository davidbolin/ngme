#include "measError.h"
#include "error_check.h"




IGMeasurementError::IGMeasurementError() : NormalVarianceMixtureBaseError(){
  nu        = 1;
  beta  = nu + 1.;
  dnu       = 0;
  ddnu      = 0;
  dnu_old = 0;
  noise = "IG";

}



void IGMeasurementError::printIter()
{
	NormalVarianceMixtureBaseError::printIter();
	Rcpp::Rcout << "\n nu = " << nu;

}
void IGMeasurementError::setupStoreTracj(const int Niter) // setups to store the tracjetory
{

	NormalVarianceMixtureBaseError::setupStoreTracj(Niter);
	nu_vec.resize(Niter);
}





Rcpp::List IGMeasurementError::toList()
{
  Rcpp::List out = NormalVarianceMixtureBaseError::toList();
  out["nu"]          = nu;

  if(store_param){
  	out["nu_vec"]    = nu_vec;
  	out["nu"]          = nu_vec[nu_vec.size() - 1];
  }
  return(out);
}
void IGMeasurementError::initFromList(Rcpp::List const &init_list)
{

  NormalVarianceMixtureBaseError::initFromList(init_list);

  if(init_list.containsElementNamed("nu"))
    nu = Rcpp::as < double >( init_list["nu"]);
  else
    nu = 1.;

    beta = nu + 1.;
    EV  = 1.; // not true but the mode is one
    EiV = nu / beta;

   npars += 1;
  digamma_nu  =  R::digamma(nu);
  trigamma_nu =  R::trigamma(nu);

 int i = 0;

 if( init_list.containsElementNamed("Vs" )){
 	Rcpp::List Vs_list = init_list["Vs"];
 	Vs.resize(Vs_list.length());
    for( Rcpp::List::iterator it = Vs_list.begin(); it != Vs_list.end(); ++it ) {
      Vs[i++] = Rcpp::as < Eigen::VectorXd >( it[0]);
    }
 }else
 	  throw("in IGMeasurementError::initFromList Vs must be set! \n");


}

double IGMeasurementError::simulate_V()
{
	return rgig.sample(-nu, 0, 2 * beta );
}

double IGMeasurementError::sample_V(const double res2_j, const int n_s)
{
	if(common_V == 0)
		return rgig.sample(-(nu + .5), 0 , res2_j + 2 * beta);

	return rgig.sample(-  (nu + .5 * n_s), 0, res2_j + 2 * beta );
}





void IGMeasurementError::gradient(const int i,
                                 const Eigen::VectorXd& res,
                                 const double weight)
{
    NormalVarianceMixtureBaseError::gradient(i, res, weight);
    Eigen::VectorXd iV = Vs[i].cwiseInverse();
    if(common_V == 0){
    	double logV = Vs[i].array().log().sum();
      // beta = nu + 1
    	dnu  += weight * ( res.size() *  (nu / beta + log(beta)  - digamma_nu) - logV -  iV.array().sum());
    	ddnu += weight * (res.size() * ( 2/ beta - nu/(beta * beta) - trigamma_nu) );
    }else{

    	double logV = log(Vs[i][0]);
    	dnu  +=  weight * ( (nu / beta + log(beta)  - digamma_nu) - logV - iV[0] );
    	ddnu +=  weight * (2/ beta - nu/(beta * beta) - trigamma_nu);
    }
}

void IGMeasurementError::step_nu(const double stepsize, const double learning_rate,const int burnin)
{
	double nu_temp = -1;
  dnu /= ddnu;
  double stepsize_temp  =stepsize;
  dnu_old = dnu_old * learning_rate + dnu;
  while(nu_temp < 0)
  {
    nu_temp = nu - stepsize_temp * dnu;
    stepsize_temp *= 0.5;
    if(stepsize_temp <= 1e-16)
        throw("in IGMeasurementError:: can't make nu it positive \n");
  }
  nu = nu_temp;
  beta = nu + 1.;
  EV  = 1.;  // not true it is the mode that is 1.
  EiV = nu / beta ;

  ddnu = 0;
  digamma_nu  =  R::digamma(nu);
  trigamma_nu =  R::trigamma(nu);

}

void IGMeasurementError::step_theta(const double stepsize,
									const double learning_rate,
									const double polyak_rate,
									const int burnin)
{
  	NormalVarianceMixtureBaseError::step_theta(stepsize, learning_rate, polyak_rate, burnin);

  	step_nu(stepsize, learning_rate,burnin);
  	clear_gradient();

	if(store_param){

		if(vec_counter ==1 || polyak_rate == -1)
			nu_vec[vec_counter-1] = nu; // -1 since NormalVarianceMixtureBaseError increase vec_counter
		else
			nu_vec[vec_counter-1] =  polyak_rate * nu + (1- polyak_rate) * nu_vec[vec_counter-2];
	}
}

void IGMeasurementError::clear_gradient()
{
	NormalVarianceMixtureBaseError::clear_gradient();
	dnu    = 0;
}
Eigen::VectorXd IGMeasurementError::get_gradient()
{
	Eigen::VectorXd g = NormalVarianceMixtureBaseError::get_gradient();
	g[1] = dnu;
	return(g);
}
