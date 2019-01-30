#include "measError.h"
#include "error_check.h"

NIGMeasurementError::NIGMeasurementError() : NormalVarianceMixtureBaseError(){
  nu        = 1;
  dnu       = 0;
  ddnu      = 0;
	dnu_old = 0;
	term1 = 0;
	term2 = 0;
  noise = "NIG";
  nsSigma = 0;

}



void NIGMeasurementError::printIter()
{
	NormalVarianceMixtureBaseError::printIter();
	Rcpp::Rcout << "\n nu = " << nu;

}
void NIGMeasurementError::setupStoreTracj(const int Niter) // setups to store the tracjetory
{

	NormalVarianceMixtureBaseError::setupStoreTracj(Niter);
	nu_vec.resize(Niter);
}





Rcpp::List NIGMeasurementError::toList()
{
  Rcpp::List out = NormalVarianceMixtureBaseError::toList();
  out["nu"]          = nu;
  out["sd"]          = sigma;

  if(store_param){
  	out["nu_vec"]    = nu_vec;
  	out["nu"]          = nu_vec(nu_vec.size() - 1);
  	}

  return(out);
}
void NIGMeasurementError::initFromList(Rcpp::List const &init_list)
{

  NormalVarianceMixtureBaseError::initFromList(init_list);
  if(init_list.containsElementNamed("nu"))
    nu = Rcpp::as < double >( init_list["nu"]);
  else
    nu = 1.;

  EV  = 1.;
  EiV = 1. + 1./nu;
  VV  = 1./nu;
	dnu_old = 0;
  npars += 1;
  int i = 0;

 if( init_list.containsElementNamed("Vs" )){
 	Rcpp::List Vs_list = init_list["Vs"];
 	Vs.resize(Vs_list.length());
    for( Rcpp::List::iterator it = Vs_list.begin(); it != Vs_list.end(); ++it ) {
      Vs[i++] = Rcpp::as < Eigen::VectorXd >( it[0]);
    }
 }else{
    Rcpp::Rcout << "in NigMeasurementError::initFromList Vs must be set! \n";
 	  throw("in NigMeasurementError::initFromList Vs must be set! \n");
  }


}

double NIGMeasurementError::simulate_V()
{
	return rgig.sample(-0.5, nu, nu);
}

double NIGMeasurementError::sample_V(const double res2_j, const int n_s, const double mu2)
{
	if(common_V == 0)
		return rgig.sample(-1., nu + mu2, res2_j + nu);

	return rgig.sample(-0.5 * (n_s + 1), nu + mu2, res2_j + nu);
}



void NIGMeasurementError::gradient(const int i,
                                 const Eigen::VectorXd& res,
                                 const double weight)
{
    NormalVarianceMixtureBaseError::gradient(i, res, weight);
    Eigen::VectorXd iV = Vs[i].cwiseInverse();
    if(common_V == 0){
    	dnu   += weight * (0.5 * ( res.size() / nu + 2 * res.size() -   (Vs[i].array().sum() + iV.array().sum()) ));
    	ddnu  += weight * (- 0.5*res.size()/( nu * nu));
    	term1 += weight * res.size();
    	term2 += weight * (Vs[i].array().sum() + iV.array().sum()- 2 * res.size());

    }else{

    	dnu   += weight * (0.5 * ( 1. / nu + 2  -   (Vs[i][0] + iV[0]) ));
    	ddnu  += weight * (- 0.5*  1. / ( nu * nu));
    	term1 += weight * 1;
    	term2 += weight *( Vs[i][0] + iV[0] -2);
    }

}

void NIGMeasurementError::step_nu(const double stepsize, const double learning_rate,const int burnin)
{
double nu_temp = -1;
  dnu /= ddnu;
  dnu_old = learning_rate  * dnu_old + dnu;
  double step_size = stepsize;
  while(nu_temp < 0)
  {
    nu_temp = nu - step_size * dnu_old;
    step_size *= 0.5;
    if(step_size <= 1e-16)
        throw("in NIGMeasurementError:: can't make nu it positive \n");
  }
  if(burnin == 1){
    nu_temp = term1/term2;
    if(nu_temp < 0){
      nu_temp = 0.1;
    }
  }

  nu = nu_temp;
  EV  = 1.;
  EiV = 1. + 1./nu;
  ddnu = 0;
  VV  = 1./nu;

}

void NIGMeasurementError::step_theta(const double stepsize,
                                     const double learning_rate,
                                     const double polyak_rate,
                                     const int burnin)
{
  	NormalVarianceMixtureBaseError::step_theta(stepsize, learning_rate, polyak_rate, burnin);

  	step_nu(stepsize, learning_rate,burnin);
  	clear_gradient();

	if(store_param){
		if(vec_counter == 1 || polyak_rate == -1)
  			nu_vec(vec_counter-1) = nu; // -1 since NormalVarianceMixtureBaseError increase vec_counter
  		else
  			nu_vec(vec_counter-1) = polyak_rate * nu + (1 - polyak_rate) * nu_vec(vec_counter - 2);

  	}

}

void NIGMeasurementError::clear_gradient()
{
	NormalVarianceMixtureBaseError::clear_gradient();
	dnu    = 0;
	term1= 0;
	term2 = 0;
}
Eigen::VectorXd NIGMeasurementError::get_gradient()
{
	Eigen::VectorXd g = NormalVarianceMixtureBaseError::get_gradient();
	g(npars-1) = dnu;
	return(g);
}


Eigen::MatrixXd NIGMeasurementError::d2Given(const int i,
                    const Eigen::VectorXd& res,
                    const double weight)
{

  Eigen::MatrixXd d2 =   NormalVarianceMixtureBaseError::d2Given(i, res, weight);
  if(common_V == 0)
    d2(npars-2, npars-2) = weight *  (0.5 *  res.size() / ( nu * nu));
  else
    d2(npars-2, npars-2) = weight *  (0.5 *   1. / ( nu * nu));
  return(d2);
}
