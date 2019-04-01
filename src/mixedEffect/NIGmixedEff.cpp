#include "MixedEffect.h"
#include "error_check.h"
#include "GHmisc.h"
#include <chrono>

NIGMixedEffect::NIGMixedEffect(){
  counter = 0;
  noise = "NIG";
  npars = 0;
  store_param  = 0;
  accept_MALA = 0;
  count_MALA  = 0;
  weight_total = 0;
  sample_MALA = 0;
  calc_grad  = 1;
  fixedV  = 0;
}


double NIGMixedEffect::get_p_GIG(){ return(-0.5); }
void   NIGMixedEffect::set_p_GIG(){}

double NIGMixedEffect::get_a_GIG(){ return(a_GIG); }
void   NIGMixedEffect::set_a_GIG(){
  a_GIG = nu;
}
double NIGMixedEffect::get_b_GIG(){ return(b_GIG); }
void   NIGMixedEffect::set_b_GIG(){ b_GIG = nu; }



void NIGMixedEffect::printIter()
{

  GHMixedEffect::printIter();
	if(Br.size() > 0)
		Rcpp::Rcout << "nu     = " << nu << "\n";
}
void NIGMixedEffect::setupStoreTracj(const int Niter)
{
  GHMixedEffect::setupStoreTracj(Niter);

	if(Br.size() > 0)
		nu_vec.resize(Niter);
}

void NIGMixedEffect::get_param(std::vector<double> & param_in ){

  GHMixedEffect::get_param(param_in);
  if(Br.size() > 0 )
    param_in.push_back(nu);
}

void NIGMixedEffect::get_param_names(Rcpp::StringVector & names){

  GHMixedEffect::get_param_names(names);
  if(Br.size() > 0 )
    names.push_back("nu");
}

Rcpp::List NIGMixedEffect::toList()
{
  Rcpp::List out = GHMixedEffect::toList();
  
  out["nu"]     = nu;
  if(store_param){
	 if(Br.size() > 0){
		out["nu_vec"]      = nu_vec;
    if(betar_vec.rows() > 1)
		  out["nu"]          = nu_vec[nu_vec.size() - 1];
	 }
  }
  return(out);
}

void NIGMixedEffect::initFromList(Rcpp::List const &init_list)
{

  GHMixedEffect::initFromList(init_list);

  if(init_list.containsElementNamed("B_random"))
  {
    if( init_list.containsElementNamed("nu"))
      nu = Rcpp::as< double > (init_list["nu"]) ;
    else
      nu = 1.;

    set_a_GIG();
    set_b_GIG();

    if( init_list.containsElementNamed("V" ) == FALSE)
       simulate();
    
    
    
    npars += 1;
	  dnu_old = 0;
    grad_nu = 0.;
    EV  = 1.;
    EiV = 1. + 1./nu;
    VV  = 1./nu;


    set_a_GIG();
    set_b_GIG();
  }
}



double NIGMixedEffect::logdensity(const Eigen::VectorXd &  U){

  return logNIG(U,
                mu,
                -mu,
                invSigma,
                nu);
}


void NIGMixedEffect::gradient(const int i,
                              const Eigen::VectorXd& res,
                              const double log_sigma2_noise,
                              const double weight,
                              const int use_EU
                              )
{
  GHMixedEffect::gradient(i, res,log_sigma2_noise, weight, use_EU);
  if(Br.size() > 0){
    // dnu
    grad_nu += weight * 0.5 * (1. / nu - V(i) - 1. / V(i) + 2. );

    term1 += weight * (V(i) + 1. / V(i) - 2.);
    term2 += weight * 1.;
  }
}




void NIGMixedEffect::gradient2(const int i,
                                 const Eigen::VectorXd& res,
                                 const Eigen::VectorXd& iV,
                                 const Eigen::VectorXd& sigmas,  // =0
                                 const double log_sigma2_noise,  // = 0
                                 const double EiV, // = 0
                                 const double weight, //  = 1
                                 const int use_EU , // =1,
                                 const int nsigma            //  = 0
                                 )
{
  GHMixedEffect::gradient2(i, res, iV, sigmas, log_sigma2_noise, EiV, weight, use_EU, nsigma);
    if(Br.size() > 0){
      // dnu
      grad_nu += weight * 0.5 * (1. / nu - V(i) - 1. / V(i) + 2. );

      term1   += weight * (V(i) + 1. / V(i) - 2.);
      term2   += weight * 1.;
    }
}


void NIGMixedEffect::step_theta(const double stepsize,
								const double learning_rate,
								const double polyak_rate,
								const int burnin)
{
  GHMixedEffect::step_theta(stepsize, learning_rate, polyak_rate, burnin);
  if(Br.size() > 0){
    step_nu(stepsize, learning_rate,burnin);
    set_a_GIG();
  }
  clear_gradient();
  store_param_function(polyak_rate);
}

void NIGMixedEffect::store_param_function(const double polyak_rate)
{
  if(Br.size() > 0)
  {
    if(vec_counter == 0 || polyak_rate == -1){
      nu_vec[vec_counter] = nu;
    }else{
      nu_vec[vec_counter]     = polyak_rate * nu + (1 - polyak_rate) * nu_vec[vec_counter - 1];
    }
  }
  GHMixedEffect::store_param_function(polyak_rate);
}

void NIGMixedEffect::step_nu(const double stepsize, const double learning_rate,const int burnin)
{
   grad_nu  *=  (nu * nu) / (2. * weight_total); //hessian

  dnu_old = learning_rate * dnu_old + grad_nu;
  double nu_temp = -1;
  double step_size_temp = stepsize;

  while(nu_temp < 0){
  	nu_temp = nu + stepsize  * grad_nu;
    step_size_temp *= 0.5;
    if(step_size_temp <= 1e-16){
        Rcpp::Rcout << "nu = \n" << nu << "\n";
        Rcpp::Rcout << "grad_nu = " << grad_nu <<"\n";
        throw("in NIGmidexeffect nu is zero \n");
    }
  }

  if(burnin == 1){
    nu_temp = term1/term2;
    if(nu_temp < 0){
      nu_temp = 0.1;
    }
  }

  if(nu_temp  > 600){
  		nu_temp =600;
      dnu_old = 0;
  } else if(nu_temp < 5e-06){
    nu_temp =5e-06;
    dnu_old = 0;
  }
	nu = nu_temp;
  //nu = 400;
  EiV = 1. + 1./nu;
  VV = 1./nu;
  set_b_GIG();
  set_a_GIG();
}



void NIGMixedEffect::clear_gradient()
{

  GHMixedEffect::clear_gradient();
  term1 = 0.;
  term2 = 0.;
	if(Br.size() > 0)
		grad_nu = 0;
}


Eigen::VectorXd NIGMixedEffect::get_gradient()
{
  Eigen::VectorXd g = GHMixedEffect::get_gradient();
	if(Br.size() > 0)
    g[npars - 1] = grad_nu;
	return(g);
}


Eigen::MatrixXd NIGMixedEffect::d2Given( const int i,
                                        const Eigen::VectorXd& res,
                                        const double log_sigma2_noise,
                                        const double weight)
{
  Eigen::MatrixXd d2 = GHMixedEffect::d2Given(i, res, log_sigma2_noise, weight);
  d2(npars - 1 , npars - 1 ) =  weight * 0.5 / pow(nu,2);
  return(d2);
}
Eigen::MatrixXd NIGMixedEffect::d2Given2(const int i,
                                        const Eigen::VectorXd& res,
                                        const Eigen::VectorXd& iV,
                                        const double log_sigma2_noise,  // = 0
                                        const double EiV, // = 0
                                        const double weight //  = 1
                                       )
{
  Eigen::MatrixXd d2 = GHMixedEffect::d2Given2(i, res, iV, log_sigma2_noise, EiV, weight);
  d2(npars - 1 , npars - 1 ) =  weight * 0.5 / pow(nu,2);
  return(d2);
}
