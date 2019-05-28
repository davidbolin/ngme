#include "MixedEffect.h"
#include "error_check.h"
#include "GHmisc.h"
#include <chrono>

tdMixedEffect::tdMixedEffect(){
  counter = 0;
  noise = "tdist";
  npars = 0;
  store_param  = 0;
  accept_MALA = 0;
  count_MALA  = 0;
  weight_total = 0;
  sample_MALA = 0;
  calc_grad = 1;
  fixedV =  0;
}

double tdMixedEffect::get_p_GIG(){ return(-nu); }
void   tdMixedEffect::set_p_GIG(){}

double tdMixedEffect::get_a_GIG(){ return(0); }
void   tdMixedEffect::set_a_GIG(){

}
double tdMixedEffect::get_b_GIG(){ return(b_GIG); }
void   tdMixedEffect::set_b_GIG(){ b_GIG = 2. * (nu + 1);}

void tdMixedEffect::printIter()
{

GHMixedEffect::printIter();
	if(Br.size() > 0)
		Rcpp::Rcout << "nu     = " << nu << "\n";
}
void tdMixedEffect::setupStoreTracj(const int Niter)
{
  GHMixedEffect::setupStoreTracj(Niter);

	if(Br.size() > 0)
		nu_vec.resize(Niter);
}

void tdMixedEffect::get_param(std::vector<double> & param_in ){

  GHMixedEffect::get_param(param_in);
  if(Br.size() > 0 )
    param_in.push_back(nu);
}

void tdMixedEffect::get_param_names(Rcpp::StringVector & names){

  GHMixedEffect::get_param_names(names);
  if(Br.size() > 0 )
    names.push_back("nu");
}
Rcpp::List tdMixedEffect::toList()
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

void tdMixedEffect::initFromList(Rcpp::List const &init_list)
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
  	EV  = 1.;  // not true it is the mode that is 1.
  	EiV = nu / (get_b_GIG()) ;
    VV  = get_b_GIG() * get_b_GIG() /( (nu - 1) * (nu - 1) * (nu - 2) );


    set_a_GIG();
    set_b_GIG();
  }
}


double tdMixedEffect::logdensity(const Eigen::VectorXd &  U){

  return logGH(U,
                mu,
                -mu,
                invSigma,
                get_p_GIG(),
                get_a_GIG(),
                get_b_GIG());
}
void tdMixedEffect::updateFisher(const int i, 
                  Eigen::MatrixXd & Fisher, 
                  Eigen::VectorXd & grad){

  double grad_nu =  log(nu + 1) + nu/( nu + 1.) - R::digamma(nu) - log(V(i)) - 1. / V(i);
  grad(npars - 1) = 0.5 * grad_nu;
  Fisher.row(npars - 1) += - grad_nu * grad;
  Fisher.col(npars - 1) += - grad_nu * grad;
  grad(npars - 1) += 0.5 * grad_nu;

}

void tdMixedEffect::gradient(const int i,
                              const Eigen::VectorXd& res,
                              const double log_sigma2_noise,
                              const double weight,
                              const int use_EU
                              )
{
  GHMixedEffect::gradient(i, res, log_sigma2_noise, weight, use_EU);
  //Rcpp::Rcout <<" V(" << i <<") = " << V(i) << " U = " << U(i) << "\n";
  if(Br.size() > 0){
    double lik = log(nu + 1) + nu/( nu + 1.) - R::digamma(nu) - log(V(i)) - 1. / V(i);
    // dnu
    grad_nu += weight * lik;
  }
}

void tdMixedEffect::gradient2(const int i,
                                 const Eigen::VectorXd& res,
                                 const Eigen::VectorXd& iV,
                                 const Eigen::VectorXd& sigmas,  // =0
                                 const double log_sigma2_noise,  // = 0
                                 const double EiV, // = 0
                                 const double weight, //  = 1
                                 const int use_EU , // =1,
                                 const int nsigma           //  = 0
                                 )
{
  GHMixedEffect::gradient2(i, res, iV, sigmas, log_sigma2_noise, EiV, weight, use_EU, nsigma);
    if(Br.size() > 0){
      double lik = log(nu + 1.) + nu/( nu + 1.) - R::digamma(nu) - log(V(i)) - 1. / V(i);
      // dnu
      grad_nu += weight * lik;
    }
}


void tdMixedEffect::step_theta( const double stepsize,
                                const double learning_rate,
                                const double polyak_rate,
                                const int burnin)
{
  GHMixedEffect::step_theta(stepsize, learning_rate, polyak_rate, burnin);
  if(Br.size() > 0){
    step_nu(stepsize, learning_rate,burnin);
    set_a_GIG();
    set_p_GIG();
  }
  clear_gradient();
  store_param_function(polyak_rate);
}

void tdMixedEffect::store_param_function(const double polyak_rate)
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

void tdMixedEffect::step_nu(const double stepsize, const double learning_rate,const int burnin)
{
   grad_nu  /= -( (2+nu)/( pow(nu+1, 2)  ) - R::trigamma(nu) ) *  weight_total; //hessian

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




 if(nu_temp < 5e-06){
    nu_temp =5e-06;
    dnu_old = 0;
  }
  if(nu_temp  > 100){
      nu_temp =100;
      dnu_old = 0;
  }
  nu = nu_temp;
  set_b_GIG();
  set_a_GIG();

  EV  = 1.;  // not true it is the mode that is 1.
  EiV = nu / get_b_GIG() ;
  VV  = get_b_GIG() * get_b_GIG()/( (nu - 1) * (nu - 1) * (nu - 2) );

}


void tdMixedEffect::clear_gradient()
{

  GHMixedEffect::clear_gradient();
  if(Br.size() > 0)
    grad_nu = 0;
}


Eigen::VectorXd tdMixedEffect::get_gradient()
{
  Eigen::VectorXd g = GHMixedEffect::get_gradient();
  if(Br.size() > 0)
    g[npars - 1] = grad_nu;
  return(g);
}


Eigen::MatrixXd tdMixedEffect::d2Given( const int i,
                                        const Eigen::VectorXd& res,
                                        const double log_sigma2_noise,
                                        const double weight)
{
  Eigen::MatrixXd d2 = GHMixedEffect::d2Given(i, res, log_sigma2_noise, weight);
  return(d2);
}
Eigen::MatrixXd tdMixedEffect::d2Given2(const int i,
                                        const Eigen::VectorXd& res,
                                        const Eigen::VectorXd& iV,
                                        const double log_sigma2_noise,  // = 0
                                        const double EiV, // = 0
                                        const double weight //  = 1
                                       )
{
  Eigen::MatrixXd d2 = GHMixedEffect::d2Given2(i, res, iV, log_sigma2_noise, EiV, weight);
  return(d2);
}
