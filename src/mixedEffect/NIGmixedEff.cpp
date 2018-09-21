#include "MixedEffect.h"
#include "error_check.h"
#include "NGIG.h"
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
}


double NIGMixedEffect::get_p_GIG(){ return(-0.5); }
void   NIGMixedEffect::set_p_GIG(){}

double NIGMixedEffect::get_a_GIG(){ return(a_GIG); }
void   NIGMixedEffect::set_a_GIG(){
  a_GIG = mu.transpose() * (invSigma *  mu);
  a_GIG += nu;
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

void NIGMixedEffect::sampleU(const int i,
                                const Eigen::VectorXd& res,
                                const double log_sigma2_noise)
{

  	if(Br.size() == 0)
  		return;
  	if(sample_MALA)
  	{
  	  sampleU_MALA(i,
                  res,
                  log_sigma2_noise);
  	}else{
  	  sampleU_Gibbs(i,
                   res,
                   log_sigma2_noise);
  	}
    sampleV(i);

    Eigen::VectorXd b   = exp( - log_sigma2_noise) * (Br[i].transpose() * res);
    Eigen::MatrixXd Q   = exp( - log_sigma2_noise)  * Br[i].transpose() * Br[i];
    Eigen::MatrixXd Qp  = invSigma;
    Qp.array() /= V(i);
    b += Qp * (- mu + V(i) * mu);
    Q += Qp;
    EU =  Q.ldlt().solve(b);

}

void NIGMixedEffect::sampleU2(const int i,
                                const Eigen::VectorXd& res,
                                const Eigen::VectorXd& iV,
                                const double log_sigma2_noise //= 0
                                )
{
    if(Br.size() == 0)
      return;

    if(sample_MALA)
    {
      sampleU2_MALA(i,
                    res,
                    iV,
                    log_sigma2_noise);
    }else{
      sampleU2_Gibbs(i,
                    res,
                    iV,
                    log_sigma2_noise);
    }
    sampleV(i);

    Eigen::VectorXd b   = exp( - log_sigma2_noise) * (Br[i].transpose() * iV.cwiseProduct(res) );
    Eigen::MatrixXd Q   = exp( - log_sigma2_noise)  * Br[i].transpose() * iV.asDiagonal() * Br[i];
    Eigen::MatrixXd Qp  = invSigma / V(i);
    b += Qp * (- mu + V(i) * mu);
    Q += Qp;
    EU =  Q.ldlt().solve(b);



}

double NIGMixedEffect::logdensity(const Eigen::VectorXd &  U){

  return logNIG(U,
                mu,
                -mu,
                invSigma,
                nu);
}


void NIGMixedEffect::sampleU_par(const int i,
                             const Eigen::VectorXd& res,
                             const double log_sigma2_noise,
                             std::mt19937 & random_engine)
{

  if(Br.size() == 0)
    return;

  Eigen::VectorXd b   = exp( - log_sigma2_noise) * (Br[i].transpose() * res);
  Eigen::MatrixXd Q   = exp( - log_sigma2_noise)  * Br[i].transpose() * Br[i];
  Eigen::MatrixXd Qp  = invSigma / V(i);
  b += Qp * (- mu + V(i) * mu);
  Q += Qp;
  U.col(i) = sample_Nc_par(b, Q,random_engine);
  sampleV(i);
  if(Sigma_epsilon){
    std::normal_distribution<double> normal;
    Eigen::VectorXd Z;
    Z.setZero(U.col(i).size());
    for(int j =0; j < U.col(i).size(); j++)
      Z[j] =  normal(random_engine);

    Z.array() *= beta_random.array().abs() * 1e-4 + 1e-14;
    U.col(i) += Z;
  }
}

void NIGMixedEffect::sampleU_Gibbs(const int i,
                                const Eigen::VectorXd& res,
                                const double log_sigma2_noise)
{
    Eigen::VectorXd b   = exp( - log_sigma2_noise) * (Br[i].transpose() * res);
    Eigen::MatrixXd Q   = exp( - log_sigma2_noise)  * Br[i].transpose() * Br[i];
    Eigen::MatrixXd Qp  = invSigma;
    Qp.array() /= V(i);
    b += Qp * (- mu + V(i) * mu);
    Q += Qp;
    U.col(i) = sample_Nc(b, Q);
    if(Sigma_epsilon){
      Eigen::VectorXd Z = Rcpp::as< Eigen::VectorXd >(Rcpp::rnorm( U.col(i).size()) );
      Z.array() *= beta_random.array().abs() * 1e-8 + 1e-14;
      U.col(i) += Z;
    }
}

void NIGMixedEffect::sampleU2_Gibbs(const int i,
                                const Eigen::VectorXd& res,
                                const Eigen::VectorXd& iV,
                                const double log_sigma2_noise //= 0
                                )
{
    Eigen::VectorXd b   = exp( - log_sigma2_noise) * (Br[i].transpose() * iV.cwiseProduct(res) );
    Eigen::MatrixXd Q   = exp( - log_sigma2_noise)  * Br[i].transpose() * iV.asDiagonal() * Br[i];
    Eigen::MatrixXd Qp  = invSigma / V(i);
    b += Qp * (- mu + V(i) * mu);
    Q += Qp;
    U.col(i) = sample_Nc(b, Q);
    if(Sigma_epsilon){
      Eigen::VectorXd Z = Rcpp::as< Eigen::VectorXd >(Rcpp::rnorm( U.col(i).size()) );
      Z.array() *= beta_random.array().abs() * 1e-8 + 1e-14;
      U.col(i) += Z;
    }
}

void NIGMixedEffect::sampleU2_par(const int i,
                              const Eigen::VectorXd& res,
                              const Eigen::VectorXd& iV,
                              std::mt19937 & random_engine,
                              const double log_sigma2_noise //= 0
)
{
  if(Br.size() == 0)
    return;

  Eigen::VectorXd b   = exp( - log_sigma2_noise) * (Br[i].transpose() * iV.cwiseProduct(res) );
  Eigen::MatrixXd Q   = exp( - log_sigma2_noise)  * Br[i].transpose() * iV.asDiagonal() * Br[i];
  Eigen::MatrixXd Qp  = invSigma / V(i);
  b += Qp * (- mu + V(i) * mu);
  Q += Qp;
  U.col(i) = sample_Nc_par(b, Q,random_engine);
  sampleV(i);
  if(Sigma_epsilon){
    std::normal_distribution<double> normal;
    Eigen::VectorXd Z;
    Z.setZero(U.col(i).size());
    for(int j =0; j < U.col(i).size(); j++)
      Z[j] =  normal(random_engine);

    Z.array() *= beta_random.array().abs() * 1e-8 + 1e-14;
    U.col(i) += Z;
  }
}

void NIGMixedEffect::sampleU_MALA(const int i,
                                const Eigen::VectorXd& res,
                                const double log_sigma2_noise)
{
    Eigen::MatrixXd Q_noise = MatrixXd::Identity(res.size(),res.size());
    Q_noise *= exp( - log_sigma2_noise);
    Eigen::VectorXd b_prior   =  (Br[i].transpose() * (Q_noise * res));

    sampleU_MALA_(i, res, b_prior, Q_noise);
}

void NIGMixedEffect::sampleU2_MALA(const int i,
                                const Eigen::VectorXd& res,
                                const Eigen::VectorXd& iV,
                                const double log_sigma2_noise //= 0
                                )
{

    Eigen::VectorXd b_prior   = exp( - log_sigma2_noise) * (Br[i].transpose() * iV.cwiseProduct(res) );
    Eigen::MatrixXd Q_noise   = exp( - log_sigma2_noise)  * iV.asDiagonal();
    sampleU_MALA_(i, res, b_prior, Q_noise);
}






void NIGMixedEffect::gradient(const int i,
                              const Eigen::VectorXd& res,
                              const double log_sigma2_noise,
                              const double weight,
                              const int use_EU // =1
                              )
{
    Eigen::VectorXd res_  = res;
    if(Br.size() > 0){

      Eigen::VectorXd U_ = U.col(i) - (-1 + V(i)) * mu;
      gradient_sigma(i, U_, weight);

      if(use_EU)
      {
        res_ -= Br[i] * EU;
        U_ = EU - (-1 + V(i)) * mu;
      }else{
        res_ -= Br[i] * U.col(i);
      }


      grad_beta_r  +=  weight * exp( - log_sigma2_noise) * (Br[i].transpose() * res_);
      grad_beta_r2 +=  weight * (invSigma * U_)/V(i);
      H_beta_random += weight * exp( - log_sigma2_noise) * (Br[i].transpose() * Br[i]);


      gradMu   += weight * ((-1 + V(i) )/V(i) ) * (invSigma * U_);

      gradMu_2 += weight * (-1 + V(i) ) * exp( - log_sigma2_noise) * (Br[i].transpose() * res_);

      term1_mu += weight * ((-1 + V(i) )/V(i) )*(-1 + V(i) );
      term2_mu += weight * ((-1 + V(i))/V(i) )*(invSigma*U_) + weight*(-1 + V(i))*exp(-log_sigma2_noise)*(Br[i].transpose()*res_);

      // dnu
      grad_nu += weight * 0.5 * (1. / nu - V(i) - 1. / V(i) + 2. );
      term1 += weight * (V(i) + 1. / V(i) - 2.);
      term2 += weight * 1.;
    }
    if(Bf.size() > 0){
      grad_beta_f   +=  weight * exp( - log_sigma2_noise) * (Bf[i].transpose() * res_);
      H_beta_fixed  +=  weight * exp( - log_sigma2_noise) * (Bf[i].transpose() * Bf[i]);
    }
    counter++;
  weight_total += weight;
}


Eigen::MatrixXd NIGMixedEffect::d2Given(const int i,
                                           const Eigen::VectorXd& res,
                                           const double log_sigma2_noise,
                                           const double weight)
{

  Eigen::VectorXd res_  = res;

  int n_s = 0;
  int n_f = 0;
  int n_r = 0;
  if(Br.size()>0){
     n_r = Br[i].cols();
     n_s = n_r * (n_r +1) /2;
   }
   if(Bf.size() > 0 )
      n_f = Bf[i].cols();

  Eigen::MatrixXd d2            = Eigen::MatrixXd::Zero(n_s + n_f + 2 * n_r + 2,n_s + n_f+ 2 * n_r + 2);
  if(Br.size()>0){
    double B_mu =  (-1 + V(i));
    Eigen::VectorXd U_ = U.col(i) - B_mu * mu;
    //beta_r
    d2.block(  n_f      , n_f       , n_r, n_r)  =  weight * exp( - log_sigma2_noise) * (Br[i].transpose() * Br[i]);

    //mu
    d2.block( n_f  + n_r, n_f   + n_r, n_r, n_r)  =  weight * ((B_mu * B_mu ) / V(i) ) * invSigma;

    res_ -= Br[i] * U.col(i);
    //Sigma
    d2.block(2 * n_r +n_f , 2 * n_r +n_f, n_s, n_s)  -=  0.5* weight * Dd.transpose() * iSkroniS * Dd;
    Eigen::MatrixXd UUT = U_ * U_.transpose();
    UUT.array() /= V(i);
    d2.block(2 * n_r +n_f , 2 * n_r + n_f, n_s, n_s)  += weight *  Dd.transpose() * kroneckerProduct(invSigma, invSigma * UUT * invSigma ) * Dd;

    // dmu dSigma

    d2.block( n_r +n_f , 2 * n_r +n_f, n_r, n_s)  =  kroneckerProduct(B_mu * invSigma, (invSigma * U_).transpose()) * Dd;
    d2.block( n_r +n_f , 2 * n_r +n_f, n_r, n_s)  *= 0.5* weight * (1 / V(i));
    d2.block(  2 * n_r +n_f, n_r +n_f , n_s, n_r)  =  d2.block( n_r +n_f , 2 * n_r +n_f, n_r, n_s).transpose();
  }
  if(Br.size() * Bf.size()>0){
    d2.block(  0      , n_f     , n_f, n_r)  =  weight * exp( - log_sigma2_noise) * (Bf[i].transpose() * Br[i]);
    d2.block(n_f      , 0       , n_r, n_f)  =  weight * exp( - log_sigma2_noise) * (Br[i].transpose() * Bf[i]);
  }
  if(Bf.size() > 0 )
    d2.block(0      , 0     , n_f, n_f)  =  weight * exp( - log_sigma2_noise) * (Bf[i].transpose() * Bf[i]);

  if(Br.size()>0){
    d2.block(n_f              , n_s + 2 * n_r +n_f + 1, n_r , 1 )     =  2 * weight * exp( - 1.5 * log_sigma2_noise)  * (Br[i].transpose() * res_);
    d2.block(n_s + 2 * n_r + n_f  + 1, n_f              , 1   , n_r ) = d2.block(n_f , n_s + 2 * n_r +n_f + 1, n_r , 1 ) .transpose();
  }
 if(Bf.size() > 0){
  // dbeta_r dsigma
    d2.block(0            , n_s + 2 * n_r + n_f + 1 , n_f , 1 )   =  2 * weight * exp( - 1.5 * log_sigma2_noise)  * (Bf[i].transpose() * res_);
    d2.block(n_s + 2 * n_r + n_f + 1, 0             , 1   , n_f ) = d2.block(0            , n_s + 2 * n_r + n_f + 1, n_f , 1 ).transpose();
  }
  d2(n_s +  2 * n_r +n_f , n_s +  2 * n_r +n_f ) =  weight * 0.5 / pow(nu,2);
  d2(n_s +  2 * n_r +n_f + 1, n_s +  2 * n_r +n_f + 1) =  3  * weight * exp( - 2   * log_sigma2_noise)  * res_.array().square().sum();
  d2(n_s +  2 * n_r +n_f + 1, n_s +  2 * n_r +n_f + 1) +=  -1 * weight * res_.size()  * exp( - log_sigma2_noise);

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

  Eigen::VectorXd res_  = res;

  int n_s = 0;
  int n_f = 0;
  int n_r = 0;
  if(Br.size()>0){
     n_r = Br[i].cols();
     n_s = n_r * (n_r +1) /2;
   }
   if(Bf.size() > 0 )
      n_f = Bf[i].cols();

  Eigen::MatrixXd d2            = Eigen::MatrixXd::Zero(n_s + n_f + 2 * n_r + 2,n_s + n_f+ 2 * n_r + 2);
  if(Br.size()>0){
    double B_mu =  (-1 + V(i));
    Eigen::VectorXd U_ = U.col(i) - B_mu * mu;
    //beta_r
    d2.block(  n_f      , n_f       , n_r, n_r)  =  weight * exp( - log_sigma2_noise) * (Br[i].transpose() * iV.asDiagonal() * Br[i]);

    //mu
    d2.block( n_f  + n_r, n_f   + n_r, n_r, n_r)  =  weight * ((B_mu * B_mu ) / V(i) ) * invSigma;
      res_ -= Br[i] * U.col(i);
    //Sigma

    d2.block(2 * n_r +n_f , 2 * n_r +n_f, n_s, n_s)  -=  0.5* weight * Dd.transpose() * iSkroniS * Dd;
    Eigen::MatrixXd UUT = U_ * U_.transpose();
    UUT.array() /= V(i);
    d2.block(2 * n_r +n_f , 2 * n_r + n_f, n_s, n_s)  += weight *  Dd.transpose() * kroneckerProduct(invSigma, invSigma * UUT * invSigma ) * Dd;
    // dmu dSigma

    d2.block( n_r +n_f , 2 * n_r +n_f, n_r, n_s)   =  kroneckerProduct(B_mu * invSigma, (invSigma * U_).transpose()) * Dd;
    d2.block( n_r +n_f , 2 * n_r +n_f, n_r, n_s)   *=  0.5* weight * (1 / V(i));
    d2.block(  2 * n_r +n_f, n_r +n_f , n_s, n_r)  =  d2.block( n_r +n_f , 2 * n_r +n_f, n_r, n_s).transpose();
  }
  if(Br.size() * Bf.size()>0){
    d2.block(  0      , n_f     , n_f, n_r)  =  weight * exp( - log_sigma2_noise) * (Bf[i].transpose() * iV.asDiagonal() * Br[i]);
    d2.block(n_f      , 0       , n_r, n_f)  =  weight * exp( - log_sigma2_noise) * (Br[i].transpose() * iV.asDiagonal() * Bf[i]);
  }
  if(Bf.size() > 0 )
    d2.block(0      , 0     , n_f, n_f)  =  weight * exp( - log_sigma2_noise) * (Bf[i].transpose() * iV.asDiagonal() * Bf[i]);

  if(Br.size()>0){
    // dbeta_r dsigma
    d2.block(n_f              , n_s + 2 * n_r +n_f + 1, n_r , 1 )   =  2 * weight * exp( - 1.5 * log_sigma2_noise)  * (Br[i].transpose() * iV.asDiagonal() * res_);
    d2.block(n_s + 2 * n_r + n_f  + 1, n_f              , 1   , n_r ) = d2.block(n_f              , n_s + 2 * n_r +n_f + 1, n_r , 1 ).transpose();
  }

 if(Bf.size() > 0){
  // dbeta_r dsigma
    d2.block(0            , n_s + 2 * n_r + n_f + 1 , n_f , 1 )   =  2 * weight * exp( - 1.5 * log_sigma2_noise)  * (Bf[i].transpose() * iV.asDiagonal()* res_);
    d2.block(n_s + 2 * n_r + n_f + 1, 0             , 1   , n_f ) = d2.block(0            , n_s + 2 * n_r + n_f  + 1, n_f , 1 ).transpose();
  }
  d2(n_s +  2 * n_r +n_f , n_s +  2 * n_r +n_f ) =  weight * 0.5 / pow(nu,2);
  d2(n_s +  2 * n_r +n_f + 1, n_s +  2 * n_r +n_f + 1) =  3  * weight * exp( - 2   * log_sigma2_noise)  * res_.dot(iV.cwiseProduct(res_));
  d2(n_s +  2 * n_r +n_f + 1, n_s +  2 * n_r +n_f + 1) +=  -1 * weight * res_.size()  * exp( - log_sigma2_noise);
  return(d2);
}


void NIGMixedEffect::gradient2(const int i,
                                 const Eigen::VectorXd& res,
                                 const Eigen::VectorXd& iV,
                                 const double log_sigma2_noise,  // = 0
                                 const double EiV, // = 0
                                 const double weight, //  = 1
                                 const int use_EU // =1
                                 )
{
   Eigen::VectorXd res_  = res;
    if(Br.size() > 0){

      Eigen::VectorXd U_ = U.col(i) - (-1 + V(i)) * mu;
      gradient_sigma(i, U_, weight);
      if(use_EU)
      {
        res_ -= Br[i] * EU;
        U_ = EU - (-1 + V(i)) * mu;
      }else{
        res_ -= Br[i] * U.col(i);
      }


      grad_beta_r  += weight * exp( - log_sigma2_noise) * (Br[i].transpose() *  iV.cwiseProduct(res_));
      grad_beta_r2 += weight *  (invSigma * U_)/V(i);
      H_beta_random +=  exp( - log_sigma2_noise) * (Br[i].transpose() * iV.asDiagonal() * Br[i]);



      gradMu   += weight * ((-1 + V(i) )/V(i) ) * (invSigma * U_);
      gradMu_2 += weight * (-1 + V(i) ) * exp( - log_sigma2_noise) * (Br[i].transpose() * res_);
      term1_mu += weight * ((-1 + V(i) )/V(i) )*(-1 + V(i) );
      term2_mu += weight * ((-1 + V(i))/V(i) )*(invSigma*U_) + weight*(-1 + V(i))*exp(-log_sigma2_noise)*(Br[i].transpose()*res_);
      // dnu
      grad_nu += weight * 0.5 * (1. / nu - V(i) - 1. / V(i) + 2. );
      term1   += weight * (V(i) + 1. / V(i) - 2.);
      term2   += weight * 1.;
    }
    if(Bf.size() > 0){
      grad_beta_f   += weight * exp( - log_sigma2_noise) * (Bf[i].transpose() *  iV.cwiseProduct(res_));
      H_beta_fixed  += weight * exp( - log_sigma2_noise) * (Bf[i].transpose() * iV.asDiagonal() * Bf[i]);
    }
    counter++;
  weight_total += weight;
}

void NIGMixedEffect::gradient_sigma(const int i, Eigen::VectorXd& U_ ,const double weight)
{
  Eigen::MatrixXd UUT =  (U_ * U_.transpose());
  UUT.array() /= V(i);
  UUt    += weight * vec( UUT);
}
void NIGMixedEffect::step_theta(const double stepsize,
								const double learning_rate,
								const double polyak_rate,
								const int burnin)
{
  if(Br.size() > 0){
    step_beta_random(stepsize, 0,burnin);
    step_mu(stepsize, 0,burnin);
    step_nu(stepsize, learning_rate,burnin);
    step_Sigma(stepsize, 0,burnin);
    a_GIG = mu.transpose() * (invSigma * mu);
    a_GIG += nu;
    H_beta_random.setZero(Br[0].cols(), Br[0].cols());
  }

  if(Bf.size() > 0)
    step_beta_fixed(stepsize, learning_rate,burnin);

  counter = 0;
  weight_total = 0;

  clear_gradient();

  if(store_param){
    if(Bf.size() > 0){
      if(vec_counter == 0 || polyak_rate == -1){
        betaf_vec.row(vec_counter)  = beta_fixed;
      } else {
        betaf_vec.row(vec_counter)  = polyak_rate * beta_fixed;
        betaf_vec.row(vec_counter) += (1 - polyak_rate) * betaf_vec.row(vec_counter - 1);
        /*
        Rcpp::Rcout << "----" << polyak_rate << "\n";
        Rcpp::Rcout << beta_fixed.transpose() << "\n";
        Rcpp::Rcout << betaf_vec.row(vec_counter - 1)  << "\n";
        Rcpp::Rcout << betaf_vec.row(vec_counter)  << "\n";
        Rcpp::Rcout << "----" << "\n";
        */
      }
    }
    if(Br.size() > 0)
    {
      betar_vec.row(vec_counter)  = beta_random;
  	  mu_vec.row(vec_counter) = mu;
		  nu_vec[vec_counter] = nu;
      Eigen::Map<Eigen::VectorXd> temp(Sigma.data(),Sigma.size());
      Sigma_vec.row(vec_counter)  = temp;
    }
    if(Br.size() > 0)
    {
      	Eigen::Map<Eigen::VectorXd> temp(Sigma.data(),Sigma.size());
    	if(vec_counter == 0 || polyak_rate == -1){
      		betar_vec.row(vec_counter)  = beta_random;
      		Sigma_vec.row(vec_counter)  = temp;
      		mu_vec.row(vec_counter) = mu;
      		nu_vec[vec_counter] = nu;
      	}else{
      		betar_vec.row(vec_counter).array()  = polyak_rate * beta_random.array();
      		betar_vec.row(vec_counter).array()  +=  (1 - polyak_rate) * betar_vec.row(vec_counter - 1).array();
      		Sigma_vec.row(vec_counter).array()  = polyak_rate * temp.array();
      		Sigma_vec.row(vec_counter).array()  += (1 - polyak_rate) * Sigma_vec.row(vec_counter - 1).array();

      		mu_vec.row(vec_counter)          = polyak_rate * mu.array();
      		mu_vec.row(vec_counter).array() +=  (1 - polyak_rate) * mu_vec.row(vec_counter - 1).array();

      		nu_vec[vec_counter]     = polyak_rate * nu + (1 - polyak_rate) * nu_vec[vec_counter - 1];
      	}
    }
    vec_counter++;
  }


}
void NIGMixedEffect::step_Sigma(const double stepsize, const double learning_rate,const int burnin)
{
  double pos_def = 0;

  UUt -= weight_total * vec(Sigma);
  dSigma_vech = 0.5 * Dd.transpose() * iSkroniS * UUt;
  ddSigma = 0.5 * weight_total * Dd.transpose() * iSkroniS * Dd;
  dSigma_vech = ddSigma.ldlt().solve(dSigma_vech);
  dSigma_vech_old.array() *= learning_rate;
  dSigma_vech_old.array() += dSigma_vech.array();

  double stepsize_temp  = stepsize;
  while(pos_def <= 0){
    Eigen::VectorXd Sigma_vech_temp = Sigma_vech;
    Sigma_vech_temp.array() += stepsize_temp * dSigma_vech_old.array();
    Eigen::VectorXd temp = Dd*Sigma_vech_temp;
    Sigma = veci(temp, Sigma.rows(), Sigma.cols());
    stepsize_temp *= 0.5;
    SelfAdjointEigenSolver<MatrixXd> eig(Sigma,EigenvaluesOnly);
    pos_def = eig.eigenvalues().minCoeff();
    if(stepsize_temp <= 1e-16){
        Rcpp::Rcout << "Sigma = \n" << Sigma << "\n";
        Rcpp::Rcout << "pos_def = " << pos_def <<"\n";
        throw("in midexeffect not pos def \n");
    }



  }
  SelfAdjointEigenSolver<MatrixXd> eig(Sigma,EigenvaluesOnly);
  pos_def = eig.eigenvalues().minCoeff();
  //Rcpp::Rcout << "Sigma = " << Sigma  << "\n";
  //Rcpp::Rcout << "eig = " << eig.eigenvalues() << "\n";

  if(pos_def <= 1e-6){
  		dSigma_vech_old *= 0;
  		gradMu_old      *= 0;
  }

  if(pos_def <= 1e-14){
		Rcpp::Rcout << "NIGMixedEffect:: Sigma almost singular \n" ;
		throw("error");

	}

  if(pos_def/mu.array().abs().maxCoeff() < 0.0005)
    sample_MALA = 1;
  else
    sample_MALA = 0;

    UUt.setZero(Sigma.cols() * Sigma.rows());
    invSigma  = Sigma.inverse();
    iSkroniS = kroneckerProduct(invSigma, invSigma);
    Sigma_vech = vech(Sigma);
}

void NIGMixedEffect::step_mu(const double stepsize, const double learning_rate,const int burnin)
{
	gradMu_old.array() *= learning_rate;
  gradMu_old += 0.5 *  H_beta_random.ldlt().solve(gradMu) / VV;
  // H_beta_random = H_mu_random
  gradMu_old += 0.5 * (Sigma * gradMu_2)/ (weight_total * (2*EiV - EV));

  if(burnin == 1){
    mu = (1/term1_mu) * (Sigma*term2_mu);
  } else {
    mu += stepsize * gradMu_old;
  }

  gradMu_2.setZero(Br[0].cols(), 1);
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

  if(nu_temp  > 100){
  		nu_temp =100;
      dnu_old = 0;
  } else if(nu_temp < 5e-06){
    nu_temp =5e-06;
    dnu_old = 0;
  }

	nu = nu_temp;

  EiV = 1. + 1./nu;
  VV = 1./nu;
  set_b_GIG();
  set_a_GIG();
}
void NIGMixedEffect::step_beta_fixed(const double stepsize, const double learning_rate,const int burnin)
{
  if(beta_fixed_constrainted.sum() > 0){
	 dbeta_f_old.array() *= learning_rate;


    solve_const_x_Ab(dbeta_f_old, 
                   beta_fixed_constrainted,
                   grad_beta_f,
                   H_beta_fixed);
    beta_fixed += stepsize *  dbeta_f_old;
  }
  H_beta_fixed.setZero(Bf[0].cols(), Bf[0].cols());
}
void NIGMixedEffect::step_beta_random(const double stepsize, const double learning_rate,const int burnin)
{

  if(beta_random_constrainted.sum() > 0){
    dbeta_r_old.array() *= learning_rate;

    solve_const_x_Ab(dbeta_r_old, 
                   beta_random_constrainted,
                   0.5 * grad_beta_r,
                   H_beta_random);
	 dbeta_r_old += 0.5 * (Sigma * grad_beta_r2)/ (weight_total * EiV);
	 beta_random += stepsize * dbeta_r_old;
  }
	grad_beta_r2.setZero(Br[0].cols());
  H_beta_random.setZero(Br[0].cols(), Br[0].cols());
  grad_beta_r2.setZero(Br[0].cols());
}




void NIGMixedEffect::clear_gradient()
{
  term1 = 0.;
  term2 = 0.;
	if(Bf.size() > 0)
		grad_beta_f.setZero(Bf[0].cols());

	if(Br.size() > 0){
		dSigma_vech.setZero(Sigma_vech.size());
		grad_beta_r.setZero(Br[0].cols());
		gradMu.setZero(Br[0].cols(), 1);
		grad_nu = 0;
    UUt.array() *= 0;
	}

  weight_total = 0;

}


Eigen::VectorXd NIGMixedEffect::get_gradient()
{
	Eigen::VectorXd g(npars);
	int start = 0;
	if(Bf.size() > 0 ){
		g.segment(0, Bf[0].cols()) = grad_beta_f;
		start += Bf[0].cols();
	}
	if(Br.size() > 0)
	{
		g.segment(start, Br[0].cols()) = grad_beta_r;
		start += Br[0].cols();

		g.segment(start, Br[0].cols()) = gradMu;
		start += Br[0].cols();


		Eigen::MatrixXd UUt_temp = UUt;
		UUt_temp -= weight_total * vec(Sigma);
  	dSigma_vech = 0.5 * Dd.transpose() * iSkroniS  * UUt_temp;
		g.segment(start, dSigma_vech.size()) = dSigma_vech;
    start += dSigma_vech.size();
    g[start] = grad_nu;
	}
	return(g);
}
