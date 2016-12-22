#include <Rcpp.h>
#include <random>
#include <chrono>
#include "GIG.h"
#include "Rmath.h"
using namespace Rcpp;


double db_EiV_GIG(double p, double a, double b) {
    double sqrt_ab = sqrt(a * b);
    double K1 = R::bessel_k(sqrt_ab, p, 2);
    double K0 = R::bessel_k(sqrt_ab, p+1, 2);
    double sqrt_a_div_b = sqrt(a/b);
    double K0dK1 = K0 / K1;

    double db_EiV = 0;
    db_EiV += - 1. - (p + 1.) * K0dK1 / sqrt_ab;
    db_EiV -=  ( - K0dK1 + (p / sqrt_ab) ) * K0dK1;
    db_EiV *= 0.5 * sqrt_a_div_b *sqrt_a_div_b;
    db_EiV -= 0.5 * K0dK1 * sqrt_a_div_b / b;
    db_EiV += (2 * p) / (b*b);  
    return db_EiV;
}


double EiV_GIG(double p, double a, double b) {
    double sqrt_ab = sqrt(a * b);
    double K1 = R::bessel_k(sqrt_ab, p, 2);
    double K0 = R::bessel_k(sqrt_ab, p+1, 2);
    double sqrt_a_div_b = sqrt(a/b);
    double EiV = K0 / K1;
    EiV    *=  sqrt_a_div_b;
    EiV    -= (2 * p) * 1./b;
  return EiV;
}

double dlambda_V(const double loglambda,
                 const Eigen::VectorXd &V, 
                 const Eigen::VectorXd &h,
                 const int GAL)
{
  double dlambda = 0;
  if(GAL){
  for(int i=0; i < h.size(); i++){
    double h_lambda = exp(loglambda) * h[i];
    //digamma(0.1) = digamma(1.1) - 1/0.1;
    if(h_lambda > 1){
        dlambda -=  h_lambda * R::digamma(h_lambda);
      }else
      {
        dlambda -=  h_lambda * R::digamma(h_lambda + 1) - 1.;
      }
    dlambda += h_lambda *  log(V(i) ) ;
  }
  }else{
    double srqt_two = pow(2, 0.5);
    for(int i=0; i < h.size(); i++){
      dlambda +=  1 -  ( pow(h(i), 2) / V(i) ) * exp( 2 * loglambda);
      dlambda += srqt_two * h(i)  * exp(loglambda);  
    }
      
  }
  
  return(dlambda); 
}

Eigen::VectorXd sampleV_pre(gig &sampler,
                        const Eigen::VectorXd &h, 
                        const double nu,
                        const std::string type)
{
  Eigen::VectorXd V(h.size());
  double Vadj = 1e-14;
  if(type == "NIG"){
    for(int i = 0; i < h.size(); i++)
    	V[i] = sampler.sample(-0.5 , nu, pow(h[i] , 2)* nu);
      
  }else if(type == "GAL"){
    for(int i = 0; i < h.size(); i++)
      V[i] = sampler.sample( h[i] * nu, 2 *nu, 0) + Vadj; 
      
  }else if(type == "CH"){
    for(int i = 0; i < h.size(); i++)
      V[i] = sampler.sample( -0.5, 0, 2*0.25 * pow(h[i], 2)); 
  	
  }else{
  	throw("sampleV_pre type must either be NIG, GAL or CH");
  }
  
  
  return(V);
}


Eigen::VectorXd sampleV_post(gig &sampler,
                        	 const Eigen::VectorXd &h, 
                        	 const Eigen::VectorXd &KX,
                        	 const double sigma,
                        	 const double mu,
                        	 const double nu,
                        	 const std::string type)
{
  Eigen::VectorXd  p, b;
  
  b = (KX +  mu * h) / sigma;
  b = b.array().square();
  double b_adj = 1e-14;
  double a  =  pow(mu / sigma, 2);
  if(type == "GAL"){
    p = h * nu;
    p.array() -= 0.5;
    a += 2 * nu;
    b.array() += b_adj;
  }else if(type == "NIG"){
    p.setOnes(h.size());
    p *= -1.;
    a += pow(nu, 2);
    b.array() += (nu * h).array().square();
  }else if(type == "CH"){
    p.setOnes(h.size());
    p.array() *= -1.;
    b.array() *= 4 * 0.5;
    b.array() += (h ).array().square();
  	b.array() *= 2 * 0.25;
  }else{
  	throw("sampleV_pre type must either be NIG, GAL or CH");
  }
  
  Eigen::VectorXd V(KX.size());
  
  
  
  for(int i = 0; i < KX.size(); i++)
  	V[i] = sampler.sample( p[i], a, b[i] ) ; 
  
   double Vadj  = 1e-13;  
    if(type == "GAL")
    	V.array() += Vadj;
  
  return(V);
}
