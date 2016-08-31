#include <Rcpp.h>
#include <random>
#include <chrono>
#include "GIG.h"
using namespace Rcpp;


double dlambda_V(const double loglambda,
                 const Eigen::VectorXd &V, 
                 const Eigen::VectorXd &h,
                 const int GAL)
{
  double dlambda = 0;
  double Vadj = 1e-12;
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
    dlambda += h_lambda *  log(V(i) - Vadj + 1e-12) ;
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
    	V[i] = sampler.sample(-0.5 , pow(nu, 2), pow(h[i] * nu, 2));
      
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
  
  b = (KX + mu * h) / sigma;
  b = b.array().square();
  double a  =  pow(mu / sigma, 2);
  if(type == "GAL"){
    p = h * nu;
    a += 2 * nu;
    p.array() -= 0.5;
  }else if(type == "NIG"){
    p.setOnes(h.size());
    p *= -1.;
    a += pow(nu, 2);
    b.array() += (h * nu).array().square();
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
  double b_adj = 1e-14;
  double Vadj  = 1e-14;
  for(int i = 0; i < KX.size(); i++)
      V[i] = sampler.sample( p[i], a, b[i] + b_adj) + Vadj; 
  
  return(V);
}
