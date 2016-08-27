#ifndef __MEFF__
#define __MEFF__
#include <vector>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <random>
#include <chrono>
#include "MatrixAlgebra.h"
#include "GIG.h"
//TODO: fixed and random effects should be updated jointly if the are very correclated!!
class MixedEffect {

  protected:

  public:
  	int store_param; // store the parameter to list
  	int vec_counter; // internal parameter counter
    virtual void printIter(){}; //print iteration data
    virtual void setupStoreTracj(const int Niter){}; // setups to store the tracjetory
  	
  	
  	int npars; // number of parameters
  	
  	
  	Eigen::MatrixXd Cov_theta;// assymptotic covariance of the parameters
    Eigen::MatrixXd Sigma;
    int Sigma_epsilon;  // added for epsilon of normal noise to U to ensure Sigma pos def
    std::string noise;
    Eigen::MatrixXd U;
    std::vector< Eigen::MatrixXd > Bf; // fixed covariates
    std::vector< Eigen::MatrixXd > Br; // mixed covariates
    Eigen::VectorXd beta_random;
    Eigen::VectorXd beta_fixed;
    virtual void initFromList(Rcpp::List const &)=0;
    virtual Rcpp::List toList()=0;
    virtual void sampleU(const int, const Eigen::VectorXd &,  const double ) = 0;
    // sampleU2 is sampling with diagonal covariance matrix for the noise
    virtual void sampleU2(const int,
      					  const Eigen::VectorXd &,
      					  const Eigen::VectorXd &,
      					  const double  ) = 0;
    virtual void remove_cov(const int , Eigen::VectorXd & )  = 0;
    virtual void add_cov(const int    , Eigen::VectorXd & )  = 0;
    virtual void add_inter(const int, Eigen::VectorXd &)     = 0;
    virtual void remove_inter(const int, Eigen::VectorXd &)  = 0;
    // gradient for fixed variance noise
    virtual void gradient(const int , const Eigen::VectorXd&, const double ) = 0;
    // gradient for variable variance noise
    virtual void gradient2(const int ,
    					   const Eigen::VectorXd&,
    					   const Eigen::VectorXd& ,
    					   const double,
    					   const double) = 0;
    					   
    // returns the gradient of all the parameters					   
    virtual Eigen::VectorXd get_gradient() = 0;
    virtual void step_theta(double stepsize) = 0;
    /*
    	simulates from the prior distribution
		putting into Y 
    */
    virtual void simulate()                                 = 0;
    virtual void simulate(std::vector< Eigen::VectorXd > &) = 0;
    virtual void simulate(Eigen::VectorXd  & ,const int )   = 0;
    /*
    	clear gradient
    */
	virtual void clear_gradient() = 0;
    /*
    	stores the covariance of the parameters 
    */
	void set_covariance(const Eigen::MatrixXd & Cov_in) {Cov_theta = Cov_in;};

};

class NormalMixedEffect  : public MixedEffect{
  private:
    Eigen::MatrixXd invSigma;
    Eigen::MatrixXd iSkroniS; // helper matrix
    int counter;
    Eigen::VectorXd UUt;
    Eigen::VectorXd dSigma_vech;
    Eigen::MatrixXd ddSigma;
    Eigen::VectorXd Sigma_vech;
    
    
    Eigen::MatrixXd betaf_vec;
    Eigen::MatrixXd betar_vec;
    Eigen::MatrixXd Sigma_vec;

    Eigen::VectorXd grad_beta_r; // gradient for random intercept
    Eigen::VectorXd grad_beta_r2; //second gradient for random intercept
    Eigen::VectorXd grad_beta_f; // gradient for fixed intercept
    Eigen::MatrixXd H_beta_random; // obsereved fisher infromation for random effect
    Eigen::MatrixXd H_beta_fixed;// obsereved fisher infromation for fixed effect
  public:

    Eigen::MatrixXi D;
    Eigen::MatrixXd Dd;

    void printIter(); //print iteration data
    void setupStoreTracj(const int Niter); // setups to store the tracjetory

    NormalMixedEffect();
    void initFromList(Rcpp::List const &);

    void remove_inter(const int i, Eigen::VectorXd & Y) {Y -= Br[i]*U.col(i);} ;
    void add_inter(const int i, Eigen::VectorXd & Y)    {Y += Br[i]*U.col(i);} ;
    void remove_cov(const int , Eigen::VectorXd & );
    void add_cov(const int    , Eigen::VectorXd & );
    /* computes gradient for the parameters
    	@param index of individual
    	@param residuals
    	@param log_sigma2_noise (logarithm of the measurement error)
	*/
    void gradient(const int  , const Eigen::VectorXd&, const double );
    void gradient2(const int i,
    			   const Eigen::VectorXd& res,
    			   const Eigen::VectorXd& iV,
    			   const double log_sigma2_noise = 0,
    			   const double EiV = 1.);

    void sampleU(const int, const Eigen::VectorXd &, const double ) ;
    virtual void sampleU2(const int i,
      					  const Eigen::VectorXd & res,
      					  const Eigen::VectorXd & iV,
      					  const double log_sigma2 = 0);


    void step_theta(double stepsize);
    void step_Sigma(double stepsize);
    void step_beta_fixed(double stepsize);
    void step_beta_random(double stepsize);
    void simulate();
    void simulate(std::vector< Eigen::VectorXd >  &);
    void simulate(Eigen::VectorXd  & ,const int );

    Rcpp::List toList();
    
    
    /*
    	clear gradient
    */
	void clear_gradient();
	
	
    // returns the gradient of all the parameters					   
    Eigen::VectorXd get_gradient();

};



class NIGMixedEffect  : public MixedEffect{
  private:
    Eigen::MatrixXd invSigma;
    Eigen::MatrixXd iSkroniS; // helper matrix
    int counter;
    Eigen::VectorXd UUt;
    Eigen::VectorXd dSigma_vech;
    Eigen::MatrixXd ddSigma;
    Eigen::VectorXd Sigma_vech;
    
    
    Eigen::MatrixXd betaf_vec;
    Eigen::MatrixXd betar_vec;
    Eigen::MatrixXd mu_vec;
    Eigen::VectorXd nu_vec;
    Eigen::MatrixXd Sigma_vec;
    


    double  grad_nu; // gradient for shape parameter
    Eigen::VectorXd gradMu; // gradient for skewness
    Eigen::VectorXd gradMu_2;// second gradient for skewness
    Eigen::VectorXd grad_beta_r; // gradient for random intercept
    Eigen::VectorXd grad_beta_r2; //second gradient for random intercept
    Eigen::VectorXd grad_beta_f; // gradient for fixed intercept
    Eigen::MatrixXd H_beta_random; // obsereved fisher infromation for random effect
    Eigen::MatrixXd H_beta_fixed;// obsereved fisher infromation for fixed effect
    double EV; //  prior expecation  of V, used for the Hessian of random effects
    double EiV; // prior expectation of 1/V used for the Hessian of random effects
    double VV; // prior variance of V used for the Hessian of random effects
    double a_GIG;
    gig rgig;

    void step_Sigma(double stepsize);
    void step_mu(double stepsize);
    void step_nu(double stepsize);

  public:
    Eigen::MatrixXi D;
    Eigen::MatrixXd Dd;
    Eigen::VectorXd mu;
    double          nu;
    Eigen::VectorXd V;

    NIGMixedEffect();
    void sampleV(const int);
    void initFromList(Rcpp::List const &);
    void sampleU(const int, const Eigen::VectorXd &, const double) ;
    virtual void sampleU2(const int i,
      					  const Eigen::VectorXd & res,
      					  const Eigen::VectorXd & iV,
      					  const double log_sigma2_noise = 0);
    void remove_inter(const int i, Eigen::VectorXd & Y) {Y -= Br[i]*U.col(i);} ;
    void add_inter(const int i, Eigen::VectorXd & Y)    {Y += Br[i]*U.col(i);} ;
    void remove_cov(const int , Eigen::VectorXd & );
    void add_cov(const int    , Eigen::VectorXd & );
    void gradient2(const int i,
    			   const Eigen::VectorXd& res,
    			   const Eigen::VectorXd& iV,
    			   const double log_sigma2_noise = 0,
    			   const double EiV = 1.);
    void gradient(const int , const Eigen::VectorXd& , const double );
    void gradient_sigma(const int , Eigen::VectorXd& );
    void step_theta(double stepsize);
    void step_beta_fixed(double stepsize);
    void step_beta_random(double stepsize);
    Rcpp::List toList();
    void simulate();
    void simulate(std::vector< Eigen::VectorXd > & );
    void simulate(Eigen::VectorXd  & ,const int );
    
    virtual void printIter(); //print iteration data
    virtual void setupStoreTracj(const int Niter); // setups to store the tracjetory


    /*
    	clear gradient
    */
	void clear_gradient();
	
	
    // returns the gradient of all the parameters					   
    Eigen::VectorXd get_gradient();

};

// mixed effect util

// sampling NormalCanoncial N_c(b, Q^) \propto e^{-x'Qx + bx}
Eigen::VectorXd sample_Nc(const Eigen::VectorXd & ,const  Eigen::MatrixXd & );
#endif