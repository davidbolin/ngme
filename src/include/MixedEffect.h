#ifndef __MEFF__
#define __MEFF__
#include <vector>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <random>
#include <chrono>
#include <string>
#include <math.h>
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
    virtual ~MixedEffect() {};


  	int npars; // number of parameters

  	Eigen::MatrixXd Cov_theta;// assymptotic covariance of the parameters
    Eigen::MatrixXd Sigma;
    int Sigma_epsilon;  // added for epsilon of normal noise to U to ensure Sigma pos def
    std::string noise;
    Eigen::MatrixXd U;
    Eigen::MatrixXd EU;
    std::vector< Eigen::MatrixXd > Bf; // fixed covariates
    std::vector< Eigen::MatrixXd > Br; // mixed covariates
    Eigen::VectorXd beta_random_constrainted;
    Eigen::VectorXd beta_random;
    Eigen::VectorXd beta_fixed_constrainted;
    Eigen::VectorXd beta_fixed;
	Eigen::VectorXd dbeta_r_old;
	Eigen::VectorXd dbeta_f_old;
    virtual void initFromList(Rcpp::List const &)=0;
    virtual void get_param( std::vector<double> & param_in)
    {
      if(Bf.size() > 0 )
        {
         for (int i = 0; i < Bf[0].cols(); i++)
            param_in.push_back(beta_fixed[i]);
        }
        if(Br.size() > 0 )
        {
         for (int i = 0; i < Br[0].cols(); i++)
            param_in.push_back(beta_random[i]);
        }


    }
    virtual void get_param_names(Rcpp::StringVector & names){
        if(Bf.size() > 0 )
        {
         for (int i = 0; i < Bf[0].cols(); ++i)
          {
            names.push_back("beta_fixed_" + std::to_string(i+1));
          }
        }
        if(Br.size() > 0 )
        {
          for (int i = 0; i < Br[0].cols(); ++i)
            names.push_back("beta_random_" + std::to_string(i+1));
        }

    };
    virtual Rcpp::List toList()=0;
    virtual void sampleU(const int, const Eigen::VectorXd &,  const double ) = 0;
    virtual void sampleU_par(const int, const Eigen::VectorXd &,  const double, std::mt19937 &) = 0;
    // sampleU2 is sampling with diagonal covariance matrix for the noise
    virtual void sampleU2(const int,
      					  const Eigen::VectorXd &,
      					  const Eigen::VectorXd &,
      					  const double  ) = 0;
    virtual void sampleU2_par(const int,
                          const Eigen::VectorXd &,
                          const Eigen::VectorXd &,
                          std::mt19937 &,
                          const double) = 0;
    void remove_inter(const int i, Eigen::VectorXd & Y) {if(Br.size()>0){
    													  Y -= Br[i]*U.col(i);}};
    void add_inter(const int i, Eigen::VectorXd & Y)    { if(Br.size()>0){
    													 Y += Br[i]*U.col(i);} };

	void remove_cov(const int i, Eigen::VectorXd & Y)
	{
  		if(Br.size() > 0 )
    		Y -= Br[i] * beta_random;
  		if(Bf.size() > 0)
    		Y -= Bf[i] * beta_fixed;
	};
	void add_cov(const int    i, Eigen::VectorXd & Y)
	{
  		if(Br.size() > 0 )
    		Y += Br[i] * beta_random;
  		if(Bf.size() > 0)
    		Y += Bf[i] * beta_fixed;
	};

    // gradient for fixed variance noise
    virtual void gradient(const int ,
    					  const Eigen::VectorXd&,
    					  const double,
    					  const double,
                          const int = 1) = 0;
    // gradient for variable variance noise
    virtual void gradient2(const int ,
    					   const Eigen::VectorXd&,
    					   const Eigen::VectorXd& ,
    					   const double,
    					   const double,
    					   const double,
                           const int = 1) = 0;
    virtual Eigen::MatrixXd d2Given(const int ,
                          const Eigen::VectorXd&,
                          const double,
                          const double ){return(Eigen::MatrixXd::Zero(0,0));}; //computes second derivates given, latent data
                                                                              // fixed variance
    virtual Eigen::MatrixXd d2Given2(const int ,
                           const Eigen::VectorXd&,
                           const Eigen::VectorXd& ,
                           const double,
                           const double,
                           const double){return(Eigen::MatrixXd::Zero(0,0));};  //computes second derivates given, latent data
                                                                               // variable variance


    // returns the gradient of all the parameters
    virtual Eigen::VectorXd get_gradient() = 0;
    virtual void step_theta(const double stepsize,
    						const double  learning_rate = 0,
    						const double polyak_rate = -1,
    						const int burnin = 0) = 0;
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
    Eigen::VectorXd dSigma_vech_old;
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
    Eigen::MatrixXd H_beta;// obsereved fisher infromation for fixed effect
	Eigen::VectorXd grad_beta;
	int n_f, n_r;
    double weight_total;
  public:

    Eigen::MatrixXi D;
    Eigen::MatrixXd Dd;

    void printIter(); //print iteration data
    void setupStoreTracj(const int Niter); // setups to store the tracjetory

    NormalMixedEffect();
    ~NormalMixedEffect();
    void initFromList(Rcpp::List const &);

    void get_param(std::vector<double> &);
    void get_param_names(Rcpp::StringVector & names);

    /* computes gradient for the parameters
    	@param index of individual
    	@param residuals
    	@param log_sigma2_noise (logarithm of the measurement error)
	*/
    void gradient(const int  ,
     			  const Eigen::VectorXd&,
     			  const double,
     			  const double,
                  const int use_EU = 1);
    void gradient2(const int i,
    			   const Eigen::VectorXd& res,
    			   const Eigen::VectorXd& iV,
    			   const double log_sigma2_noise = 0,
    			   const double EiV = 1.,
    			   const double weight = 1.,
             const int use_EU = 1);
    Eigen::MatrixXd d2Given(const int ,
                          const Eigen::VectorXd&,
                          const double,
                          const double ); //computes second derivates given, latent data
                                                                              // fixed variance
    virtual Eigen::MatrixXd d2Given2(const int ,
                           const Eigen::VectorXd&,
                           const Eigen::VectorXd& ,
                           const double,
                           const double,
                           const double);  //computes second derivates given, latent data


    void sampleU(const int, const Eigen::VectorXd &, const double ) ;
    void sampleU_par(const int, const Eigen::VectorXd &,  const double, std::mt19937 &);
    void sampleU2(const int i,
      					  const Eigen::VectorXd & res,
      					  const Eigen::VectorXd & iV,
      					  const double log_sigma2 = 0);
    void sampleU2_par(const int i,
                  const Eigen::VectorXd & res,
                  const Eigen::VectorXd & iV,
                  std::mt19937 & random_engine,
                  const double log_sigma2 = 0);


    void step_theta(const double stepsize,
    				const double learning_Rate  = 0,
    				const double polyak_rate   = -1,
    				const int burnin = 0);
    void step_Sigma(const double, const double,const int);
    void step_beta_fixed(const double, const double,const int );
    void step_beta_random(const double, const double,const int );
    void step_beta(const double ,const double ,const int );
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
    Eigen::VectorXd dSigma_vech_old;


    Eigen::MatrixXd betaf_vec;
    Eigen::MatrixXd betar_vec;
    Eigen::MatrixXd mu_vec;
    Eigen::VectorXd nu_vec;
    Eigen::MatrixXd Sigma_vec;
    double weight_total;


	double dnu_old;
    double  grad_nu; // gradient for shape parameter
    Eigen::VectorXd gradMu; // gradient for skewness
    Eigen::VectorXd gradMu_old;
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

    void step_Sigma(const double, const double,const int);
    void step_mu(const double , const double,const int);
    void step_nu(const double, const double,const int);
    double term1,term2,term1_mu;
    Eigen::VectorXd term2_mu;
  public:

    int count_MALA;
    int accept_MALA;
    int sample_MALA;
    Eigen::MatrixXi D;
    Eigen::MatrixXd Dd;
    Eigen::VectorXd mu;
    double          nu;
    Eigen::VectorXd V;

    NIGMixedEffect();
    ~NIGMixedEffect(){};

    void get_param(std::vector<double> &);
    void get_param_names(Rcpp::StringVector & );
    void sampleV(const int);
    void initFromList(Rcpp::List const &);
    void sampleU(const int, const Eigen::VectorXd &, const double) ;
    void sampleU_MALA(const int, const Eigen::VectorXd &, const double) ;
    void sampleU_Gibbs(const int, const Eigen::VectorXd &, const double) ;
    void sampleU_par(const int, const Eigen::VectorXd &,  const double, std::mt19937 &);
    void sampleU_MALA_(const int ,
                                const Eigen::VectorXd & ,
                                const Eigen::VectorXd & ,
                                const Eigen::MatrixXd   & );
    void sampleU2(const int i,
      					const Eigen::VectorXd & res,
      					const Eigen::VectorXd & iV,
      					const double log_sigma2_noise = 0);
    void sampleU2_Gibbs( const int i,
                        const Eigen::VectorXd & res,
                        const Eigen::VectorXd & iV,
                        const double log_sigma2_noise = 0);
    void sampleU2_MALA( const int i,
                        const Eigen::VectorXd & res,
                        const Eigen::VectorXd & iV,
                        const double log_sigma2_noise = 0);
    void sampleU2_par(const int i,
                  const Eigen::VectorXd & res,
                  const Eigen::VectorXd & iV,
                  std::mt19937 & random_engine,
                  const double log_sigma2_noise = 0);
    void gradient2(const int i,
    			   const Eigen::VectorXd& res,
    			   const Eigen::VectorXd& iV,
    			   const double log_sigma2_noise = 0,
    			   const double EiV = 1.,
    			   const double weight = 1.,
                   const int use_EU = 1);


    Eigen::MatrixXd d2Given(const int ,
                  const Eigen::VectorXd& ,
                  const double ,
                  const double);

    Eigen::MatrixXd d2Given2(const int ,
                           const Eigen::VectorXd&,
                           const Eigen::VectorXd& ,
                           const double,
                           const double,
                           const double);
    void gradient(const int ,
    			  const Eigen::VectorXd& ,
    			  const double ,
    			  const double,
                  const int use_EU = 1);
    void gradient_sigma(const int ,
    					Eigen::VectorXd& ,
    					const double);
    void step_theta(const double stepsize,
    				const double learning_Rate  = 0,
    				const double polyak_rate   = -1,
    				const int burnin = 0);
    void step_beta_fixed(const double stepsize, const double,const int);
    void step_beta_random(const double stepsize, const double,const int);
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

//solve constrained
// solves x[constrained == 1] += A[constrained==1, H_constrained == 1]^{-1} b
void solve_const_x_Ab(Eigen::VectorXd & , const Eigen::VectorXd & , const Eigen::VectorXd & ,const  Eigen::MatrixXd & );

// sampling NormalCanoncial N_c(b, Q^) \propto e^{-x'Qx + bx}
Eigen::VectorXd sample_Nc(const Eigen::VectorXd & ,const  Eigen::MatrixXd & );
Eigen::VectorXd sample_Nc_par(const Eigen::VectorXd & ,const  Eigen::MatrixXd &,std::mt19937 &);
// computing expectation 1/ of the variance component for Normal Generalised inverse Gaussian mixture (E[V^-1]).
// for the posterior distribution
double EiV_NGIG(const Eigen::VectorXd & ,
                const Eigen::MatrixXd & ,
                const Eigen::VectorXd & ,
                const Eigen::VectorXd & ,
                const double ,
                const double ,
                const double );
// d/dU E[V^-1|U]
Eigen::VectorXd dU_EiV_NGIG(const Eigen::VectorXd & ,
                const Eigen::MatrixXd & ,
                const Eigen::VectorXd & ,
                const Eigen::VectorXd & ,
                const double ,
                const double ,
                const double );
void dU_ddU_NIG(
              Eigen::VectorXd &,
              Eigen::MatrixXd &,
        const Eigen::VectorXd &,
        const Eigen::MatrixXd &,
        const Eigen::VectorXd &,
        const Eigen::VectorXd &,
        const double ,
        const double ,
        const double ,
        const Eigen::VectorXd & ,
        const Eigen::MatrixXd & ,
        const Eigen::MatrixXd & );
#endif
