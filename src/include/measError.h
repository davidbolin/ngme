#ifndef __MEASE__
#define __MEASE__
#include <vector>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <random>
#include <chrono>
#include "MatrixAlgebra.h"
#include "GIG.h"
class MeasurementError {
  protected:

  public:

  	int store_param; // store the parameter to list
  	int vec_counter; // internal parameter counter
    virtual void printIter(){}; //print iteration data
    virtual void setupStoreTracj(const int Niter) = 0; // setups to store the tracjetory
    virtual ~MeasurementError(){};
    MeasurementError() {sigmas.resize(1);};
	int npars; // number of parameters
  	Eigen::MatrixXd Cov_theta;// assymptotic covariance of the parameters
   	double EV;  // if there the random variance in the Noise E[V]
    double EiV; // if there is random varoance in the noise E[V^-1]
  	double sigma;

  	int nsSigma; // non stationary sigma

  	std::vector< Eigen::VectorXd > sigmas; // if sigma is non stationary stored in sigmas
  	std::vector< Eigen::VectorXd > sigmas_pred;
  	std::vector< Eigen::VectorXd > Vs;
    std::string noise;
    virtual void get_param(std::vector<double> & param_in){param_in.push_back(sigma); };
    virtual void get_param_names(Rcpp::StringVector & names){names.push_back("sigma_error"); };
    virtual void gradient(const int ,
    					  const Eigen::VectorXd& ,
    					  const double) = 0;
    virtual void step_theta(const double stepsize,const double learning_rate = 0,
    						const double polyak_rate = -1,
    						const int burnin = 0) = 0;
    virtual void initFromList(Rcpp::List const &)=0;
    virtual Rcpp::List toList()=0;
    virtual void sampleV(const int , const Eigen::VectorXd&, int = -1) {};

    // sampling from the prior model
  	virtual std::vector< Eigen::VectorXd > simulate( const std::vector< Eigen::VectorXd > )  = 0;

  	virtual Eigen::VectorXd  simulate( const Eigen::VectorXd &)  = 0;
  	virtual Eigen::VectorXd  simulate_par( const Eigen::VectorXd &, std::mt19937 &)  = 0;
  	
  	// Sample the entire noise vector for the i:th subject
  	//inputs subject number, random number generator, size of vector
  	virtual Eigen::VectorXd  simulate_par( const int, std::mt19937 &, const int)  = 0;
    /*
    	clear gradient
    */
	virtual void clear_gradient() = 0;
    /*
    	stores the covariance of the parameters
    */

    virtual Eigen::MatrixXd d2Given(const int ,
    					  const Eigen::VectorXd& ,
    					  const double){return(Eigen::MatrixXd::Zero(0,0));}; //computes second derivates given, latent data
                                                                              // fixed variance

	void set_covariance(const Eigen::MatrixXd & Cov_in) {Cov_theta = Cov_in;};

    /*
     returns the gradient of all the parameters
      */
    virtual Eigen::VectorXd get_gradient()  =0;

    /*
		Function active only for varaince-mixture error
    */
    
    virtual void remove_asym(const int i, Eigen::VectorXd & Y) {};
    virtual void add_asym(const int i, Eigen::VectorXd & Y)    {};

    /*
      Get the diagonal of the covaraince matrix of the measrurment error

    */
    virtual Eigen::VectorXd  getSigma(const int i,  const int n = 0) = 0; 
};


class GaussianMeasurementError : public MeasurementError{
	private:
		double dsigma;
		double ddsigma;
		double dsigma_old;
    	double counter;
    	Eigen::VectorXd sigma_vec;

	public:
		GaussianMeasurementError();
	  ~GaussianMeasurementError(){};
		void gradient(const int ,
					  const Eigen::VectorXd&,
					  const double);
		void step_theta(const double stepsize,
						const double learning_rate = 0,
						const double polyak_rate = -1,
						const int burnin = 0);
		void initFromList(Rcpp::List const &);
		void sampleV(const int i, const Eigen::VectorXd& res, int n_s = -1) {};
		Rcpp::List toList();
		std::vector< Eigen::VectorXd > simulate( const std::vector< Eigen::VectorXd >);
		Eigen::VectorXd  simulate( const Eigen::VectorXd &);
		Eigen::VectorXd  simulate_par( const Eigen::VectorXd &, std::mt19937 &);
		Eigen::VectorXd  simulate_par( const int, std::mt19937 &, const int);

		void clear_gradient();

		Eigen::VectorXd get_gradient();

		void printIter(); //print iteration data
        void setupStoreTracj(const int ); // setups to store the tracjetory

    Eigen::VectorXd  getSigma(const int i,  const int n = 0)  { Eigen::VectorXd Res; 
                                                                Res.setOnes(n);
                                                                Res *= pow(sigma, 2);
                                                                return(Res);}; 
};

class nsGaussianMeasurementError : public MeasurementError{
private:
  Eigen::VectorXd dtheta, theta;
  Eigen::MatrixXd ddtheta;
  Eigen::VectorXd dtheta_old, step;
  std::vector< Eigen::MatrixXd > B;
  std::vector< Eigen::MatrixXd > Bpred;
  double counter;
  Eigen::MatrixXd theta_vec;
  int nrep;
  Rcpp::List out_list;
  int predict;
public:
  nsGaussianMeasurementError();
  ~nsGaussianMeasurementError(){};
  void gradient(const int ,
                const Eigen::VectorXd&,
                const double);
  void step_theta(const double stepsize,
                  const double learning_rate = 0,
                  const double polyak_rate = -1,
                  const int burnin = 0);
  void initFromList(Rcpp::List const &);
  void sampleV(const int i, const Eigen::VectorXd& res, int n_s = -1) {};
  Rcpp::List toList();
  std::vector< Eigen::VectorXd > simulate( const std::vector< Eigen::VectorXd >);
  Eigen::VectorXd  simulate( const Eigen::VectorXd &);
  Eigen::VectorXd  simulate_par( const Eigen::VectorXd &, std::mt19937 &);
  Eigen::VectorXd  simulate_par( const int, std::mt19937 &, const int);
  
  void clear_gradient();
  
  Eigen::VectorXd get_gradient();
  
  void printIter(); //print iteration data
  void setupStoreTracj(const int ); // setups to store the tracjetory
  
    Eigen::VectorXd  getSigma(const int i,  const int n = 0)  { return(sigmas[i].array().square());}; 
  
};


class NormalVarianceMixtureBaseError : public MeasurementError{

	public:
		double dsigma;
		double ddsigma;
    	double dsigma_old;
		gig rgig;
    	double counter;
    	Eigen::VectorXd sigma_vec;
    	Eigen::VectorXd nu_vec;

    	int assymetric;
    	double mu;
    	double dmu;
    	double dmu_old;
    	double ddmu;
    	Eigen::VectorXd mu_vec;
    	double weight_total;
    	double VV;

	  int common_V;
		double nu;
		NormalVarianceMixtureBaseError();
		~NormalVarianceMixtureBaseError(){};
		void gradient(const int ,
					  const Eigen::VectorXd& ,
					  const double);
		void step_theta(const double stepsize,
						const double learning_rate = 0,
						const double polyak_rate = -1,
						const int burnin = 0);
		void step_sigma(const double , const double,const int, const double);
		void step_mu(const double , const double,const int, const double);
		void initFromList(Rcpp::List const &);
		void sampleV(const int , const Eigen::VectorXd& , int = -1);
		Rcpp::List toList();
		std::vector< Eigen::VectorXd > simulate( const std::vector< Eigen::VectorXd >);

		Eigen::VectorXd  simulate( const Eigen::VectorXd &);
		Eigen::VectorXd  simulate_par( const Eigen::VectorXd &, std::mt19937 &);
		Eigen::VectorXd  simulate_par( const int, std::mt19937 &, const int);
		void clear_gradient();

		Eigen::VectorXd get_gradient();

		void printIter(); //print iteration data
        void setupStoreTracj(const int ); // setups to store the tracjetory
        virtual double simulate_V();
        /*
			sample_V - sampling variance component
					res2_j        - residual value
					n_s           - number of samples 
					mu^2/sigma^2  - assymetric parameter
        */
        virtual double sample_V(const double, const int, const double) {return -1;};
    	virtual void get_param(std::vector<double> & param_in);
        virtual void get_param_names(Rcpp::StringVector & names);


     	void remove_asym(const int i, Eigen::VectorXd & Y) ;
    	void add_asym(const int i, Eigen::VectorXd & Y)    ;
      virtual Eigen::MatrixXd d2Given(const int ,
                    const Eigen::VectorXd& ,
                    const double);


    Eigen::VectorXd  getSigma(const int i,  const int n = 0); 

};

class NIGMeasurementError : public NormalVarianceMixtureBaseError{

	private:
		double dnu_old;
 	public:
 	  double term1,term2,term3;
		double dnu;
		double ddnu;
 		NIGMeasurementError();
 		~NIGMeasurementError(){};
 		void printIter();
 		void setupStoreTracj(const int ) ;
 		Rcpp::List toList();
		void initFromList(Rcpp::List const &);
		double simulate_V();
		double sample_V(const double, const int, const double);
		void gradient(const int ,
					  const Eigen::VectorXd& ,
					  const double );
		void step_nu(const double, const double,const int );
		void step_theta(const double stepsize,
						const double learning_rate = 0,
						const double polyak_rate = -1,
						const int burnin = 0);
		void clear_gradient();
		Eigen::VectorXd get_gradient();

		Eigen::MatrixXd d2Given(const int ,
    					  		const Eigen::VectorXd& ,
    					  		const double);

		void get_param(std::vector<double> & param_in){ NormalVarianceMixtureBaseError::get_param(param_in); param_in.push_back(nu); };
        void get_param_names(Rcpp::StringVector & names){NormalVarianceMixtureBaseError::get_param_names(names); names.push_back("nu_error"); };


};

class IGMeasurementError : public NormalVarianceMixtureBaseError{


 	public:
 		double digamma_nu;
 		double trigamma_nu;
		double beta;
		double dnu;
		double dnu_old;
		double ddnu;
 		IGMeasurementError();
 		~IGMeasurementError(){};
 		void printIter();
 		void setupStoreTracj(const int ) ;
 		Rcpp::List toList();
		void initFromList(Rcpp::List const &);
		double simulate_V();
		double sample_V(const double, const int, const double);
		void gradient(const int ,
					  const Eigen::VectorXd&,
					  const double);
		void step_nu(const double ,const double,const int);
		void step_theta(const double stepsize,
						const double learning_rate = 0,
						const double polyak_rate = -1,
						const int burnin = 0);
		void clear_gradient();
		Eigen::VectorXd get_gradient();

    	void get_param(std::vector<double> & param_in){ NormalVarianceMixtureBaseError::get_param(param_in); param_in.push_back(nu); };
        void get_param_names(Rcpp::StringVector & names){NormalVarianceMixtureBaseError::get_param_names(names); names.push_back("nu_error"); };

        
    Eigen::MatrixXd d2Given(const int ,
                    const Eigen::VectorXd& ,
                    const double);
};


#endif
