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


	int npars; // number of parameters
  	Eigen::MatrixXd Cov_theta;// assymptotic covariance of the parameters
   	double EV;  // if there the random variance in the Noise E[V]
    double EiV; // if there is random varoance in the noise E[V^-1]
  	double sigma;
  	std::vector< Eigen::VectorXd > Vs;
    std::string noise;
    virtual void gradient(const int , const Eigen::VectorXd& ) = 0;
    virtual void step_theta(const double stepsize,const double learning_rate = 0,
    						const double polyak_rate = -1) = 0;
    virtual void initFromList(Rcpp::List const &)=0;
    virtual Rcpp::List toList()=0;
    virtual void sampleV(const int , const Eigen::VectorXd&, int = -1) = 0;

    // sampling from the prior model
  	virtual std::vector< Eigen::VectorXd > simulate( const std::vector< Eigen::VectorXd > )  = 0;

  	virtual Eigen::VectorXd  simulate( const Eigen::VectorXd &)  = 0;
  	virtual Eigen::VectorXd  simulate_par( const Eigen::VectorXd &, std::mt19937 &)  = 0;
    /*
    	clear gradient
    */
	virtual void clear_gradient() = 0;
    /*
    	stores the covariance of the parameters
    */
	void set_covariance(const Eigen::MatrixXd & Cov_in) {Cov_theta = Cov_in;};

    /*
     returns the gradient of all the parameters
      */
    virtual Eigen::VectorXd get_gradient()  =0;

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
		void gradient(const int , const Eigen::VectorXd&);
		void step_theta(const double stepsize,
						const double learning_rate = 0,
						const double polyak_rate = -1);
		void initFromList(Rcpp::List const &);
		void sampleV(const int i, const Eigen::VectorXd& res, int n_s = -1) {};
		Rcpp::List toList();
		std::vector< Eigen::VectorXd > simulate( const std::vector< Eigen::VectorXd >);
		Eigen::VectorXd  simulate( const Eigen::VectorXd &);
		Eigen::VectorXd  simulate_par( const Eigen::VectorXd &, std::mt19937 &);

		void clear_gradient();

		Eigen::VectorXd get_gradient();

		void printIter(); //print iteration data
        void setupStoreTracj(const int ); // setups to store the tracjetory


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


	  int common_V;
		double nu;
		NormalVarianceMixtureBaseError();
		void gradient(const int , const Eigen::VectorXd&);
		void step_theta(const double stepsize,
						const double learning_rate = 0,
						const double polyak_rate = -1);
		void step_sigma(const double , const double);
		void initFromList(Rcpp::List const &);
		void sampleV(const int , const Eigen::VectorXd& , int = -1);
		Rcpp::List toList();
		std::vector< Eigen::VectorXd > simulate( const std::vector< Eigen::VectorXd >);

		Eigen::VectorXd  simulate( const Eigen::VectorXd &);
		Eigen::VectorXd  simulate_par( const Eigen::VectorXd &, std::mt19937 &);
		void clear_gradient();

		Eigen::VectorXd get_gradient();

		void printIter(); //print iteration data
        void setupStoreTracj(const int ); // setups to store the tracjetory
        virtual double simulate_V();
        virtual double sample_V(const double, const int) {return -1;};

};

class NIGMeasurementError : public NormalVarianceMixtureBaseError{

	private:
		double dnu_old;
 	public:
		double dnu;
		double ddnu;
 		NIGMeasurementError();
 		void printIter();
 		void setupStoreTracj(const int ) ;
 		Rcpp::List toList();
		void initFromList(Rcpp::List const &);
		double simulate_V();
		double sample_V(const double, const int);
		void gradient(const int , const Eigen::VectorXd&);
		void step_nu(const double, const double );
		void step_theta(const double stepsize,
						const double learning_rate = 0,
						const double polyak_rate = -1);
		void clear_gradient();
		Eigen::VectorXd get_gradient();

};

class IGMeasurementError : public NormalVarianceMixtureBaseError{


 	public:
 		double digamma_nu;
 		double trigamma_nu;
		double dnu;
		double dnu_old;
		double ddnu;
 		IGMeasurementError();
 		void printIter();
 		void setupStoreTracj(const int ) ;
 		Rcpp::List toList();
		void initFromList(Rcpp::List const &);
		double simulate_V();
		double sample_V(const double, const int);
		void gradient(const int , const Eigen::VectorXd&);
		void step_nu(const double ,const double);
		void step_theta(const double stepsize,
						const double learning_rate = 0,
						const double polyak_rate = -1);
		void clear_gradient();
		Eigen::VectorXd get_gradient();

};


#endif
