#ifndef __PROCESS__
#define __PROCESS__


#include <vector>
#include <Rcpp.h>
#include <RcppEigen.h>
#include "operatorMatrix.h"
#include "solver.h"
#include "GIG.h"
#include "rgig.h"

//object for storing the processes object for several classes
class Process {

  public:
	int npars;
  	int store_param;
  	int nindv; // number indiviuals
  	std::vector< Eigen::VectorXd > Xs;
  	std::vector< Eigen::VectorXd >  Vs;
  	Eigen::SparseMatrix<double,0,int>  Q;
  	Eigen::VectorXd  h;
  	Eigen::VectorXd  iV;
  	std::string type_process;
    Process() {};
    Eigen::MatrixXd Cov_theta;// assymptotic covariance of the parameters
    
    virtual Eigen::VectorXd  get_gradient() { Eigen::VectorXd temp; return(temp);};
    virtual void  clear_gradient() {};


	//print iteration data
    virtual void printIter(){};
    // setups to store the tracjetory
    virtual void setupStoreTracj(const int Niter){};
    virtual ~Process(){};
    virtual Rcpp::List toList() {};
    virtual void gradient( const int i ,
			   			  const Eigen::SparseMatrix<double,0,int> & K,
			   			  const Eigen::SparseMatrix<double,0,int> & A,
			   			  const Eigen::VectorXd& res,
			   			  const double sigma,
			   			  const double trace_var){};
    virtual void step_theta(const double step){};
    virtual void sample_V(const int,
    					            gig &,
                          const Eigen::SparseMatrix<double,0,int> &){};

	
	virtual void simulate_V(const int,
    			  gig &){};
    virtual void initFromList(const Rcpp::List &, const Eigen::VectorXd  &){};
    // sampling where the measurement noise is normal
    virtual void sample_X(const int i,
              Eigen::VectorXd & Z,
              const Eigen::VectorXd & Y,
              const Eigen::SparseMatrix<double,0,int> & Q,
              const Eigen::SparseMatrix<double,0,int> & K,
              const Eigen::SparseMatrix<double,0,int> & A,
              const double sigma,
              cholesky_solver       & solver){};
    // sampling where the measurement noise is non-normal
    virtual void sample_Xv2(  const int i,
              Eigen::VectorXd & Z,
              const Eigen::VectorXd & Y,
              const Eigen::SparseMatrix<double,0,int> & Q,
              const Eigen::SparseMatrix<double,0,int> & K,
              const Eigen::SparseMatrix<double,0,int> & A,
              const double sigma,
              cholesky_solver       & solver,
              const Eigen::VectorXd & iV_noise ){};



    virtual void gradient_v2( const int i ,
			   			  const Eigen::SparseMatrix<double,0,int> & K,
			   			  const Eigen::SparseMatrix<double,0,int> & A,
			   			  const Eigen::VectorXd& res,
			   			  const double sigma,
			   			  const Eigen::VectorXd& iV_noise,
			   			  const double EiV_noise,
			   			  const double trace_var) {};

	//simulate from prior distribution
    virtual  void simulate(const int ,
			  Eigen::VectorXd & ,
			  const Eigen::SparseMatrix<double,0,int> & , 
              const Eigen::SparseMatrix<double,0,int> & ,
			  Eigen::VectorXd& ,
			 cholesky_solver  & ) = 0;
	
	
	/*
    	stores the covariance of the parameters 
    */
	void set_covariance(const Eigen::MatrixXd & Cov_in) {Cov_theta = Cov_in;};
};

class GaussianProcess : public Process{

	void initFromList(const Rcpp::List  &,const  Eigen::VectorXd &);
	void sample_X( const int i,
              Eigen::VectorXd & Z,
              const Eigen::VectorXd & Y,
              const Eigen::SparseMatrix<double,0,int> & Q,
              const Eigen::SparseMatrix<double,0,int> & K,
              const Eigen::SparseMatrix<double,0,int> & A,
              const double sigma,
              cholesky_solver       & solver);
    void sample_Xv2(  const int i,
              Eigen::VectorXd & Z,
              const Eigen::VectorXd & Y,
              const Eigen::SparseMatrix<double,0,int> & Q,
              const Eigen::SparseMatrix<double,0,int> & K,
              const Eigen::SparseMatrix<double,0,int> & A,
              const double sigma,
              cholesky_solver       & solver,
              const Eigen::VectorXd & iV_noise);

    Rcpp::List toList();
    
    
	//simulate from prior distribution
    void  simulate(const int ,
			  Eigen::VectorXd & ,
			  const Eigen::SparseMatrix<double,0,int> & , 
              const Eigen::SparseMatrix<double,0,int> & ,
			  Eigen::VectorXd&,
			  cholesky_solver  &  );
};



class GHProcess : public Process{


	private:

  		Eigen::VectorXd  h2;
  		double h_sum;
  		double h_min;
  		double h_digamma;
  		double h_trigamma;
		double dmu ;
		double dnu, ddnu ;
		Eigen::VectorXd EV;
		Eigen::VectorXd EiV;
		Eigen::VectorXd mu_vec;
		Eigen::VectorXd nu_vec;
		int vec_counter;
		double counter;
		double H_mu, ddmu_1, ddmu_2;
		double Vv_mean;
		double h3_mean;

		void update_nu();
		void grad_nu(const int);
		void gradient_mu_centered(const int ,
								  const Eigen::SparseMatrix<double,0,int> & );

	public:


	double mu;
	double nu;
	void initFromList(const Rcpp::List  &, const Eigen::VectorXd &);
	void sample_X( const int i,
              Eigen::VectorXd & Z,
              const Eigen::VectorXd & Y,
              const Eigen::SparseMatrix<double,0,int> & Q,
              const Eigen::SparseMatrix<double,0,int> & K,
              const Eigen::SparseMatrix<double,0,int> & A,
              const double sigma,
              cholesky_solver       & solver);


    void sample_Xv2(  const int i,
              Eigen::VectorXd & Z,
              const Eigen::VectorXd & Y,
              const Eigen::SparseMatrix<double,0,int> & Q,
              const Eigen::SparseMatrix<double,0,int> & K,
              const Eigen::SparseMatrix<double,0,int> & A,
              const double sigma,
              cholesky_solver       & solver,
              const Eigen::VectorXd & iV_noise);

    void gradient( const int i ,
			   			  const Eigen::SparseMatrix<double,0,int> & K,
			   			  const Eigen::SparseMatrix<double,0,int> & A,
			   			  const Eigen::VectorXd& res,
			   			  const double sigma,
			   			  const double trace_var);
	void gradient_v2( const int i ,
			   			  const Eigen::SparseMatrix<double,0,int> & K,
			   			  const Eigen::SparseMatrix<double,0,int> & A,
			   			  const Eigen::VectorXd& res,
			   			  const double sigma,
			   			  const Eigen::VectorXd& iV_noise,
			   			  const double EiV_noise,
			   			  const double trace_var);
    void step_theta(const double );
    void step_mu(const double );
    void step_nu(const double );
    void printIter();
    Rcpp::List toList();
    void setupStoreTracj(const int);
    void sample_V(const int,
    			  gig &,
                  const Eigen::SparseMatrix<double,0,int> &);
	
	//simulate from prior distribution
    void simulate(const int ,
			  Eigen::VectorXd & ,
			  const Eigen::SparseMatrix<double,0,int> & , 
              const Eigen::SparseMatrix<double,0,int> & ,
			  Eigen::VectorXd& ,
			  cholesky_solver  &  );
	void simulate_V(const int,
    			  gig &);
    			  
    Eigen::VectorXd  get_gradient();
    void  clear_gradient();
};



#endif