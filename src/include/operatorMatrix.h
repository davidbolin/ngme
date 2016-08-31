#ifndef __OPERATOR__MATRIX__
#define __OPERATOR__MATRIX__
#include <string.h>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "localprint.h"
#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>
#include <Eigen/OrderingMethods>
#include <Rcpp.h>
#include <RcppEigen.h>
#include "MatrixAlgebra.h"
#include "solver.h"
#ifdef TIME_IT
	#include "TimeIt.h"
#endif

#ifdef _OPENMP
	#include<omp.h>
#endif

class operatorMatrix {
  protected:
    solver * Qsolver;
  public:
  
  
  
  	Eigen::MatrixXd Cov_theta;// assymptotic covariance of the parameters
  	
	int npars; // number of parameters
    Eigen::VectorXd  loc; // location of the position
    operatorMatrix() {Qsolver = NULL;};
    virtual ~operatorMatrix(){delete Qsolver;};
    int d; //dimension
    Eigen::SparseMatrix<double,0,int> Q; // the generic matrix object
    Eigen::MatrixXd K;                   // the generic matrix object if Q is full!
  
  
  
      
    
    virtual Eigen::VectorXd  get_gradient() { Eigen::VectorXd temp; return(temp);};
    virtual void  clear_gradient() {};

  
    virtual void initFromList(Rcpp::List const &)=0;
    virtual void initFromList(Rcpp::List const &, Rcpp::List const &) {Rcpp::Rcout << "initFromList(list1,list2) not implimented in operatorMatrix\n";};

    virtual Rcpp::List output_list() = 0;
    virtual void gradient( const Eigen::VectorXd &, const Eigen::VectorXd & ){};
    virtual void gradient_init( const int, const int){};
    virtual void gradient_add( const Eigen::VectorXd & , 
							   const Eigen::VectorXd & ,
							   const Eigen::VectorXd & ){};
    virtual void step_theta(const double ){};
    virtual void print_parameters( ){};
    virtual double trace_variance( const Eigen::SparseMatrix<double,0,int> & ){return 1;};
    double tau;
    Eigen::VectorXd  tauVec;
    int counter;
    
    
	/*
    	stores the covariance of the parameters 
    */
	void set_covariance(const Eigen::MatrixXd & Cov_in) {Cov_theta = Cov_in;};
};

class constMatrix : public operatorMatrix{
  protected:
    Eigen::SparseMatrix<double,0,int> m;
    Eigen::VectorXd v;
    double m_loc;
    Eigen::VectorXd  h; 
    double h_average;
  public:

	void gradient(const Eigen::VectorXd &, const Eigen::VectorXd & );
  void gradient_init(const int, const int);
  void gradient_add( const Eigen::VectorXd & , 
							   const Eigen::VectorXd & ,
							   const Eigen::VectorXd & );
	void step_theta(const double);
  	double  dtau;
  	double ddtau;
    void initFromList(Rcpp::List const &);
    void initFromList(Rcpp::List const &, Rcpp::List const &);
    Rcpp::List output_list();
    void print_parameters();
    
        
    Eigen::VectorXd  get_gradient() { Eigen::VectorXd g(1); g[0] = dtau; return(g);};
    void  clear_gradient() {dtau = 0;};
    
};

class fd2Operator : public constMatrix {

public:
  Rcpp::List output_list();
  double trace_variance( const Eigen::SparseMatrix<double,0,int> & ) ;
};

class MaternOperator : public operatorMatrix{
  protected:
    double ldet;
    
    Eigen::VectorXd  h; 
    double h_average;
    Eigen::VectorXd g,p;
    Eigen::SparseMatrix<double,0,int> G, C;
    double kappa, dkappa, ddkappa, dtau, ddtau;
    bool use_chol;
    double counter;
    int calc_det;
    solver * Qepssolver;
    void set_matrices();
    SparseMatrix<double,0,int> d2tauQ, dtauQ, dkappaQ, d2kappaQ;
    double tau_trace, tau_trace2, kappa_trace, kappa_trace2;
  public:

  	double tau;
  	int n;
    MaternOperator(){ counter = 0;};
    void initFromList(Rcpp::List const &){Rcpp::Rcout << "Supply solver list when using initFromlist";};
    void initFromList(Rcpp::List const &, Rcpp::List const &);

    Rcpp::List output_list();
    void gradient( const Eigen::VectorXd &, const Eigen::VectorXd & );
    void gradient_init(const int, const int);
    void gradient_add( const Eigen::VectorXd & , 
							   const Eigen::VectorXd & ,
							   const Eigen::VectorXd & );
    void step_theta(const double );
    Eigen::VectorXd kappaVec;
    void print_parameters();
    double trace_variance( const Eigen::SparseMatrix<double,0,int> & );
    
    Eigen::VectorXd  get_gradient();
    void  clear_gradient();

};


#endif