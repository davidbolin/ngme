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
    solver ** Qsolver;
    int Q_act;
    Rcpp::List  out_list;
  public:
    int nop;
    std::vector<Eigen::VectorXd >  h;
  	Eigen::MatrixXd Cov_theta;// assymptotic covariance of the parameters
	  int npars; // number of parameters
    std::vector<Eigen::VectorXd >  loc; // location of the position
    std::vector<int > d; //dimension
    Eigen::SparseMatrix<double,0,int> *Q; // the generic matrix object
    //std::vector< Eigen::SparseMatrix<double,0,int> > Q;
    std::vector<Eigen::MatrixXd> K;                   // the generic matrix object if Q is full!



    double tau;
    Eigen::VectorXd  tauVec;
    int counter;

    operatorMatrix() {Qsolver = NULL;};
    virtual ~operatorMatrix() = default;

    virtual Eigen::VectorXd  get_gradient() { Eigen::VectorXd temp; return(temp);};
    virtual void  clear_gradient() {};
    virtual void initFromList(Rcpp::List const &)=0;
    virtual void initFromList(Rcpp::List const &, Rcpp::List const &) {Rcpp::Rcout << "initFromList(list1,list2) not implimented in operatorMatrix\n";};
    virtual Rcpp::List output_list() = 0;
    virtual void gradient( const Eigen::VectorXd &, const Eigen::VectorXd & ){};
    virtual void gradient_init( const int, const int){};
    virtual void gradient_add( const Eigen::VectorXd & ,
							   const Eigen::VectorXd & ,
							   const Eigen::VectorXd & ,
							   int,
							   const double){};
    virtual Eigen::MatrixXd d2Given( const Eigen::VectorXd & ,
                               const Eigen::VectorXd & ,
                               const Eigen::VectorXd & ,
                               int,
                               const double){return(Eigen::MatrixXd::Zero(0,0));};
    virtual void step_theta(const double stepsize,
    						const double learning_rate = 0,
    						const double polyak_rate   = -1,
    						const int burnin = 0){};
    virtual void print_parameters(){};
    virtual void get_param_names(Rcpp::StringVector & names){};
    virtual void get_param(std::vector<double> & ){};
    virtual double trace_variance( const Eigen::SparseMatrix<double,0,int> &, int ){return 1;};

	/*
    	stores the covariance of the parameters
    */
	void set_covariance(const Eigen::MatrixXd & Cov_in) {Cov_theta = Cov_in;};
};

class constMatrix : public operatorMatrix{
  protected:
    Eigen::SparseMatrix<double,0,int> m;
    Eigen::VectorXd v;
    std::vector<double>  m_loc;
    std::vector<double>  h_average;
  public:
  double term1, term2, term3;

  void get_param_names(Rcpp::StringVector & names){
    names.push_back("tau_operator");
  };
  void get_param(std::vector<double> & param_in){
    param_in.push_back(tau);
  };
	~constMatrix();
    Eigen::MatrixXd d2Given( const Eigen::VectorXd & ,
                               const Eigen::VectorXd & ,
                               const Eigen::VectorXd & ,
                               int,
                               const double);
	void gradient(const Eigen::VectorXd &, const Eigen::VectorXd & );
  void gradient_init(const int, const int);
  void gradient_add( const Eigen::VectorXd & ,
							   const Eigen::VectorXd & ,
							   const Eigen::VectorXd & ,
							   int,
							   const double);
	void step_theta(const double stepsize,
					const double learning_rate = 0,
					const double polyak_rate   = -1,
					const int burnin = 0);
  	double  dtau;
  	double ddtau;
	double dtau_old;
    void initFromList(Rcpp::List const &);
    void initFromList(Rcpp::List const &, Rcpp::List const &);
    Rcpp::List output_list();
    void print_parameters();

    Eigen::VectorXd  get_gradient() { Eigen::VectorXd g(1); g[0] = dtau; return(g);};
    void  clear_gradient() {dtau = 0;};
    double trace_variance( const Eigen::SparseMatrix<double,0,int> &, int  ) ;


};

class fd2Operator : public constMatrix {

public:
  //fd2Operator(){};
  //~fd2Operator();

  double trace_variance( const Eigen::SparseMatrix<double,0,int> &, int  ) ;
};

class MaternOperator : public operatorMatrix{
  protected:
    double ldet;
    int is_initialized = 0;
    std::vector<double>  h_average,tau_trace, tau_trace2, kappa_trace, kappa_trace2;
    Eigen::VectorXd g,p;
    Eigen::SparseMatrix<double,0,int> *G, *C;
    double kappa, dkappa, ddkappa, dtau, ddtau, ddtaukappa;
    double dtau_old, dkappa_old;
    bool use_chol;
    std::vector<int> matrix_set;
    double counter;
    int calc_det;
    solver ** Qepssolver;
    SparseMatrix<double,0,int> *d2tauQ, *dtauQ, *dkappaQ, *d2kappaQ;
    void set_matrices();
    void set_matrix(const int);
  public:


    void get_param(std::vector<double> & param_in){
      param_in.push_back(tau);
      param_in.push_back(kappa);
    };
    void get_param_names(Rcpp::StringVector & names){
      names.push_back("tau_operator");
      names.push_back("kappa_operator");
    };
  	double tau;

    MaternOperator(){ counter = 0;};
    ~MaternOperator();

    void initFromList(Rcpp::List const &){Rcpp::Rcout << "Supply solver list when using initFromlist";};
    void initFromList(Rcpp::List const &, Rcpp::List const &);

    Rcpp::List output_list();
    void gradient( const Eigen::VectorXd &, const Eigen::VectorXd & );
    void gradient_init(const int, const int);
    void gradient_add( const Eigen::VectorXd & ,
							   const Eigen::VectorXd & ,
							   const Eigen::VectorXd & ,
							   int,
							   const double);

    Eigen::MatrixXd d2Given( const Eigen::VectorXd & ,
                               const Eigen::VectorXd & ,
                               const Eigen::VectorXd & ,
                               int,
                               const double);
    void step_theta(const double stepsize,
    				const double learning_rate = 0,
    				const double polyak_rate   = -1,
    				const int burnin = 0);
    Eigen::VectorXd kappaVec;
    void print_parameters();
    double trace_variance( const Eigen::SparseMatrix<double,0,int> &, int );

    Eigen::VectorXd  get_gradient();
    void  clear_gradient();

};

class ExponentialOperator : public operatorMatrix{
protected:
  double ldet;
  int is_initialized = 0;
  std::vector<double>  h_average,tau_trace, tau_trace2, kappa_trace, kappa_trace2;
  Eigen::VectorXd g,p;
  Eigen::SparseMatrix<double,0,int> *G, *C;
  double kappa, dkappa, ddkappa, dtau, ddtau, ddtaukappa;
  double dtau_old, dkappa_old;
  bool use_chol;
  std::vector<int> matrix_set;
  double counter;
  int calc_det;
  solver ** Qepssolver;
  SparseMatrix<double,0,int> *d2tauQ, *dtauQ, *dkappaQ, *d2kappaQ;
  void set_matrices();
  void set_matrix(const int);
public:
  
  
  void get_param(std::vector<double> & param_in){
    param_in.push_back(tau);
    param_in.push_back(kappa);
  };
  void get_param_names(Rcpp::StringVector & names){
    names.push_back("tau_operator");
    names.push_back("kappa_operator");
  };
  double tau;
  
  ExponentialOperator(){ counter = 0;};
  ~ExponentialOperator();
  
  void initFromList(Rcpp::List const &){Rcpp::Rcout << "Supply solver list when using initFromlist";};
  void initFromList(Rcpp::List const &, Rcpp::List const &);
  
  Rcpp::List output_list();
  void gradient( const Eigen::VectorXd &, const Eigen::VectorXd & );
  void gradient_init(const int, const int);
  void gradient_add( const Eigen::VectorXd & ,
                     const Eigen::VectorXd & ,
                     const Eigen::VectorXd & ,
                     int,
                     const double);
  
  Eigen::MatrixXd d2Given( const Eigen::VectorXd & ,
                           const Eigen::VectorXd & ,
                           const Eigen::VectorXd & ,
                           int,
                           const double);
  void step_theta(const double stepsize,
                  const double learning_rate = 0,
                  const double polyak_rate   = -1,
                  const int burnin = 0);
  Eigen::VectorXd kappaVec;
  void print_parameters();
  double trace_variance( const Eigen::SparseMatrix<double,0,int> &, int );
  
  Eigen::VectorXd  get_gradient();
  void  clear_gradient();
  
};

class MaternOperator2D : public operatorMatrix{
protected:
  double ldet;
  int is_initialized = 0;
  std::vector<double>  h_average,tau1_trace, tau1_trace2, kappa1_trace, kappa1_trace2;
  std::vector<double>  tau2_trace, tau2_trace2, kappa2_trace, kappa2_trace2;
  std::vector<double>  rho_trace, rho_trace2, theta_trace, theta_trace2;
  Eigen::VectorXd g,p;
  Eigen::SparseMatrix<double,0,int> *G, *C;
  double kappa1, dkappa1, ddkappa1, dtau1, ddtau1;
  double kappa2, dkappa2, ddkappa2, dtau2, ddtau2;
  double rho, drho, ddrho, theta, dtheta, ddtheta;
  double dtau1_old, dkappa1_old, dtau2_old, dkappa2_old,drho_old, dtheta_old;;
  bool use_chol, estimate_theta;
  std::vector<int> matrix_set;
  double counter;
  int calc_det;
  solver ** Qepssolver;
  SparseMatrix<double,0,int> *d2tau1Q, *dtau1Q, *dkappa1Q, *d2kappa1Q;
  SparseMatrix<double,0,int> *d2tau2Q, *dtau2Q, *dkappa2Q, *d2kappa2Q;
  SparseMatrix<double,0,int> *d2rhoQ, *drhoQ, *dthetaQ, *d2thetaQ;
  void set_matrices();
  void set_matrix(const int);
public:
  double tau1, tau2;
  Eigen::VectorXd kappa1Vec, kappa2Vec, tau1Vec, tau2Vec, rhoVec, thetaVec;
  
  void get_param(std::vector<double> & param_in){
    param_in.push_back(tau1);
    param_in.push_back(tau2);
    param_in.push_back(kappa1);
    param_in.push_back(kappa2);
    param_in.push_back(rho);
    param_in.push_back(theta);
  };
  void get_param_names(Rcpp::StringVector & names){
    names.push_back("tau1_operator");
    names.push_back("tau2_operator");
    names.push_back("kappa1_operator");
    names.push_back("kappa2_operator");
    names.push_back("rho_operator");
    names.push_back("theta_operator");
  };
  
  MaternOperator2D(){ counter = 0;};
  ~MaternOperator2D();
  
  void initFromList(Rcpp::List const &){Rcpp::Rcout << "Supply solver list when using initFromlist";};
  void initFromList(Rcpp::List const &, Rcpp::List const &);
  
  Rcpp::List output_list();
  void gradient( const Eigen::VectorXd &, const Eigen::VectorXd & );
  void gradient_init(const int, const int);
  void gradient_add( const Eigen::VectorXd & ,
                     const Eigen::VectorXd & ,
                     const Eigen::VectorXd & ,
                     int,
                     const double);
  
  Eigen::MatrixXd d2Given( const Eigen::VectorXd & ,
                           const Eigen::VectorXd & ,
                           const Eigen::VectorXd & ,
                           int,
                           const double);
  void step_theta(const double stepsize,
                  const double learning_rate = 0,
                  const double polyak_rate   = -1,
                  const int burnin = 0);
  
  void print_parameters();
  double trace_variance( const Eigen::SparseMatrix<double,0,int> &, int );
  
  Eigen::VectorXd  get_gradient();
  void  clear_gradient();
  
};

#endif
