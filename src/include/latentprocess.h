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

  int useEV; // for GAL use EV rather EV when computing gradient
  int npars;
  int store_param;
  int nindv; // number indiviuals
  std::vector< Eigen::VectorXd > ElogV_post;
  std::vector< Eigen::VectorXd > EV_post;
  std::vector< Eigen::VectorXd > EV;
  std::vector< Eigen::VectorXd > Xs;
  std::vector< Eigen::VectorXd > Ws;
  std::vector< Eigen::VectorXd >  Vs;
  std::vector< Eigen::VectorXd >  mu0;
  Eigen::SparseMatrix<double,0,int>  Q;
  std::vector < Eigen::VectorXd > h;
  Eigen::VectorXd  iV;
  std::string type_process;
  Process() {};
  Eigen::MatrixXd Cov_theta;// assymptotic covariance of the parameters

  virtual Eigen::VectorXd get_mean_prior(const int i, const Eigen::SparseMatrix<double,0,int> & K){
  Eigen::VectorXd res;
  res.setZero(h[i].size());
  return(res);
  };
  virtual Eigen::VectorXd  get_gradient() { Eigen::VectorXd temp; return(temp);};
  virtual void  clear_gradient() {};
  virtual void get_param(std::vector<double> & param_in){};
  virtual void get_param_names(Rcpp::StringVector & names){
  };

  //print iteration data
  virtual void printIter(){};
  // setups to store the tracjetory
  virtual void setupStoreTracj(const int Niter){};
  virtual ~Process(){};
  virtual Rcpp::List toList() {};
  virtual Eigen::VectorXd  mean_X(const int i) {Eigen::VectorXd temp(h[i].size()); return temp.setZero(h[i].size());};
  virtual Eigen::VectorXd d2Given_cross(  const int i ,
                                          const Eigen::SparseMatrix<double,0,int> & K,
                                          const Eigen::SparseMatrix<double,0,int> & A,
                                          const Eigen::VectorXd& res,
                                          const double sigma,
                                          const Eigen::MatrixXd Bf,
                                          const Eigen::MatrixXd Br,
                                          const double weight){return(Eigen::VectorXd::Zero(0));};
  virtual  Eigen::VectorXd d2Given_v2_cross(  const int i ,
                                              const Eigen::SparseMatrix<double,0,int> & K,
                                              const Eigen::SparseMatrix<double,0,int> & A,
                                              const Eigen::VectorXd& res,
                                              const double sigma,
                                              const Eigen::MatrixXd Bf,
                                              const Eigen::MatrixXd Br,
                                              const Eigen::VectorXd& iV_noise,
                                              const double weight){return(Eigen::VectorXd::Zero(0));};;
  virtual Eigen::MatrixXd d2Given(  const int i ,
                                    const Eigen::SparseMatrix<double,0,int> & K,
                                    const Eigen::SparseMatrix<double,0,int> & A,
                                    const Eigen::VectorXd& res,
                                    const double sigma,
                                    const double trace_var,
                                    const double weight){return(Eigen::MatrixXd::Zero(0,0));};

  virtual void gradient( const int i ,
                         const Eigen::SparseMatrix<double,0,int> & K,
                         const Eigen::SparseMatrix<double,0,int> & A,
                         const Eigen::VectorXd& res,
                         const double sigma,
                         const double trace_var,
                         const double weight){};
  virtual void step_theta(const double stepsize,
                          const double learning_rate = 0,
                          const double polyak_rate   =  -1,
                          const int burnin = 0){};
  virtual void sample_V(const int,
                        gig &,
                        const Eigen::SparseMatrix<double,0,int> &){};


  virtual void simulate_V(const int,
                          gig &){};
  virtual void initFromList(const Rcpp::List &, const std::vector <Eigen::VectorXd > &){};
  // sampling where the measurement noise is normal
  virtual void sample_X(const int i,
                        Eigen::VectorXd & Z,
                        const Eigen::VectorXd & Y,
                        const Eigen::SparseMatrix<double,0,int> & Q,
                        const Eigen::SparseMatrix<double,0,int> & K,
                        const Eigen::SparseMatrix<double,0,int> & A,
                        const double sigma,
                        solver       & solver){};
  // sampling where the measurement noise is non-normal
  virtual void sample_Xv2(  const int i,
                            Eigen::VectorXd & Z,
                            const Eigen::VectorXd & Y,
                            const Eigen::SparseMatrix<double,0,int> & Q,
                            const Eigen::SparseMatrix<double,0,int> & K,
                            const Eigen::SparseMatrix<double,0,int> & A,
                            const double sigma,
                            solver       & solver,
                            const Eigen::VectorXd & iV_noise ){};


  virtual Eigen::MatrixXd d2Given_v2( const int i ,
                                      const Eigen::SparseMatrix<double,0,int> & K,
                                      const Eigen::SparseMatrix<double,0,int> & A,
                                      const Eigen::VectorXd& res,
                                      const double sigma,
                                      const Eigen::VectorXd& iV_noise,
                                      const double EiV_noise,
                                      const double trace_var,
                                      const double weight){return(Eigen::MatrixXd::Zero(0,0));};
  virtual void gradient_v2( const int i ,
                            const Eigen::SparseMatrix<double,0,int> & K,
                            const Eigen::SparseMatrix<double,0,int> & A,
                            const Eigen::VectorXd& res,
                            const double sigma,
                            const Eigen::VectorXd& iV_noise,
                            const double EiV_noise,
                            const double trace_var,
                            const double weight) {};

  //simulate from prior distribution
  virtual  void simulate(const int ,
                         Eigen::VectorXd & ,
                         const Eigen::SparseMatrix<double,0,int> & ,
                         const Eigen::SparseMatrix<double,0,int> & ,
                         Eigen::VectorXd& ,
                         solver  & ) = 0;


  /*
  stores the covariance of the parameters
  */
  void set_covariance(const Eigen::MatrixXd & Cov_in) {Cov_theta = Cov_in;};
};

class GaussianProcess : public Process{
  ~GaussianProcess();
  void initFromList(const Rcpp::List  &, const std::vector <Eigen::VectorXd > &);
  void sample_X( const int i,
                 Eigen::VectorXd & Z,
                 const Eigen::VectorXd & Y,
                 const Eigen::SparseMatrix<double,0,int> & Q,
                 const Eigen::SparseMatrix<double,0,int> & K,
                 const Eigen::SparseMatrix<double,0,int> & A,
                 const double sigma,
                 solver       & solver);
  void sample_Xv2(  const int i,
                    Eigen::VectorXd & Z,
                    const Eigen::VectorXd & Y,
                    const Eigen::SparseMatrix<double,0,int> & Q,
                    const Eigen::SparseMatrix<double,0,int> & K,
                    const Eigen::SparseMatrix<double,0,int> & A,
                    const double sigma,
                    solver       & solver,
                    const Eigen::VectorXd & iV_noise);

  Rcpp::List toList();


  //simulate from prior distribution
  void  simulate(const int ,
                 Eigen::VectorXd & ,
                 const Eigen::SparseMatrix<double,0,int> & ,
                 const Eigen::SparseMatrix<double,0,int> & ,
                 Eigen::VectorXd&,
                 solver  &  );
};





class GHProcessBase : public Process{


public:
  std::vector< Eigen::VectorXd > sampleVbool;
  std::vector< Eigen::VectorXd > h2;
  std::vector< double > h_sum;
  std::vector< double >  h_min;
  double h_MIN;
  std::vector< double >  h_digamma;
  std::vector< double >  h_trigamma;
  std::vector< double >  h3_mean;
  double dmu ;
  double dnu, ddnu ;
  double dnu_prev;
  double dmu_prev;
  std::vector< Eigen::VectorXd > EiV;
  Eigen::MatrixXd mu_vec;
  Eigen::MatrixXd nu_vec;
  int vec_counter;
  double counter;
  double  ddmu_1, ddmu_2;
  std::vector< double > H_mu;
  std::vector<double> Vv_mean;
  std::vector< Eigen::VectorXd > toSampleV;

  Eigen::VectorXd get_mean_prior(const int, const Eigen::SparseMatrix<double,0,int> & );

  virtual void update_nu();
  virtual void grad_nu(const int, const double);
  virtual void gradient_mu_centered(const int ,
                            const Eigen::SparseMatrix<double,0,int> & ,
                            const double);


  double term1,term2;

  double mu;
  double nu;
  virtual void initFromList(const Rcpp::List  &, const std::vector <Eigen::VectorXd > &) = 0;
  void sample_X( const int i,
                 Eigen::VectorXd & Z,
                 const Eigen::VectorXd & Y,
                 const Eigen::SparseMatrix<double,0,int> & Q,
                 const Eigen::SparseMatrix<double,0,int> & K,
                 const Eigen::SparseMatrix<double,0,int> & A,
                 const double sigma,
                 solver       & solver);

  void get_param(std::vector<double> & param_in){
    if( type_process != "CH"){
      param_in.push_back(mu);
      param_in.push_back(nu);
    }else{
      param_in.push_back(0);
    }
  };
  void get_param_names(Rcpp::StringVector & names){

    names.push_back("mu_process");
    if( type_process != "CH"){
      names.push_back("nu_process");
    }
  };

  void sample_Xv2(  const int i,
                    Eigen::VectorXd & Z,
                    const Eigen::VectorXd & Y,
                    const Eigen::SparseMatrix<double,0,int> & Q,
                    const Eigen::SparseMatrix<double,0,int> & K,
                    const Eigen::SparseMatrix<double,0,int> & A,
                    const double sigma,
                    solver       & solver,
                    const Eigen::VectorXd & iV_noise);

  void gradient( const int i ,
                 const Eigen::SparseMatrix<double,0,int> & K,
                 const Eigen::SparseMatrix<double,0,int> & A,
                 const Eigen::VectorXd& res,
                 const double sigma,
                 const double trace_var,
                 const double weight);

  Eigen::MatrixXd d2Given(  const int i ,
                            const Eigen::SparseMatrix<double,0,int> & K,
                            const Eigen::SparseMatrix<double,0,int> & A,
                            const Eigen::VectorXd& res,
                            const double sigma,
                            const double trace_var,
                            const double weight);

  Eigen::VectorXd d2Given_cross(  const int i ,
                                  const Eigen::SparseMatrix<double,0,int> & K,
                                  const Eigen::SparseMatrix<double,0,int> & A,
                                  const Eigen::VectorXd& res,
                                  const double sigma,
                                  Eigen::MatrixXd Bf,
                                  Eigen::MatrixXd Br,
                                  const double weight);
  Eigen::VectorXd d2Given_v2_cross(  const int i ,
                                     const Eigen::SparseMatrix<double,0,int> & K,
                                     const Eigen::SparseMatrix<double,0,int> & A,
                                     const Eigen::VectorXd& res,
                                     const double sigma,
                                     Eigen::MatrixXd Bf,
                                     Eigen::MatrixXd Br,
                                     const Eigen::VectorXd& iV_noise,
                                     const double weight);

  Eigen::MatrixXd d2Given_v2( const int i ,
                              const Eigen::SparseMatrix<double,0,int> & K,
                              const Eigen::SparseMatrix<double,0,int> & A,
                              const Eigen::VectorXd& res,
                              const double sigma,
                              const Eigen::VectorXd& iV_noise,
                              const double EiV_noise,
                              const double trace_var,
                              const double weight);

  void gradient_v2( const int i ,
                    const Eigen::SparseMatrix<double,0,int> & K,
                    const Eigen::SparseMatrix<double,0,int> & A,
                    const Eigen::VectorXd& res,
                    const double sigma,
                    const Eigen::VectorXd& iV_noise,
                    const double EiV_noise,
                    const double trace_var,
                    const double weight);

  void step_theta(const double stepsize, const double learning_rate = 0, const double polyak_rate = -1, const int burnin = 0);
  virtual void step_mu(const double, const double, const int, const double);
  virtual void step_nu(const double, const double, const int,const double);
  void printIter();
  virtual Rcpp::List toList();
  virtual void setupStoreTracj(const int);
  virtual void sample_V(const int,
                gig &,
                const Eigen::SparseMatrix<double,0,int> &);

  //simulate from prior distribution
  void simulate(const int ,
                Eigen::VectorXd & ,
                const Eigen::SparseMatrix<double,0,int> & ,
                const Eigen::SparseMatrix<double,0,int> & ,
                Eigen::VectorXd& ,
                solver  &  );
  void simulate_V(const int,
                  gig &);

  virtual Eigen::VectorXd  get_gradient();
  virtual void  clear_gradient();
  virtual Eigen::VectorXd  mean_X(const int );
};



class GHProcess : public GHProcessBase{



public:
 
  void initFromList(const Rcpp::List  &, const std::vector <Eigen::VectorXd > &);
};

class MGHProcess : public GHProcessBase{



public:
  void gradient_mu_centered(const int ,
                            const Eigen::SparseMatrix<double,0,int> & ,
                            const double);
  void step_mu(const double, const double, const int, const double);
  void grad_nu(const int i, const double weight);
  void step_nu(const double, const double, const int,const double);
  void update_nu();
  void clear_gradient();
  Eigen::VectorXd  Nu;
  Eigen::VectorXd  Mu;
  Eigen::VectorXd  dNu;
  Eigen::VectorXd  dMu;
  Eigen::VectorXd  dMu_prev;
  Eigen::VectorXd  dNu_prev;
  Eigen::MatrixXd  ddNu;
  Eigen::MatrixXd  ddMu;


  std::vector< Eigen::MatrixXd>  Bmu;
  std::vector< Eigen::MatrixXd>  Bnu;
  void setupStoreTracj(const int);
  void initFromList(const Rcpp::List  &, const std::vector <Eigen::VectorXd > &);
  Rcpp::List toList();
  Eigen::VectorXd  get_gradient();
  virtual Eigen::VectorXd  mean_X(const int );
  void sample_V(const int ,
                gig & ,
                const Eigen::SparseMatrix<double,0,int> & );
};


#endif