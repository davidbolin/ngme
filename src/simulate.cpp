
#include <Rcpp.h>
#include <random>
#include <chrono>
#include <string>
#include "operatorMatrix.h"
#include "solver.h"
#include "operator_helper.h"
#include "GIG.h"
#include "MixedEffect.h"
#include "measError.h"
using namespace Rcpp;



/*
	Simulating from the prior model
*/
// [[Rcpp::export]]
List simulateLongGH_cpp(Rcpp::List in_list)
{

 	//**********************************
	//setting up the main data
	//**********************************
	Rcpp::List obs_list  = Rcpp::as<Rcpp::List> (in_list["obs_list"]);
	int nindv = obs_list.length(); //count number of patients
	std::vector< Eigen::SparseMatrix<double,0,int> > As( nindv);
	std::vector< Eigen::VectorXd > Ys( nindv);
	int counter = 0;
	for( List::iterator it = obs_list.begin(); it != obs_list.end(); ++it ) {
    	List obs_tmp = Rcpp::as<Rcpp::List>(*it);
    	As[counter] = Rcpp::as<Eigen::SparseMatrix<double,0,int> >(obs_tmp["A"]);
    	Ys[counter++] = Rcpp::as<Eigen::VectorXd>(obs_tmp["Y"]);
    }


	//**********************************
	//operator setup
	//***********************************
	Rcpp::List operator_list  = Rcpp::as<Rcpp::List> (in_list["operator_list"]);
	operator_list["nIter"] = 1;
	std::string type_operator = Rcpp::as<std::string>(operator_list["type"]);
	operatorMatrix* Kobj;
	operator_select(type_operator, &Kobj);
	Kobj->initFromList(operator_list, List::create(Rcpp::Named("use.chol") = 1));
	Eigen::VectorXd h = Rcpp::as<Eigen::VectorXd>( operator_list["h"]);
	//Prior solver
	cholesky_solver Qsolver;
	Qsolver.init( Kobj->d, 0, 0, 0);
	Qsolver.analyze( Kobj->Q);
	Qsolver.compute( Kobj->Q);


	//Create solvers for each patient
	std::vector<  cholesky_solver >  Solver( nindv);
	Eigen::SparseMatrix<double, 0, int> Q;

  counter = 0;
	for( List::iterator it = obs_list.begin(); it != obs_list.end(); ++it ) {
    	List obs_tmp = Rcpp::as<Rcpp::List>( *it);
    	Solver[counter].init(Kobj->d, 0, 0, 0);
    	Q = Eigen::SparseMatrix<double,0,int>(Kobj->Q.transpose());
    	Q = Q  * Kobj->Q;
    	Q = Q + As[counter].transpose()*As[counter];
    	Solver[counter].analyze(Q);
    	Solver[counter++].compute(Q);
  }

	//**********************************
	// mixed effect setup
	//***********************************
	Rcpp::List mixedEffect_list  = Rcpp::as<Rcpp::List> (in_list["mixedEffect_list"]);
	std::string type_mixedEffect = Rcpp::as<std::string> (mixedEffect_list["noise"]);
	MixedEffect *mixobj;
	if(type_mixedEffect == "Normal")
    	mixobj = new NormalMixedEffect;
	else
		mixobj   = new NIGMixedEffect;

	mixobj->initFromList(mixedEffect_list);


  //***********************************
	// measurement error setup
	//***********************************
  Rcpp::List MeasureError_list  = Rcpp::as<Rcpp::List> (in_list["measurment_list"]);
  MeasurementError *errObj;
  std::string MeasureNoise = Rcpp::as <std::string> (MeasureError_list["noise"]);
  if(MeasureNoise == "Normal")
    errObj = new GaussianMeasurementError;
  else
    errObj = new NIGMeasurementError;

  errObj->initFromList(MeasureError_list);


	//***********************************
	// stochastic processes setup
	//***********************************
	Rcpp::List processes_list   = Rcpp::as<Rcpp::List>  (in_list["processes_list"]);
	Rcpp::List V_list           = Rcpp::as<Rcpp::List>  (processes_list["V"]);
	std::string type_processes  = Rcpp::as<std::string> (processes_list["noise"]);
  double nu = -1;
  double mu = 0;
  if(type_processes != "Normal"){
    nu  = Rcpp::as< double > (processes_list["nu"]);
    mu  = Rcpp::as< double > (processes_list["mu"]);
  }
	std::vector< Eigen::VectorXd >   Vs( nindv);
  std::vector< Eigen::VectorXd > Xs( nindv);
  std::vector< Eigen::VectorXd > Zs( nindv);
  for(int i = 0; i < nindv; i++ ){
    	Xs[i].resize( Kobj->d );
    	Vs[i] = h;
    	Zs[i].resize(Kobj->d);
  }
  	/*
  	Simulation objects
  	*/
  	std::mt19937 random_engine;
  	std::normal_distribution<double> normal;
  	std::default_random_engine gammagenerator;
  	random_engine.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());
  	gig rgig;
  	rgig.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());
  	Eigen::VectorXd  z;
  	z.setZero(Kobj->d);

  	Eigen::VectorXd b;
  	b.setZero(Kobj->d);

    //*********************************************
    //        simulating the measurement error
    //*********************************************
    std::vector< Eigen::VectorXd > Ysim = errObj->simulate( Ys);

    //*********************************************
    //        simulating the mixed effect
    //*********************************************
    mixobj->simulate();
    for(int i = 0; i < Ysim.size(); i++)
    {
 		mixobj->add_inter(i, Ysim[i]);
		mixobj->add_cov(i, Ysim[i]);
 	  }
    //*********************************************
    //        simulating the processes
    //*********************************************
    Eigen::VectorXd iV;
    for(int i = 0; i < Ysim.size(); i++) {

      if(type_processes != "Normal")
        Vs[i] = sampleV_pre(rgig, h, nu, type_processes );

      iV.resize(Vs[i].size());
      iV.array() = Vs[i].array().inverse();

      for(int ii =0; ii < Kobj->d; ii++){
        z[ii] =   sqrt(Vs[i][ii]) * normal(random_engine);
        if(type_processes != "Normal")
          z[ii] += - mu * h[ii] + Vs[i][ii] * mu;
        z[ii];
      }

      Eigen::SparseMatrix<double,0,int> K = Eigen::SparseMatrix<double,0,int>(Kobj->Q);
      //	K *= sqrt(tau);

      Eigen::SparseLU< Eigen::SparseMatrix<double,0,int> > chol(K);  // performs a Cholesky factorization of A


      Xs[i] = chol.solve(z);         // use the factorization to solve for the given right hand side
      Zs[i] = z;
      Ysim[i] += As[i] * Xs[i];
  }


  Rcpp::List out_list;
  out_list["Y"]    = Ysim;
  out_list["U"]    = mixobj->U;
  out_list["X"]    = Xs;
  out_list["Z"]    = Zs;
  Eigen::SparseMatrix<double,0,int> K = Eigen::SparseMatrix<double,0,int>(Kobj->Q);
      	//K *= sqrt(tau);
  out_list["K"]    = K;
  if(type_processes != "Normal")
    out_list["V"] = Vs;

  if(errObj->noise != "Normal")
    out_list["V_noise"] = errObj->Vs;
  return(out_list);

}

