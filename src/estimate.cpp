#include <Rcpp.h>
#include <RcppEigen.h>
#include <random>
#include <chrono>
#include <vector>
#include <string>
#include <math.h>
#include "operatorMatrix.h"
#include "operator_helper.h"
#include "solver.h"
#include "GIG.h"
#include "MixedEffect.h"
#include "measError.h"
#include "process.h"
using namespace Rcpp;

double estDigamma(double x)
{
  return(R::digamma(x));
}

// [[Rcpp::export]]
List estimateLong_cpp(Rcpp::List in_list)
{
	//**********************************
	//      basic parameter
	//**********************************

	double pSubsample = Rcpp::as< double > (in_list["pSubsample"]);
	int nIter      = Rcpp::as< double > (in_list["nIter"]);
	int nSim       = Rcpp::as< double > (in_list["nSim"]);
  	int nBurnin    = Rcpp::as< double > (in_list["nBurnin"] );
  	int silent     = Rcpp::as< int    > (in_list["silent"]);
  	double alpha     = Rcpp::as< double    > (in_list["alpha"]);
  	double step0     = Rcpp::as< double    > (in_list["step0"]);
	//**********************************
	//     setting up the main data
	//**********************************
	Rcpp::List obs_list  = Rcpp::as<Rcpp::List> (in_list["obs_list"]);
	int nindv = obs_list.length(); //count number of patients
  int nSubsample = ceil(pSubsample * nindv);
	std::vector< Eigen::SparseMatrix<double,0,int> > As( nindv);
	std::vector< Eigen::VectorXd > Ys( nindv);
	int count;
	count = 0;
	for( List::iterator it = obs_list.begin(); it != obs_list.end(); ++it ) {
    	List obs_tmp = Rcpp::as<Rcpp::List>(*it);
    	As[count] = Rcpp::as<Eigen::SparseMatrix<double,0,int> >(obs_tmp["A"]);
    	Ys[count] = Rcpp::as<Eigen::VectorXd>(obs_tmp["Y"]);
      count++;
    }

	//**********************************
	//operator setup
	//***********************************
	Rcpp::List operator_list  = Rcpp::as<Rcpp::List> (in_list["operator_list"]);
	operator_list["nIter"] = nIter;
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

	count = 0;
	for( List::iterator it = obs_list.begin(); it != obs_list.end(); ++it ) {
    List obs_tmp = Rcpp::as<Rcpp::List>( *it);
    Solver[count].init(Kobj->d, 0, 0, 0);
  	Q = Eigen::SparseMatrix<double,0,int>(Kobj->Q.transpose());
    Q = Q  * Kobj->Q;
  	Q = Q + As[count].transpose()*As[count];
    Solver[count].analyze(Q);
  	Solver[count].compute(Q);
    count++;
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
	mixobj->setupStoreTracj(nIter);

  //**********************************
	// measurement setup
	//***********************************
  MeasurementError *errObj;
  Rcpp::List measurementError_list  = Rcpp::as<Rcpp::List> (in_list["measurementError_list"]);
  std::string type_MeasurementError= Rcpp::as <std::string> (measurementError_list["noise"]);
  if(type_MeasurementError == "Normal")
    errObj = new GaussianMeasurementError;
  else
    errObj = new NIGMeasurementError;

  errObj->initFromList(measurementError_list);
  errObj->setupStoreTracj(nIter);

	//**********************************
	// stochastic processes setup
	//***********************************
	Rcpp::List processes_list   = Rcpp::as<Rcpp::List>  (in_list["processes_list"]);
	Rcpp::List V_list           = Rcpp::as<Rcpp::List>  (processes_list["V"]);
	std::string type_processes  = Rcpp::as<std::string> (processes_list["noise"]);

  Process *process;

  if (type_processes != "Normal"){
  	process  = new GHProcess;
  }else{ process  = new GaussianProcess;}

  process->initFromList(processes_list, h);
  process->setupStoreTracj(nIter);
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

  Eigen::VectorXd b, Ysim;
	b.setZero(Kobj->d);

  std::vector<int> longInd;
  for (int i=0; i< nindv; i++) longInd.push_back(i);

  for(int iter=0; iter < nIter + nBurnin; iter++){
    /*
      printing output
    */
    if(silent == 0){
      Rcpp::Rcout << "i = " << iter << ": \n";
      process->printIter();
      Rcpp::Rcout << "\n";
      Kobj->print_parameters();
      mixobj->printIter();
      errObj->printIter();
      Rcpp::Rcout << "\n";
    }

    Eigen::SparseMatrix<double,0,int> K = Eigen::SparseMatrix<double,0,int>(Kobj->Q);

    // subsampling
    //sample Nlong values without replacement from 1:nrep
    std::random_shuffle(longInd.begin(), longInd.end());
    Kobj->gradient_init(nSubsample,nSim);
    for(int ilong = 0; ilong < nSubsample; ilong++ )
    {
      int i = longInd[ilong];
      Eigen::SparseMatrix<double,0,int> A = As[i];
      Eigen::VectorXd  Y = Ys[i];
      for(int ii = 0; ii < nSim; ii ++)
      {
      	Eigen::VectorXd  res = Y;
      	//***************************************
      	//***************************************
      	//   building the residuals and sampling
      	//***************************************
      	//***************************************

      	// removing fixed effect from Y
      	mixobj->remove_cov(i, res);

    		res -= A * process->Xs[i];
  			//***********************************
      	// mixobj sampling
    		//***********************************
      	if(type_MeasurementError == "Normal")
    			mixobj->sampleU( i, res, 2 * log(errObj->sigma));
  			else
  			  mixobj->sampleU2( i, res, errObj->Vs[i].cwiseInverse(), 2 * log(errObj->sigma));

        mixobj->remove_inter(i, res);

      	//***********************************
    		// sampling processes
  			//***********************************
      	Eigen::SparseMatrix<double,0,int> Q = Eigen::SparseMatrix<double,0,int>(K.transpose());
    		Eigen::VectorXd iV(process->Vs[i].size());
  			iV.array() = process->Vs[i].array().inverse();
      	Q =  Q * iV.asDiagonal();
    		Q =  Q * K;
        for(int j =0; j < Kobj->d; j++)
    			z[j] =  normal(random_engine);

      	res += A * process->Xs[i];
      	//Sample X|Y, V, sigma

        if(type_MeasurementError == "Normal")
      		process->sample_X(i, z, res, Q, K, A, errObj->sigma, Solver[i]);
        else
          process->sample_Xv2( i, z, res, Q, K, A, errObj->sigma, Solver[i], errObj->Vs[i].cwiseInverse());

        res -= A * process->Xs[i];

        if(res.cwiseAbs().sum() > 1e16){
          throw("res outof bound\n");
        }

        // sample V| X
        process->sample_V(i, rgig, K);

        // random variance noise sampling
     		if(type_MeasurementError != "Normal"){
      	  errObj->sampleV(i, res);
     		}

      	//***************************************
      	//  computing gradients
      	//***************************************
      	if(iter >= nBurnin){

      		// mixobj gradient
      		mixobj->add_inter(i, res);
      	  if(type_MeasurementError != "Normal")
      		  mixobj->gradient2(i, res, errObj->Vs[i].cwiseInverse(), 2 * log(errObj->sigma), errObj->EiV);
          else
            mixobj->gradient(i, res, 2 * log(errObj->sigma));

      		mixobj->remove_inter(i, res);


      	  // measurent error  gradient
  			  errObj->gradient(i, res);

      	  // operator gradient
      		Kobj->gradient_add( process->Xs[i], process->Vs[i].cwiseInverse());

      		// process gradient
          if(type_MeasurementError != "Normal"){
              
              process->gradient_v2(i,
                                K,
                                A,
                                res,
                                errObj->sigma,
                                errObj->Vs[i].cwiseInverse(),
                                errObj->EiV,
                                Kobj->trace_variance(A));
            }else{
              process->gradient(i,
                                K,
                                A,
                                res,
                                errObj->sigma,
                                Kobj->trace_variance(A));
            }
      		}
      }
    }
    //**********************************
  	//  gradient step
	//***********************************
    if(iter >= nBurnin){
      double stepsize = step0 / pow(iter - nBurnin + 1, alpha);
      mixobj->step_theta(stepsize);
      errObj->step_theta(stepsize);
      Kobj->step_theta(stepsize);
      process->step_theta(stepsize);
    }
  }
  if(silent == 0)
  	Rcpp::Rcout << "Done, storing results\n";
  // storing the results
  Rcpp::List out_list;
  out_list["pSubsample"]       = pSubsample;
  out_list["nIter"]            = nIter;
  out_list["nSim"]             = nSim;
  out_list["nBurnin"]          = nBurnin;
  out_list["silent"]           = silent;
  out_list["step0"]            = step0;
  out_list["alpha"]            = alpha;
  out_list["obs_list"]         = obs_list;
  out_list["Xs"]               = process->Xs;
  out_list["Vs"]               = process->Vs;

  Rcpp::List mixobj_list       = mixobj->toList();
  out_list["mixedEffect_list"] = mixobj_list;


  Rcpp::List errobj_list            = errObj->toList();
  out_list["measurementError_list"] = errobj_list;

  Rcpp::List process_list           = process->toList();
  out_list["processes_list"]        = process_list;

  Rcpp::List olist          = Kobj->output_list();
  out_list["operator_list"] = olist;
    return(out_list);
}

// [[Rcpp::export]]
List estimateFisher(Rcpp::List in_list)
{
	//**********************************
	//      basic parameter
	//**********************************

	double pSubsample = Rcpp::as< double > (in_list["pSubsample"]);
	int nIter      = Rcpp::as< double > (in_list["nIter"]);
	int nSim       = Rcpp::as< double > (in_list["nSim"]);
  	int nBurnin    = Rcpp::as< double > (in_list["nBurnin"] );
  	int silent     = Rcpp::as< int    > (in_list["silent"]);
	//**********************************
	//     setting up the main data
	//**********************************
	Rcpp::List obs_list  = Rcpp::as<Rcpp::List> (in_list["obs_list"]);
	int nindv = obs_list.length(); //count number of patients
  	int nSubsample = ceil(pSubsample * nindv);
	std::vector< Eigen::SparseMatrix<double,0,int> > As( nindv);
	std::vector< Eigen::VectorXd > Ys( nindv);
	int count;
	count = 0;
	for( List::iterator it = obs_list.begin(); it != obs_list.end(); ++it ) {
    	List obs_tmp = Rcpp::as<Rcpp::List>(*it);
    	As[count] = Rcpp::as<Eigen::SparseMatrix<double,0,int> >(obs_tmp["A"]);
    	Ys[count] = Rcpp::as<Eigen::VectorXd>(obs_tmp["Y"]);
      count++;
    }

	//**********************************
	//operator setup
	//***********************************
	Rcpp::List operator_list  = Rcpp::as<Rcpp::List> (in_list["operator_list"]);
	operator_list["nIter"] = nIter;
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

	count = 0;
	for( List::iterator it = obs_list.begin(); it != obs_list.end(); ++it ) {
    	List obs_tmp = Rcpp::as<Rcpp::List>( *it);
    	Solver[count].init(Kobj->d, 0, 0, 0);
  		Q = Eigen::SparseMatrix<double,0,int>(Kobj->Q.transpose());
    	Q = Q  * Kobj->Q;
  		Q = Q + As[count].transpose()*As[count];
    	Solver[count].analyze(Q);
  		Solver[count].compute(Q);
    	count++;
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
  //**********************************
	// measurement setup
	//***********************************
  MeasurementError *errObj;
  Rcpp::List measurementError_list  = Rcpp::as<Rcpp::List> (in_list["measurementError_list"]);
  std::string type_MeasurementError= Rcpp::as <std::string> (measurementError_list["noise"]);
  if(type_MeasurementError == "Normal")
    errObj = new GaussianMeasurementError;
  else
    errObj = new NIGMeasurementError;

  errObj->initFromList(measurementError_list);
	//**********************************
	// stochastic processes setup
	//***********************************
	Rcpp::List processes_list   = Rcpp::as<Rcpp::List>  (in_list["processes_list"]);
	Rcpp::List V_list           = Rcpp::as<Rcpp::List>  (processes_list["V"]);
	std::string type_processes  = Rcpp::as<std::string> (processes_list["noise"]);
  Process *process;

  if (type_processes != "Normal"){
  	process  = new GHProcess;
  }else{ process  = new GaussianProcess;}

  process->initFromList(processes_list, h);
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

  Eigen::VectorXd b, Ysim;
	b.setZero(Kobj->d);

  std::vector<int> longInd;
  
  for (int i=0; i< nindv; i++) longInd.push_back(i);


	int npars =   mixobj->npars 
    			+ errObj->npars  
    			+ process->npars 
    			+ Kobj->npars;
  Eigen::MatrixXd Egrad2;
  Egrad2.setZero(npars, npars);
  Eigen::VectorXd Egrad;
  Egrad.setZero(npars);

  for(int iter=0; iter < nIter ; iter++){


    Eigen::SparseMatrix<double,0,int> K = Eigen::SparseMatrix<double,0,int>(Kobj->Q);

    // subsampling
    //sample Nlong values without replacement from 1:nrep
    std::random_shuffle(longInd.begin(), longInd.end());
    Kobj->gradient_init(nSubsample, nSim);
    for(int ilong = 0; ilong < nSubsample; ilong++ )
    {
      int i = longInd[ilong];
      Eigen::SparseMatrix<double,0,int> A = As[i];
      
      Eigen::VectorXd  Y = errObj->simulate( Ys[i]);
	  mixobj->simulate(Y, i);
	  
	 for(int j =0; j < Kobj->d; j++)
    			z[j] =  normal(random_engine);
     
     process->simulate_V(i, rgig);
	 process->simulate(i, z, A, K, Y, Solver[i]);
	  
      for(int ii = 0; ii < nSim + nBurnin; ii ++)
      {
      	Eigen::VectorXd  res = Y;
      	//***************************************
      	//***************************************
      	//   building the residuals and sampling
      	//***************************************
      	//***************************************

      	// removing fixed effect from Y
      	mixobj->remove_cov(i, res);

    		res -= A * process->Xs[i];
  			//***********************************
      	// mixobj sampling
    		//***********************************
      	if(type_MeasurementError == "Normal")
    			mixobj->sampleU( i, res, 2 * log(errObj->sigma));
  			else
  			  mixobj->sampleU2( i, res, errObj->Vs[i].cwiseInverse(), 2 * log(errObj->sigma));

        mixobj->remove_inter(i, res);

      	//***********************************
    		// sampling processes
  			//***********************************
      	Eigen::SparseMatrix<double,0,int> Q = Eigen::SparseMatrix<double,0,int>(K.transpose());
    		Eigen::VectorXd iV(process->Vs[i].size());
  			iV.array() = process->Vs[i].array().inverse();
      	Q =  Q * iV.asDiagonal();
    		Q =  Q * K;
        for(int j =0; j < Kobj->d; j++)
    			z[j] =  normal(random_engine);

      	res += A * process->Xs[i];
      	//Sample X|Y, V, sigma

        if(type_MeasurementError == "Normal")
      		process->sample_X(i, z, res, Q, K, A, errObj->sigma, Solver[i]);
        else
          process->sample_Xv2( i, z, res, Q, K, A, errObj->sigma, Solver[i], errObj->Vs[i].cwiseInverse());

        res -= A * process->Xs[i];

        if(res.cwiseAbs().sum() > 1e16){
          throw("res outof bound\n");
        }

        // sample V| X
        process->sample_V(i, rgig, K);

        // random variance noise sampling
     		if(type_MeasurementError != "Normal"){
      	  errObj->sampleV(i, res);
     		}

      	//***************************************
      	//  computing gradients
      	//***************************************
      	if(ii >= nBurnin){

      		// mixobj gradient
      		mixobj->add_inter(i, res);
      	  if(type_MeasurementError != "Normal")
      		  mixobj->gradient2(i, res, errObj->Vs[i].cwiseInverse(), 2 * log(errObj->sigma), errObj->EiV);
          else
            mixobj->gradient(i, res, 2 * log(errObj->sigma));

      		mixobj->remove_inter(i, res);


      	  // measurent error  gradient
  			  errObj->gradient(i, res);

      	  // operator gradient
      		Kobj->gradient_add( process->Xs[i], process->Vs[i].cwiseInverse());

      		// process gradient
          if(type_MeasurementError != "Normal"){
              process->gradient_v2(i,
                                K,
                                A,
                                res,
                                errObj->sigma,
                                errObj->Vs[i].cwiseInverse(),
                                errObj->EiV,
                                Kobj->trace_variance(A));
            }else{
              process->gradient(i,
                                K,
                                A,
                                res,
                                errObj->sigma,
                                Kobj->trace_variance(A));
            }
      		}
      }
    }
    
    Eigen::VectorXd grad(npars);
    grad.segment(0, mixobj->npars) = mixobj->get_gradient();
    grad.segment(mixobj->npars, errObj->npars) = errObj->get_gradient();
    grad.segment(mixobj->npars + errObj->npars, process->npars) =  process -> get_gradient();
    grad.segment(mixobj->npars + errObj->npars + process->npars , Kobj->npars) =  Kobj -> get_gradient();
    if(silent == 0)
    	Rcpp::Rcout << "grad  = " << grad.transpose() << "\n";
    
    mixobj->clear_gradient();
	errObj->clear_gradient();
	process -> clear_gradient();
    Kobj ->  clear_gradient();
    grad.array() /= nSim;
    // correct for the samples
    grad.array() *= Ys.size() /  nSubsample; 
    
    
    Egrad2 += grad * grad.transpose();
    Egrad += grad;
  
  }
  
  
  Egrad2.array() /= nIter;
  Egrad.array()  /= nIter; 
  Eigen::MatrixXd Vgrad = Egrad2 - Egrad * Egrad.transpose();
  //Rcpp::Rcout << "Vgrad = \n" << Vgrad <<"\n";
  Eigen::MatrixXd invVgrad  = Vgrad.inverse();
  mixobj->set_covariance(invVgrad.block(0, 0, mixobj->npars, mixobj->npars));
  errObj->set_covariance(invVgrad.block(mixobj->npars, mixobj->npars, errObj->npars, errObj->npars));
  process->set_covariance(invVgrad.block(mixobj->npars +  errObj->npars , mixobj->npars +  errObj->npars,
  										  process->npars, process->npars));
  Kobj->set_covariance(invVgrad.block(mixobj->npars +  errObj->npars + process->npars , 
  									  mixobj->npars +  errObj->npars + process->npars,
  									  Kobj->npars,
  									  Kobj->npars));  
  
  
  
  
  if(silent == 0){
  	Rcpp::Rcout << " std = " << invVgrad.diagonal().cwiseSqrt() << "\n";
  	Rcpp::Rcout << "Done, storing results\n";
  	}
  // storing the results
  Rcpp::List out_list;
  out_list["pSubsample"]       = pSubsample;
  out_list["nIter"]            = nIter;
  out_list["nSim"]             = nSim;
  out_list["nBurnin"]          = nBurnin;
  out_list["silent"]           = silent;
  out_list["obs_list"]         = obs_list;
  out_list["Xs"]               = process->Xs;
  out_list["Vs"]               = process->Vs;

  Rcpp::List mixobj_list       = mixobj->toList();
  out_list["mixedEffect_list"] = mixobj_list;

  Rcpp::List errobj_list            = errObj->toList();
  out_list["measurementError_list"] = errobj_list;

  Rcpp::List process_list           = process->toList();
  out_list["processes_list"]        = process_list;

  Rcpp::List olist          = Kobj->output_list();
  out_list["operator_list"] = olist;
  return(out_list);
}
