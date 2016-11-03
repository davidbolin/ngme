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
#include "latentprocess.h"
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
	int nIter      = Rcpp::as< int > (in_list["nIter"]);
	int nSim       = Rcpp::as< int > (in_list["nSim"]);
  int nBurnin    = Rcpp::as< int > (in_list["nBurnin"] );
  int nBurnin_base    = Rcpp::as< int > (in_list["nBurnin_base"] );
  	int nBurnin_learningrate = nBurnin;
  	if(in_list.containsElementNamed("nBurnin_learningrate"))
  		nBurnin_learningrate    = Rcpp::as< int > (in_list["nBurnin_learningrate"] );
  	int silent     = Rcpp::as< int    > (in_list["silent"]);
  	double alpha     = Rcpp::as< double    > (in_list["alpha"]);
  	double step0     = Rcpp::as< double    > (in_list["step0"]);
  	int subsample_type = Rcpp::as< int    > (in_list["subsample_type"]);
  	unsigned long seed = 0;
  	if(in_list.containsElementNamed("seed"))
  		seed = Rcpp::as< unsigned long    > (in_list["seed"]);
  	double learning_rate = 0;
  	int process_active = 0;
  	if(in_list.containsElementNamed("processes_list"))
  		process_active = 1;
  	if(in_list.containsElementNamed("learning_rate"))
  		learning_rate = Rcpp::as< double    > (in_list["learning_rate"]);
	
  	double polyak_rate = -1;
  	if(in_list.containsElementNamed("polyak_rate"))
  		polyak_rate = Rcpp::as< double    > (in_list["polyak_rate"]);

  	int debug = 0;
	//**********************************
	//     setting up the main data
	//**********************************
	if(silent == 0){
	  Rcpp::Rcout << " Setup data\n";
	}
	Rcpp::List obs_list  = Rcpp::as<Rcpp::List> (in_list["obs_list"]);
	int nindv = obs_list.length(); //count number of patients
  int nSubsample = ceil(pSubsample * nindv);
	std::vector< Eigen::SparseMatrix<double,0,int> > As( nindv);
	std::vector< Eigen::VectorXd > Ys( nindv);
	std::vector< double > Ysize(nindv);
	std::vector< int > burnin_done(nindv);
	Eigen::VectorXd sampling_weights(nindv);
	int count;
	count = 0;
	for( List::iterator it = obs_list.begin(); it != obs_list.end(); ++it ) {
    	List obs_tmp = Rcpp::as<Rcpp::List>(*it);
    	if(process_active)
    		As[count] = Rcpp::as<Eigen::SparseMatrix<double,0,int> >(obs_tmp["A"]);
    	Ys[count] = Rcpp::as<Eigen::VectorXd>(obs_tmp["Y"]);
    	Ysize[count] = (double) Ys[count].size();
    	burnin_done[count] = 0;
    	if(subsample_type == 1){
    	 sampling_weights[count] = 1.0;
    	} else if(subsample_type == 2){
  	  		sampling_weights[count] = Ysize[count];
  		} else if (subsample_type == 3){ //Biased sampling
    		sampling_weights[count] = 1.0;
    	}
    	count++;
  }
	sampling_weights /= sampling_weights.sum();

	//***********************************
	//Debug setup
	//***********************************
	int sampleV = 1;
	if(in_list.containsElementNamed("sampleV"))
		sampleV = Rcpp::as<int> (in_list["sampleV"]);

	int sampleX = 1;
	if(in_list.containsElementNamed("sampleX"))
		sampleX = Rcpp::as<int> (in_list["sampleX"]);

	//**********************************
	//operator setup
	//***********************************
	if(silent == 0){
	  Rcpp::Rcout << " Setup operator\n";
	}
	operatorMatrix* Kobj;
	//Eigen::VectorXd h = Rcpp::as<Eigen::VectorXd>( operator_list["h"]);
	//Create solvers for each patient
	std::vector<  cholesky_solver >  Solver( nindv);
	Eigen::SparseMatrix<double, 0, int> Q, K;

	int common_grid = 1;
	if(process_active){
		Rcpp::List operator_list  = Rcpp::as<Rcpp::List> (in_list["operator_list"]);
		operator_list["nIter"] = nIter;
		std::string type_operator = Rcpp::as<std::string>(operator_list["type"]);

		operator_select(type_operator, &Kobj);
		Kobj->initFromList(operator_list, List::create(Rcpp::Named("use.chol") = 1));

		if(Kobj->nop>1)
	  		common_grid = 0;


		count = 0;
		for( List::iterator it = obs_list.begin(); it != obs_list.end(); ++it ) {
    		List obs_tmp = Rcpp::as<Rcpp::List>( *it);

    		if(common_grid == 1){
      			Solver[count].init(Kobj->d[0], 0, 0, 0);
      			K = Eigen::SparseMatrix<double,0,int>(Kobj->Q[0]);
    		} else {
      			Solver[count].init(Kobj->d[count], 0, 0, 0);
      			K = Eigen::SparseMatrix<double,0,int>(Kobj->Q[count]);
    		}
    		Q = K.transpose();
    		Q = Q * K;
  			Q = Q + As[count].transpose()*As[count];
    		Solver[count].analyze(Q);
  			Solver[count].compute(Q);
    		count++;
  		}
  }
	//**********************************
	// mixed effect setup
	//***********************************
	if(silent == 0){
	  Rcpp::Rcout << " Setup mixed effect\n";
	}
	Rcpp::List mixedEffect_list  = Rcpp::as<Rcpp::List> (in_list["mixedEffect_list"]);
	std::string type_mixedEffect = Rcpp::as<std::string> (mixedEffect_list["noise"]);
	MixedEffect *mixobj = NULL;
	if(type_mixedEffect == "Normal")
    mixobj = new NormalMixedEffect;
	else if(type_mixedEffect == "NIG")
		mixobj   = new NIGMixedEffect;
	mixobj->initFromList(mixedEffect_list);
	mixobj->setupStoreTracj(nIter);


  //**********************************
	// measurement setup
	//***********************************
	if(silent == 0){
	  Rcpp::Rcout << " Setup noise\n";
	}
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
	Process *process;
	if(process_active){
		if(silent == 0){
	  		Rcpp::Rcout << " Setup process\n";
		}
		Rcpp::List processes_list   = Rcpp::as<Rcpp::List>  (in_list["processes_list"]);
		Rcpp::List V_list           = Rcpp::as<Rcpp::List>  (processes_list["V"]);
		std::string type_processes  = Rcpp::as<std::string> (processes_list["noise"]);



  		if (type_processes != "Normal"){
  			process  = new GHProcess;
  		}else{ process  = new GaussianProcess;}

  		process->initFromList(processes_list, Kobj->h);
  		process->setupStoreTracj(nIter);
  	/*
  		Simulation objects
  	*/
  	}
  	std::mt19937 random_engine;
	std::normal_distribution<double> normal;
    std::default_random_engine subsample_generator;
  std::default_random_engine gammagenerator;
  gig rgig;

  if(in_list.containsElementNamed("seed")){
  	//rgig.seed(seed);
  	random_engine.seed(seed);
  }else{
  	random_engine.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());
  	//rgig.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());
  }

  rgig.seed(random_engine());
  gammagenerator.seed(random_engine());
  subsample_generator.seed(random_engine());
  Eigen::VectorXd  z;
  Eigen::VectorXd b, Ysim;
  if(process_active){
	z.setZero(Kobj->d[0]);
	b.setZero(Kobj->d[0]);
	}

  std::vector<int> longInd;
  std::vector<Eigen::VectorXd> Vmean;
  Eigen::VectorXd count_vec(nindv);
  Vmean.resize(nindv);
  count_vec.setZero(nindv);
  if(process_active){
  for(int i = 0; i < nindv; i++ ){
    	if(common_grid){
      	Vmean[i].setZero(Kobj->h[0].size());
    	} else {
      	Vmean[i].setZero(Kobj->h[i].size());
    	}
  	}
  }

 for (int i=0; i< nindv; i++) longInd.push_back(i);

  int npars =   mixobj->npars + errObj->npars;
  if(process_active)
    npars += process->npars; //+ Kobj->npars;

  for(int iter=0; iter < nIter; iter++){
    /*
      printing output
    */
    if(silent == 0){
      Rcpp::Rcout << "i = " << iter << ": \n";
      if(process_active){
      	process->printIter();
      	Rcpp::Rcout << "\n";
      	Kobj->print_parameters();
      	Rcpp::Rcout << "\n";
      }
      mixobj->printIter();
      Rcpp::Rcout << "\n";
      errObj->printIter();
      Rcpp::Rcout << "\n";
    }


    // subsampling
    if(subsample_type == 1){
      //Uniform sampling without replacement from 1:nrep
      std::shuffle(longInd.begin(), longInd.end(), gammagenerator);
    } else if(subsample_type == 2){
      //weighted by number of samples
      std::discrete_distribution<int> distribution(Ysize.begin(), Ysize.end());
      for (int i=0; i<nSubsample; ++i) {
        longInd[i] = distribution(gammagenerator);
      }
      //std::cout << std::endl << "P: ";
      //for (double x:distribution.probabilities()) std::cout << x << " ";
    }

    double burnin_rate = 0;
    if(process_active)
    	Kobj->gradient_init(nSubsample,nSim);

    	Eigen::VectorXd  grad_last(npars); // grad_last helper vector
      grad_last.setZero(npars);

    Eigen::MatrixXd grad_outer(nSubsample, npars); // outside person variation (subsampling)
    Eigen::MatrixXd Vgrad_inner(npars, npars); // outside person variation (subsampling)
    Eigen::VectorXd Ebias_inner(npars);
    Ebias_inner.array() *= 0;
    Vgrad_inner.setZero(npars, npars);
    for(int ilong = 0; ilong < nSubsample; ilong++ )
    {
      int i = longInd[ilong];

      double w = sampling_weights[i]; //TODO:: ADD WEIGHTING BY W TO ALL GRADIENTS!
      Eigen::SparseMatrix<double,0,int> A;
      if(process_active)
      	A = As[i];
      Eigen::VectorXd  Y = Ys[i];
      if(process_active){
      	z.setZero(Kobj->d[0]);
      	if(common_grid == 0)
      		z.setZero(Kobj->d[i]);
      }

      int n_simulations = nSim;
      int burnin_done_i = nBurnin_base;
      if(burnin_done[i] == 0){
        burnin_rate +=1;
        n_simulations += nBurnin;
        burnin_done_i = nBurnin;
        burnin_done[i] = 1;
      }else {n_simulations += nBurnin_base;}

      int count_inner = 0;
      Eigen::MatrixXd grad_inner(npars, nSim); // within person variation  (Gibbs)
      for(int ii = 0; ii < n_simulations; ii ++)
      {
      	Eigen::VectorXd  res = Y;
      	//***************************************
      	//***************************************
      	//   building the residuals and sampling
      	//***************************************
      	//***************************************


      	// removing fixed effect from Y
      	mixobj->remove_cov(i, res);
      	if(process_active)
    	  res -= A * process->Xs[i];
  		//***********************************
      	// mixobj sampling
    	//***********************************

       if(debug)
       		Rcpp::Rcout << "estimate::sample mix \n";
      	if(type_MeasurementError == "Normal")
    			mixobj->sampleU( i, res, 2 * log(errObj->sigma));
  			else
  			  mixobj->sampleU2( i, res, errObj->Vs[i].cwiseInverse(), 2 * log(errObj->sigma));
		//Rcpp::Rcout << "V[" << i << "," << ii <<  "] = " << ((NIGMixedEffect*) mixobj)->V[i] << "\n";
        mixobj->remove_inter(i, res);
      	//***********************************
    		// sampling processes
  			//***********************************
  			if(process_active){
  			Eigen::VectorXd iV(process->Vs[i].size());
  			iV.array() = process->Vs[i].array().inverse();
  			if(common_grid){
  			  K = Eigen::SparseMatrix<double,0,int>(Kobj->Q[0]);
  			} else {
  			  K = Eigen::SparseMatrix<double,0,int>(Kobj->Q[i]);
  			}
    		Q = K.transpose();
      		Q =  Q * iV.asDiagonal();
    		Q =  Q * K;
        for(int j =0; j < K.rows(); j++)
    			z[j] =  normal(random_engine);

      	res += A * process->Xs[i];
		if(debug)
       	  Rcpp::Rcout << "estimate::sample X\n";
      	//Sample X|Y, V, sigma
		    if(sampleX){
          		if(type_MeasurementError == "Normal")
      				process->sample_X(i, z, res, Q, K, A, errObj->sigma, Solver[i]);
        		else
          		process->sample_Xv2( i, z, res, Q, K, A, errObj->sigma, Solver[i], errObj->Vs[i].cwiseInverse());
		    	}
        	res -= A * process->Xs[i];
			if(res.cwiseAbs().sum() > 1e16){
        		Rcpp::Rcout << "MAX(process->Vs[i]^-1) = " << process->Vs[i].cwiseInverse().maxCoeff() << "\n";
        		Rcpp::Rcout << "Max process->Xs[i]= " << process->Xs[i].maxCoeff() << "\n";
        		Rcpp::Rcout << "Min process->Xs[i]= " << process->Xs[i].minCoeff() << "\n";

        		Rcpp::Rcout << "res out of bound\n";
        		throw("res outof bound\n");
        	}

        	if(debug)
       	  	Rcpp::Rcout << "estimate::sample V\n";
        	// sample V| X
        	if(sampleV)
         	 process->sample_V(i, rgig, K);
		}
        if(debug)
         	Rcpp::Rcout << "estimate::sample err V\n";
        // random variance noise sampling
     	  if(type_MeasurementError != "Normal"){
      	  errObj->sampleV(i, res);
     		}

      	//***************************************
      	//  computing gradients
      	//***************************************

        if(debug)
       		Rcpp::Rcout << "estimate::gradient calc \n";
      	if(ii >= burnin_done_i){

      	  //TODO:: ADDD SCALING WITH W FOR MIX GRADIENT
      		// mixobj gradient
      		mixobj->add_inter(i, res);
      	  if(type_MeasurementError != "Normal")
      		  mixobj->gradient2(i, res, errObj->Vs[i].cwiseInverse(), 2 * log(errObj->sigma), errObj->EiV);
          else
            mixobj->gradient(i, res, 2 * log(errObj->sigma));

      		mixobj->remove_inter(i, res);

		  // measurent error  gradient
      	  //TODO:: ADDD SCALING WITH W FOR ERROR GRADIENT
  			  errObj->gradient(i, res);




		  if(process_active){

		  	// operator gradient
      		Kobj->gradient_add( process->Xs[i],
      							process->Vs[i].cwiseInverse(),
      							process->mean_X(i),i);

      		// process gradient
      		res += A * process->Xs[i];

          if(type_MeasurementError != "Normal"){
                //TODO:: ADDD SCALING WITH W FOR PROCESS GRADIENT
              process->gradient_v2(i,
                                K,
                                A,
                                res,
                                errObj->sigma,
                                errObj->Vs[i].cwiseInverse(),
                                errObj->EiV,
                                Kobj->trace_variance(A, i));
            }else{
              process->gradient(i,
                                K,
                                A,
                                res,
                                errObj->sigma,
                                Kobj->trace_variance(A, i));
            }
        }
		    // collects the gradient for computing estimate of variances
		    grad_inner.block(0, count_inner,
                          mixobj->npars, 1)  = mixobj->get_gradient();
		    grad_inner.block(mixobj->npars,count_inner,
                          errObj->npars, 1) = errObj->get_gradient();
		    if(process_active){
		      grad_inner.block(mixobj->npars + errObj->npars,count_inner,
                         process->npars, 1).array() = process->get_gradient().array();
		      //grad_inner.block(mixobj->npars + errObj->npars + process->npars,count_inner,
              //              Kobj->npars, 1) = Kobj -> get_gradient();
			  }
			  // adjusting so we get last gradient not cumsum
			  Eigen::VectorXd grad_last_temp = grad_inner.col(count_inner);
			  grad_inner.col(count_inner).array() -= grad_last.array();
			  grad_last = grad_last_temp;
		    count_inner++;
      }

      }

      Eigen::VectorXd Mgrad_inner = grad_inner.rowwise().mean();
      grad_outer.row(ilong) = Mgrad_inner;

      Eigen::MatrixXd centered = grad_inner.colwise() - Mgrad_inner;

      Ebias_inner.array() += centered.col(nSim-1).array();
      Ebias_inner.array() -= centered.col(0).array();
      Eigen::MatrixXd cov = (centered * centered.transpose()) / double(grad_inner.cols() - 1);
      Vgrad_inner.array() += cov.array();

      if(process_active)
        Vmean[i] += process->Vs[i];

  		count_vec[i] += 1;

    }
    //Compute variance due to subsampling:
    Vgrad_inner.array() /=  double(grad_outer.rows() - 1);
    Eigen::VectorXd std_grad = Vgrad_inner.diagonal();
    std_grad.array() /= (nSim * nSubsample);

    //Compute Bias due to not converged
    Ebias_inner.array() /=  double(grad_outer.rows() );
    Eigen::MatrixXd centered = grad_outer.rowwise() - grad_outer.colwise().mean();


	  //Compute MC variance
	  Eigen::MatrixXd cov = (centered.transpose() * centered) / double(grad_outer.rows() - 1);
	  Eigen::VectorXd std_grad_outer = cov.diagonal();

	  //std_grad_outer.array() /= nSubsample;
	  std_grad_outer.array() -= std_grad.array();

	  std_grad.array() /= ( nSubsample);
	  std_grad_outer.array() /= ( nSubsample);
	  double sub_factor = (1-nSubsample/nindv);
	  std_grad_outer.array() *= sub_factor;

	  //compute total variance
	  Eigen::VectorXd grad_var = std_grad;
	  grad_var.array() += std_grad_outer.array();


	  std_grad_outer.array()  = std_grad_outer.array().sqrt();
	  std_grad.array()  = std_grad.array().sqrt();
	  grad_var.array()  = grad_var.array().sqrt();

	  if(silent == 0){
	    Rcpp::Rcout << "Gibbs std = " << std_grad.transpose() << "\n";
	    Rcpp::Rcout << "Outer std = " << std_grad_outer.transpose() << "\n";
	    Rcpp::Rcout << "Total std = " << grad_var.transpose() << "\n";
	    //Rcpp::Rcout << "E[bias]    = " << Ebias_inner.transpose()   << "\n";
	  }


    if(silent == 0)
      Rcpp::Rcout << "Burnin percentage " << burnin_rate/nSubsample << "\n";

    //**********************************
  	//  gradient step
	//***********************************

    if(debug)
      Rcpp::Rcout << "estimate::theta  step\n";
    double stepsize = step0 / pow(iter + 1, alpha);

    double polyak_rate_temp = polyak_rate;
    double learning_rate_temp  =learning_rate;
    if(polyak_rate == 0)
    	polyak_rate_temp = 1./ (iter + 1);
    if(iter < nBurnin_learningrate)
    	learning_rate_temp = 0;
    if(debug)
    	Rcpp::Rcout << "polyak_rate_temp = " << polyak_rate_temp <<"\n";

    mixobj->step_theta(stepsize,  learning_rate_temp, polyak_rate_temp);
    errObj->step_theta(stepsize,                   0, polyak_rate_temp);
    if(process_active){
    	Kobj->step_theta(stepsize,    learning_rate_temp, polyak_rate_temp);
    	process->step_theta(stepsize, learning_rate_temp, polyak_rate_temp);
    }


    if(debug)
      Rcpp::Rcout << "estimate::theta step done\n";

  }

  for(int i = 0; i < nindv; i++ )
  	Vmean[i].array() /= count_vec[i];

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
  if(process_active){
  	out_list["Xs"]               = process->Xs;
  	out_list["Vs"]               = process->Vs;
  	out_list["Vmean"]            = Vmean;
  Rcpp::List process_list           = process->toList();
  out_list["processes_list"]        = process_list;
  Rcpp::List olist          = Kobj->output_list();
  out_list["operator_list"] = olist;
  }
  Rcpp::List mixobj_list       = mixobj->toList();
  out_list["mixedEffect_list"] = mixobj_list;


  Rcpp::List errobj_list            = errObj->toList();
  out_list["measurementError_list"] = errobj_list;


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
  int common_grid = Rcpp::as< int > (in_list["common.grid"]);
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

	//Create solvers for each patient
	std::vector<  cholesky_solver >  Solver( nindv);
	Eigen::SparseMatrix<double, 0, int> Q,K;

	count = 0;
	for( List::iterator it = obs_list.begin(); it != obs_list.end(); ++it ) {
    	List obs_tmp = Rcpp::as<Rcpp::List>( *it);
    	Solver[count].init(Kobj->d[count], 0, 0, 0);
      if(common_grid==1){
        K = Eigen::SparseMatrix<double,0,int>(Kobj->Q[0]);
      } else {
        K = Eigen::SparseMatrix<double,0,int>(Kobj->Q[count]);
      }
  		Q = K.transpose();
    	Q = Q  * K;
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

  process->initFromList(processes_list, Kobj->h);
  /*
  Simulation objects
  */
  std::mt19937 random_engine;
	std::normal_distribution<double> normal;
  std::default_random_engine gammagenerator;
  gig rgig;
	rgig.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());
	random_engine.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());
	gammagenerator.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());
  Eigen::VectorXd  z;

  Eigen::VectorXd  Ysim;

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


    // subsampling
    //sample Nlong values without replacement from 1:nrep
    std::shuffle(longInd.begin(), longInd.end(), random_engine);
    Kobj->gradient_init(nSubsample, nSim);
    for(int ilong = 0; ilong < nSubsample; ilong++ )
    {
      int i = longInd[ilong];
      Eigen::SparseMatrix<double,0,int> A = As[i];
      if(common_grid == 1){
        K = Eigen::SparseMatrix<double,0,int>(Kobj->Q[0]);
      } else {
        K = Eigen::SparseMatrix<double,0,int>(Kobj->Q[i]);
      }

      Eigen::VectorXd  Y = errObj->simulate( Ys[i]);
	  mixobj->simulate(Y, i);

   z.setZero(Kobj->d[i]);
	 for(int j =0; j < K.rows(); j++)
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
      	Q = K.transpose();
    		Eigen::VectorXd iV(process->Vs[i].size());
  			iV.array() = process->Vs[i].array().inverse();
      	Q =  Q * iV.asDiagonal();
    		Q =  Q * K;
        for(int j =0; j < Kobj->d[i]; j++)
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
      		Kobj->gradient_add( process->Xs[i],
                              process->Vs[i].cwiseInverse(),
                              process->mean_X(i),
                              i);

      		// process gradient
          if(type_MeasurementError != "Normal"){
              process->gradient_v2(i,
                                K,
                                A,
                                res,
                                errObj->sigma,
                                errObj->Vs[i].cwiseInverse(),
                                errObj->EiV,
                                Kobj->trace_variance(A, i));
            }else{
              process->gradient(i,
                                K,
                                A,
                                res,
                                errObj->sigma,
                                Kobj->trace_variance(A, i));
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
  out_list["obs_list"]         = obs_list;
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
