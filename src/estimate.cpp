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
#include "subSampleDiagnostic.h"
#include "sample.h"
using namespace Rcpp;

double estDigamma(double x)
{
  return(R::digamma(x));
}



/*
	internal function handling all sampling
	Y - observation for indivual i

*/
Eigen::VectorXd GibbsSampling(int i,
				   Eigen::VectorXd&  Y,
				   Eigen::SparseMatrix<double,0,int>& A,
				   int               sampleX,
				   int               sampleV,
				   int               process_active,
				   int common_grid,
				   MixedEffect       & mixobj,
				   operatorMatrix    & Kobj,
				   MeasurementError  & errObj,
				   Process           & process,
				   int debug,
				   Eigen::VectorXd & z,
				   gig & rgig,
				   std::vector<  cholesky_solver > & Solver)
{
  Eigen::VectorXd b;
  if(process_active){
	  b.setZero(Kobj.d[0]);
	 }

  Eigen::VectorXd  res = Y;
  //***************************************
  //***************************************
  //   building the residuals and sampling
  //***************************************
  //***************************************

  // removing fixed effect from Y
  mixobj.remove_cov(i, res);
  if(process_active)
    res -= A * process.Xs[i];
  //***********************************
  // mixobj sampling
  //***********************************

  if(debug)
    Rcpp::Rcout << "estimate::sample mix \n";
  if(errObj.noise == "Normal")
  	mixobj.sampleU( i, res, 2 * log(errObj.sigma));
  else
	  mixobj.sampleU2( i, res, errObj.Vs[i].cwiseInverse(), 2 * log(errObj.sigma));
	//Rcpp::Rcout << "V[" << i << "," << ii <<  "] = " << ((NIGMixedEffect*) mixobj)->V[i] << "\n";
  mixobj.remove_inter(i, res);
  //***********************************
  // sampling processes
  //***********************************
	Eigen::SparseMatrix<double, 0, int> Q, K;
  if(process_active){
  	Eigen::VectorXd iV(process.Vs[i].size());
  	iV.array() = process.Vs[i].array().inverse();
  	if(common_grid){
			K = Eigen::SparseMatrix<double,0,int>(Kobj.Q[0]);
  	} else {
  		K = Eigen::SparseMatrix<double,0,int>(Kobj.Q[i]);
		}
    Q = K.transpose();
    Q =  Q * iV.asDiagonal();
  	Q =  Q * K;

    res += A * process.Xs[i];
	  if(debug)
     	Rcpp::Rcout << "estimate::sample X\n";

    //Sample X|Y, V, sigma
	  if(sampleX){
      if(errObj.noise == "Normal")
      	process.sample_X(i, z, res, Q, K, A, errObj.sigma, Solver[i]);
      else
        process.sample_Xv2( i, z, res, Q, K, A, errObj.sigma, Solver[i], errObj.Vs[i].cwiseInverse());
		}
    res -= A * process.Xs[i];
		if(res.cwiseAbs().sum() > 1e16){
      Rcpp::Rcout << "MAX(process->Vs[i]^-1) = " << process.Vs[i].cwiseInverse().maxCoeff() << "\n";
      Rcpp::Rcout << "Max process->Xs[i]= " << process.Xs[i].maxCoeff() << "\n";
    	Rcpp::Rcout << "Min process->Xs[i]= " << process.Xs[i].minCoeff() << "\n";

    	Rcpp::Rcout << "res out of bound\n";
    	throw("res outof bound\n");
    }

    if(debug)
      Rcpp::Rcout << "estimate::sample V\n";
    // sample V| X
    if(sampleV)
      process.sample_V(i, rgig, K);
	}
  if(debug)
    Rcpp::Rcout << "estimate::sample err V\n";
  // random variance noise sampling
  if(errObj.noise != "Normal"){
    errObj.sampleV(i, res);
  }
  return(res);
}

/*
	internal function handling all gradient caculations

*/

void grad_caculations(int i,
					            Eigen::VectorXd&  res,
				   	          Eigen::SparseMatrix<double,0,int>& A,
				              double w,
				              int common_grid,
				              int process_active,
				              MixedEffect       & mixobj,
				              operatorMatrix    & Kobj,
				              MeasurementError  & errObj,
				              Process           & process,
                      const int           estimate_fisher,
                      Eigen::MatrixXd   & Fisher_information,
                      int debug)
  {
    if(debug)
      Rcpp::Rcout << "estimate::gradcalc::mix \n";
  // mixobj gradient
  mixobj.add_inter(i, res);
  int use_EU = 1;
  if(estimate_fisher>0)
    use_EU = 0;
  if(errObj.noise != "Normal"){
    mixobj.gradient2(i,res,errObj.Vs[i].cwiseInverse(), 2 * log(errObj.sigma),errObj.EiV,w, use_EU);
      if(estimate_fisher)
        Fisher_information.block(0, 0, mixobj.npars + 1, mixobj.npars + 1) += mixobj.d2Given2(i,res,errObj.Vs[i].cwiseInverse(), 2 * log(errObj.sigma),errObj.EiV,w);
  }else{
    mixobj.gradient(i,res,2 * log(errObj.sigma),w, use_EU);
    if(estimate_fisher){
      Fisher_information.block(0, 0, mixobj.npars + 1, mixobj.npars + 1) += mixobj.d2Given(i,res,2 * log(errObj.sigma),w);
      }
  }

  mixobj.remove_inter(i, res);

    if(debug)
      Rcpp::Rcout << "estimate::gradcalc::errObj \n";
	// measurent error  gradient
  //TODO:: ADDD SCALING WITH W FOR ERROR GRADIENT
  errObj.gradient(i, res, w);
  if( (errObj.noise != "Normal") && (estimate_fisher > 0) && (errObj.npars > 1) )
    Fisher_information.block(mixobj.npars + 1, mixobj.npars + 1, errObj.npars - 1, errObj.npars - 1) += errObj.d2Given(i, res, w);


	if(process_active){

    if(debug)
      Rcpp::Rcout << "estimate::gradcalc::operator \n";
	  // operator gradient
    Kobj.gradient_add( process.Xs[i], process.Vs[i].cwiseInverse(), process.mean_X(i), i, w);

    if(estimate_fisher > 0 ){
      Fisher_information.block(mixobj.npars + errObj.npars, mixobj.npars + errObj.npars, Kobj.npars, Kobj.npars ) +=
      Kobj.d2Given(process.Xs[i], process.Vs[i].cwiseInverse(), process.mean_X(i), i, w);
    }
    // process gradient
    res += A * process.Xs[i];

    if(debug)
      Rcpp::Rcout << "estimate::gradcalc::process \n";
		Eigen::SparseMatrix<double, 0, int> K;
  	Eigen::VectorXd iV(process.Vs[i].size());
  	iV.array() = process.Vs[i].array().inverse();
  	if(common_grid){
  		K = Eigen::SparseMatrix<double,0,int>(Kobj.Q[0]);
  	} else {
  	  K = Eigen::SparseMatrix<double,0,int>(Kobj.Q[i]);
  	}
    if(process.type_process != "Normal"){
      if(errObj.noise != "Normal"){
        //TODO:: ADDD SCALING WITH W FOR PROCESS GRADIENT
        process.gradient_v2(i,K,A,res,errObj.sigma,
                            errObj.Vs[i].cwiseInverse(),
                            errObj.EiV,
                            Kobj.trace_variance(A, i),
                            w);

        if(estimate_fisher > 0 ){
          Fisher_information.block(mixobj.npars + errObj.npars + Kobj.npars,
                                 mixobj.npars + errObj.npars + Kobj.npars,
                                 process.npars,
                                 process.npars ) += process.d2Given_v2(i,
                                                                       K,
                                                                       A,
                                                                       res,
                                                                       errObj.sigma,
                                                                       errObj.Vs[i].cwiseInverse(),
                                                                       0,
                                                                       0,
                                                                       w);

             Eigen::MatrixXd Bf_t = Eigen::MatrixXd::Zero(0,0);
             if(mixobj.Bf.size() > 0)
                Bf_t = mixobj.Bf[i];
             Eigen::MatrixXd Br_t = Eigen::MatrixXd::Zero(0,0);
             if(mixobj.Br.size() > 0)
                Br_t = mixobj.Br[i];
             Eigen::VectorXd cross = process.d2Given_v2_cross(i,
                                                   K,
                                                   A,
                                                   res,
                                                   errObj.sigma,
                                                   Bf_t,
                                                   Br_t,
                                                   errObj.Vs[i].cwiseInverse(),
                                                   w);
             

            if(mixobj.Bf.size() > 0)
              Fisher_information.block(mixobj.npars + errObj.npars + Kobj.npars,
                                       0,
                                       1,
                                       Bf_t.cols()) += cross.head(Bf_t.cols()).transpose();
            if(mixobj.Br.size() > 0)
              Fisher_information.block(mixobj.npars + errObj.npars + Kobj.npars,
                                       Bf_t.cols(),
                                       1,
                                       Br_t.cols()) += cross.segment(Bf_t.cols(), Br_t.cols()).transpose();

              Fisher_information(mixobj.npars + errObj.npars + Kobj.npars,
                                       mixobj.npars) += cross(Bf_t.cols() + Br_t.cols());
              Fisher_information.block(0,
                                       mixobj.npars + errObj.npars + Kobj.npars,
                                       Bf_t.cols() + Br_t.cols() + 1,
                                       1
                                       ) =
                                      Fisher_information.block(mixobj.npars + errObj.npars + Kobj.npars,
                                        0,
                                       1,
                                       Bf_t.cols() + Br_t.cols() + 1).transpose();
        }
      }else{
        
        process.gradient(i,K,A,res,errObj.sigma, Kobj.trace_variance(A, i),w);
        if(estimate_fisher > 0 ){

          Fisher_information.block(mixobj.npars + errObj.npars + Kobj.npars,
                                 mixobj.npars + errObj.npars + Kobj.npars,
                                 process.npars,
                                 process.npars ) += process.d2Given(i,
                                                                    K,
                                                                    A,
                                                                    res,
                                                                    errObj.sigma,
                                                                    0.,
                                                                    w);

             Eigen::MatrixXd Bf_t = Eigen::MatrixXd::Zero(0,0);
             if(mixobj.Bf.size() > 0)
                Bf_t = mixobj.Bf[i];
             Eigen::MatrixXd Br_t = Eigen::MatrixXd::Zero(0,0);
             if(mixobj.Br.size() > 0)
                Br_t = mixobj.Br[i];
             Eigen::VectorXd cross = process.d2Given_cross(i,
                                                   K,
                                                   A,
                                                   res,
                                                   errObj.sigma,
                                                   Bf_t,
                                                   Br_t,
                                                   w);
             if(mixobj.Bf.size() > 0)
              Fisher_information.block(mixobj.npars + errObj.npars + Kobj.npars,
                                       0,
                                       1,
                                       Bf_t.cols()) += cross.head(Bf_t.cols()).transpose();

            if(mixobj.Br.size() > 0)
              Fisher_information.block(mixobj.npars + errObj.npars + Kobj.npars,
                                       Bf_t.cols(),
                                       1,
                                       Br_t.cols()) += cross.segment(Bf_t.cols(), Br_t.cols()).transpose();

              Fisher_information(mixobj.npars + errObj.npars + Kobj.npars,
                                       mixobj.npars) += cross(Bf_t.cols() + Br_t.cols());

              Fisher_information.block(0,
                                       mixobj.npars + errObj.npars + Kobj.npars,
                                       Bf_t.cols() + Br_t.cols() + 1,
                                       1
                                       ) =
                                      Fisher_information.block(mixobj.npars + errObj.npars + Kobj.npars,
                                        0,
                                       1,
                                       Bf_t.cols() + Br_t.cols() + 1).transpose();
          }
      }
    }
  }

}

// [[Rcpp::export]]
List estimateLong_cpp(Rcpp::List in_list)
{
	//**********************************
	//      basic parameter
	//**********************************

	int debug = 0;
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


  if(in_list.containsElementNamed("nBurnin_learningrate"))
  	nBurnin_learningrate    = Rcpp::as< int > (in_list["nBurnin_learningrate"] );

  int nPar_burnin = 0;
  if(in_list.containsElementNamed("nPar_burnin"))
    nPar_burnin    = Rcpp::as< int > (in_list["nPar_burnin"] );


	int estimate_fisher = 0;
	if(in_list.containsElementNamed("estimate_fisher"))
  	estimate_fisher    = Rcpp::as< int > (in_list["estimate_fisher"] );

	double polyak_rate = -1;
	if(in_list.containsElementNamed("polyak_rate"))
	  polyak_rate = Rcpp::as< double    > (in_list["polyak_rate"]);


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

	int count = 0;
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

  int nSubsample_group[2];
	double pSubsample2 = 0;
	Eigen::VectorXd free;

	Rcpp::List group_list;
	count = 0;
	int ngroup;
	std::vector<Eigen::VectorXd> groups;
	if(subsample_type == 3){
	  pSubsample2 = Rcpp::as< double > (in_list["pSubsample"]); // BUG? should it be pSubsample???
	} else if(subsample_type == 4){


    group_list = Rcpp::as<Rcpp::List> (in_list["group_list"]);
    ngroup = group_list.length();
    groups.resize(ngroup);
	  free = Rcpp::as<Eigen::VectorXd>(in_list["free_samples"]);

	  for( List::iterator it = group_list.begin(); it != group_list.end(); ++it ) {
	    groups[count] = Rcpp::as<Eigen::VectorXd>(*it);
	    Rcpp::Rcout << "group: "<< groups[count] << "\n";
	    count++;
	  }
	  int gsum = 0;
	  for(int i=0;i<groups.size();i++){
	    int ngroup = groups[i].size();
	    gsum += ngroup;
	  }
	  double gmean = 0;
	  if(groups.size()>0){
	    gmean = gsum/groups.size();
	  }

	  int nfree = free.size();
	  for(int i=0;i<nfree;i++){
	    Rcpp::Rcout << "free[i]: "<< free[i] << "\n";
	  }
	}
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
	Eigen::VectorXd  z;
	int common_grid = 1;
	if(process_active){
		Rcpp::List operator_list  = Rcpp::as<Rcpp::List> (in_list["operator_list"]);
		operator_list["nIter"] = nIter;
		std::string type_operator = Rcpp::as<std::string>(operator_list["type"]);

		operator_select(type_operator, &Kobj);
		Kobj->initFromList(operator_list, List::create(Rcpp::Named("use.chol") = 1));

		if(Kobj->nop>1)
	  	common_grid = 0;

		z.setZero(Kobj->d[0]);
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
	  Rcpp::Rcout << " Setup mixed effect: ";
	}
	Rcpp::List mixedEffect_list  = Rcpp::as<Rcpp::List> (in_list["mixedEffect_list"]);
	std::string type_mixedEffect = Rcpp::as<std::string> (mixedEffect_list["noise"]);
	MixedEffect *mixobj = NULL;
  //Rcpp::Rcout << "type_mixedEffect = " << type_mixedEffect << "\n";
	if(type_mixedEffect == "Normal"){
	  mixobj = new NormalMixedEffect;
	} else if(type_mixedEffect == "NIG") {
	  mixobj   = new NIGMixedEffect;
	} else {
	  Rcpp::Rcout << "Wrong mixed effect distribution";
	}



	if(silent == 0){
	  Rcpp::Rcout << " init";
	}
	mixobj->initFromList(mixedEffect_list);
	mixobj->setupStoreTracj(nIter);
	if(silent == 0){
	  Rcpp::Rcout << ", done\n";
	}

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
  else if(type_MeasurementError == "tdist")
    errObj = new IGMeasurementError;
  else
    errObj = new NIGMeasurementError;

  errObj->initFromList(measurementError_list);
  errObj->setupStoreTracj(nIter);
	//**********************************
	// stochastic processes setup
	//***********************************
	Process *process = NULL;
	std::string type_processes;
	if(process_active){
		if(silent == 0){
	  		Rcpp::Rcout << " Setup process\n";
		}
		Rcpp::List processes_list   = Rcpp::as<Rcpp::List>  (in_list["processes_list"]);
		Rcpp::List V_list           = Rcpp::as<Rcpp::List>  (processes_list["V"]);
		type_processes= Rcpp::as<std::string> (processes_list["noise"]);

  	if (type_processes != "Normal"){
  		process  = new GHProcess;
  	}else{
  	  process  = new GaussianProcess;
  	 }

  	process->initFromList(processes_list, Kobj->h);
    if(estimate_fisher > 0)
      process->useEV = 0;
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
  gig *rgig_pointer = &rgig;
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

  // For sampling we have the following:
  /*
  	weight   -> the weight each observation has given it is selected
  			  (note this not 1/p)
  	p[i]     -> the probability of selecting indvual [i]
  	m        -> expected number of samples we want for weighted sampling
  	selected -> 1, 0 variable keeping track of which individuals being selected
  */
  Eigen::VectorXd weight(nindv), p_inv(nindv);
  p_inv.setZero(nindv);
  weight.setZero(nindv);
  weight.array() += nindv / ((double) nSubsample);
  Eigen::VectorXd p(nindv);
  p.setZero(nindv);
  Eigen::VectorXd p_N(nindv);
  p_N.setZero(nindv);
  p_N.array() += 1. / nindv;
  std::vector<int> selected(nindv, 0);

  int par_burnin = 0;

 for (int i=0; i< nindv; i++) longInd.push_back(i);

  if(subsample_type == 4){
    groupSampling_weights (nSubsample,
                           groups,
                           free,
                           weight,
                           nSubsample_group);
  }


  int npars =   mixobj->npars + errObj->npars;
  if(process_active)
    npars += process->npars + Kobj->npars;

  Eigen::MatrixXd Fisher_information(npars, npars); // If computing the fisher information
	Fisher_information.setZero(npars, npars);

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
    if(iter < nPar_burnin)
      par_burnin = 1;

    if(debug)
      Rcpp::Rcout << "estimate::subsample \n";
	  int nSubsample_i = nSubsample;
    // subsampling
    if(subsample_type == 1){
      //Uniform sampling without replacement from 1:nrep
      std::shuffle(longInd.begin(), longInd.end(), gammagenerator);
    } else if(subsample_type == 2){
      //weighted by number of samples
      std::discrete_distribution<int> distribution(Ysize.begin(), Ysize.end());
      for (int i=0; i<nSubsample_i; ++i) {
        longInd[i] = distribution(gammagenerator);
      }
      //std::cout << std::endl << "P: ";
      //for (double x:distribution.probabilities()) std::cout << x << " ";
    }else if(subsample_type == 3){
    	//longInd
    	longInd.resize(nSubsample, 0);
    	std::fill(longInd.begin(), longInd.end(), 0);
    	std::fill(selected.begin(), selected.end(), 0);
    	int m = int(pSubsample2 * nindv); // expected number of samples from the first part
    	weight.setZero(nindv);
    	if(iter <= 10){
    		ProbSampleNoReplace(m+nSubsample, p_N, longInd, selected);
    		weight.array() += nindv / ((double) (nSubsample + m));
    		nSubsample_i = nSubsample + m;
    	}else{
    		ProbSampleNoReplace(nSubsample, p_N, longInd, selected);
    		weight.array() += nindv / ((double) nSubsample);
    		nSubsample_i = nSubsample;
    	}

    	p = p_N;
    	if(iter > 10){
    		nSubsample_i += poissonSampling_internal( m,
					  								  p_inv,
					  								  weight,
					  								  longInd,
					  								  selected);
			}
			double w_sum = 0;
			for(int ilong = 0; ilong < nSubsample_i; ilong++ )
				w_sum += weight[longInd[ilong]];

    } else if(subsample_type ==4){
      nSubsample_i = nSubsample;

      longInd.clear();
      std::fill(longInd.begin(), longInd.end(), 0);

      //groupSampling_internal(groups, free, longInd,gammagenerator);
      //Rcpp::Rcout << "nSubsample_group = (" << nSubsample_group[0] << "," << nSubsample_group[1]<< ")\n";
      groupSampling_sampling(nSubsample_group,
                             groups,
                             free,
                             longInd,
                             gammagenerator);
      nSubsample_i = longInd.size();
    }

    double burnin_rate = 0;
    if(process_active)
    	Kobj->gradient_init(nSubsample_i,nSim);

    Eigen::VectorXd  grad_last(npars); // grad_last helper vector
    grad_last.setZero(npars);

    Eigen::MatrixXd grad_outer(nSubsample_i, npars); // outside person variation (subsampling)
    Eigen::MatrixXd grad_outer_unweighted(nSubsample_i, npars); // outside person variation (subsampling)

    Eigen::MatrixXd Vgrad_inner(npars, npars); // outside person variation (subsampling)
    Eigen::VectorXd Ebias_inner(npars);
    Ebias_inner.array() *= 0;
    Vgrad_inner.setZero(npars, npars);
    if(debug)
      Rcpp::Rcout << "estimate::start patient loop \n";

    for(int ilong = 0; ilong < nSubsample_i; ilong++ )
    {
      int i = longInd[ilong];
      if(debug)
        Rcpp::Rcout << "i = " << i << " "<< weight[i] << "\n";

      Eigen::SparseMatrix<double,0,int> A;
      if(process_active)
      	A = As[i];
      Eigen::VectorXd  Y = Ys[i];
      if(estimate_fisher == 1){
      	if(process_active){
      		if(common_grid == 1){
        		K = Eigen::SparseMatrix<double,0,int>(Kobj->Q[0]);
      		} else {
        		K = Eigen::SparseMatrix<double,0,int>(Kobj->Q[i]);
      		}
      	}
      	if(debug)
      	  Rcpp::Rcout << "estimate::simulate error\n";

        Y = errObj->simulate( Ys[i]);
        if(debug)
          Rcpp::Rcout << "estimate::simulate mixed\n";
	      mixobj->simulate(Y, i);

	      if(process_active){
	        if(debug)
	          Rcpp::Rcout << "estimate::simulate process";

	        z.setZero(K.rows());
	        for(int j =0; j < K.rows(); j++)
    			  z[j] =  normal(random_engine);
	        if (type_processes != "Normal"){
	          if(debug)
	            Rcpp::Rcout << ", noise ";
     	      process->simulate_V(i, rgig);
	        }
	        if(debug)
	          Rcpp::Rcout << ", X ";

          process->simulate(i, z, A, K, Y, Solver[i]);
	      }
	      if(debug)
	        Rcpp::Rcout << ", done \n ";
      }

      int n_simulations = nSim;
      int burnin_done_i = nBurnin_base;
      if(burnin_done[i] == 0){
        burnin_rate +=1;
        burnin_done_i = nBurnin;
        n_simulations += burnin_done_i;
        burnin_done[i] = 1;
      }else {burnin_done_i = 0;}

	  if(estimate_fisher)
	  	burnin_done[i] = 0;
    int count_inner = 0;
    Eigen::MatrixXd grad_inner(npars, nSim); // within person variation  (Gibbs)

    for(int ii = 0; ii < n_simulations; ii ++)
    {
      for(int j =0; j < K.rows(); j++)
    		z[j] =  normal(random_engine);

     	  Eigen::VectorXd res =  GibbsSampling(i,
				   					  Y,
				   					  A,
				  					  sampleX,
				  					  sampleV,
				   					  process_active,
				   					  common_grid,
				   					  *mixobj,
				   					  *Kobj,
				   					  *errObj,
				   					  *process,
				   					  debug,
				   					  z,
				   					  *rgig_pointer,
				   					  Solver);

      	//***************************************
      	//  computing gradients
      	//***************************************
        if(debug)
       		Rcpp::Rcout << "estimate::gradient calc \n";
      	if(ii >= burnin_done_i){
	    	  grad_caculations(i,
			              				res,
				   	    		        A,
				     	  	          weight[i],
				      	  	        common_grid,
				      		          process_active,
				   			            *mixobj,
				   			            *Kobj,
				  			            *errObj,
				   			            *process,
                            estimate_fisher,
                            Fisher_information,
                            debug);

		      // collects the gradient for computing estimate of variances
		      grad_inner.block(0, count_inner,
                            mixobj->npars, 1)  = mixobj->get_gradient();
		      grad_inner.block(mixobj->npars,count_inner,
                            errObj->npars, 1) = errObj->get_gradient();
		      if(process_active){
		        grad_inner.block(mixobj->npars + errObj->npars,count_inner,
                             process->npars, 1).array() = process->get_gradient().array();
		        grad_inner.block(mixobj->npars + errObj->npars + process->npars,count_inner,
                            Kobj->npars, 1) = Kobj -> get_gradient();

			    }
			    // adjusting so we get last gradient not cumsum
			    Eigen::VectorXd grad_last_temp = grad_inner.col(count_inner);
			    grad_inner.col(count_inner).array() -= grad_last.array();
          Eigen::MatrixXd Fisher_temp  = 0.5 * grad_inner.col(count_inner)*grad_inner.col(count_inner).transpose() /  weight[i];
			    Fisher_information -=  Fisher_temp + Fisher_temp.transpose();
          grad_last = grad_last_temp;
		      count_inner++;
        }
      }
      // TODO::redo and move the above note i is the same always, and
	    // we need both weighted unweighted
	    // grad_inner.array() /= weight[i]; // dont use weight here(?)
	    // do grad_outer_unweighted
      Eigen::VectorXd Mgrad_inner = grad_inner.rowwise().mean();
      grad_outer.row(ilong) = Mgrad_inner;
	    grad_outer_unweighted.row(ilong) = Mgrad_inner;
	    grad_outer_unweighted.row(ilong) /= weight[i];
      Eigen::MatrixXd Fisher_add  = 0.5 * nSim  * (Mgrad_inner/weight[i]) * Mgrad_inner.transpose();
      Fisher_information +=  Fisher_add + Fisher_add.transpose() ;

      Eigen::MatrixXd centered = grad_inner.colwise() - Mgrad_inner;
      Ebias_inner.array() += centered.col(nSim-1).array();
      Ebias_inner.array() -= centered.col(0).array();
      Eigen::MatrixXd cov = (centered * centered.transpose()) /  (weight[i]* double(grad_inner.cols() - 1));
      Vgrad_inner.array() += cov.array();
      if(process_active)
        Vmean[i] += process->Vs[i];
  		count_vec[i] += 1;

    }

    // update weights given the gradient
    // change here to unweighted gradient!
    if(subsample_type == 3){
    	Eigen::VectorXd W = gradientWeight(grad_outer_unweighted);
		  for( int id = 0; id < nSubsample_i; id++)
			  p_inv[longInd[id]] = W[id];
	  }
 	  if(silent == 0){
		  subSampleDiag(Vgrad_inner,
			        		  grad_outer,
					          Ebias_inner,
				            nSim * nSubsample_i,
					          nSubsample_i / nindv);
	  }


    if(silent == 0)
      Rcpp::Rcout << "Burnin percentage " << burnin_rate/nSubsample_i << "\n";

    //**********************************
  	//  gradient step
	  //***********************************

	  if(estimate_fisher == 0){
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

    	//mixobj->step_theta(stepsize,  learning_rate_temp, polyak_rate_temp);
    	mixobj->step_theta(stepsize,                   0, polyak_rate_temp);
    	errObj->step_theta(stepsize,                   0, polyak_rate_temp);
    	if(process_active){
    		Kobj->step_theta(stepsize,    learning_rate_temp, polyak_rate_temp);
    		process->step_theta(stepsize, learning_rate_temp, polyak_rate_temp);
    	}
    	if(debug)
      Rcpp::Rcout << "estimate::theta step done\n";
	  }else{
	 	  mixobj->clear_gradient();
		  errObj->clear_gradient();
		  if(process_active){
			  process -> clear_gradient();
    	  Kobj ->  clear_gradient();
      }
	  }
  }

  for(int i = 0; i < nindv; i++ )
  	Vmean[i].array() /= count_vec[i];

  if(silent == 0)
  	Rcpp::Rcout << "Done, storing results\n";
  // storing the results
  	Fisher_information.array()  /= (nIter * nSim);

  if(estimate_fisher > 0 ){

    Eigen::MatrixXd cov_est  = Fisher_information.inverse();
    mixobj->set_covariance(cov_est.block(0, 0, mixobj->npars, mixobj->npars));
  	errObj->set_covariance(cov_est.block(mixobj->npars, mixobj->npars, errObj->npars, errObj->npars));
  	if(process_active){
  	process->set_covariance(cov_est.block(mixobj->npars +  errObj->npars , mixobj->npars +  errObj->npars,
  										  process->npars, process->npars));
  		//Kobj->set_covariance(cov_est.block(mixobj->npars +  errObj->npars + process->npars ,
  		//								  mixobj->npars +  errObj->npars + process->npars,
  		//								  Kobj->npars,
  		//								  Kobj->npars));
  	}

  }
  Rcpp::List out_list;

  std::vector<double> parameter;
  mixobj->get_param(parameter);

  errObj->get_param(parameter);

  
  Rcpp::StringVector param_names;
  mixobj->get_param_names(param_names);
  errObj->get_param_names(param_names);
  if(process_active)
  {
      Kobj->get_param_names(param_names);
      Kobj->get_param(parameter);
      process->get_param_names(param_names);
      process->get_param(parameter);
  }
  Rcpp::NumericVector parameter_out(Rcpp::wrap(parameter));
  parameter_out.names()  = param_names;
  out_list["parameter"]     = parameter_out;
  if(estimate_fisher){
    Rcpp::NumericMatrix F(Rcpp::wrap(Fisher_information));
    
    Rcpp::rownames(F) = param_names;
    Rcpp::colnames(F) = param_names;
    out_list["FisherMatrix"]     = F;

  }

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
	Rcpp::Rcout << "estimateFisher is depricated \n";
	throw("function is depericated");
  Rcpp::List out_list;
  return(out_list);
}
