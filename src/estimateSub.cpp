#include "estimate_util.h"
using namespace Rcpp;

// [[Rcpp::export]]
List estimateProcesses_cpp(Rcpp::List in_list)
{
  //**********************************
  //      basic parameter
  //**********************************

  int nIter      = Rcpp::as< int > (in_list["nIter"]);
  int nSim       = Rcpp::as< int > (in_list["nSim"]);
  int nBurnin    = Rcpp::as< int > (in_list["nBurnin"] );

  double step0     = Rcpp::as< double    > (in_list["step0"]);
  
  unsigned long seed = 0;
  if(in_list.containsElementNamed("seed"))
    seed = Rcpp::as< unsigned long    > (in_list["seed"]);
  double learning_rate = 0;
  int process_active = 0;

    Rcpp::List operator_list  = Rcpp::as<Rcpp::List> (in_list["operator_list"]);
    operator_list["nIter"] = nIter;
    std::string type_operator = Rcpp::as<std::string>(operator_list["type"]);
    
    operatorMatrix* Kobj;
    operator_select(type_operator, &Kobj);
    Kobj->initFromList(operator_list, List::create(Rcpp::Named("use.chol") = 1));

    int count = 0;
    Eigen::SparseMatrix<double,0,int> A;


  //**********************************
  // stochastic processes setup
  //***********************************
  Process *process = NULL;
  std::string type_processes;
    Rcpp::List processes_list   = Rcpp::as<Rcpp::List>  (in_list["processes_list"]);
    Rcpp::List V_list           = Rcpp::as<Rcpp::List>  (processes_list["V"]);
    type_processes= Rcpp::as<std::string> (processes_list["noise"]);

    if (type_processes != "Normal"){
      process  = new GHProcess;
    }else{
      process  = new GaussianProcess;
    }

    process->initFromList(processes_list, Kobj->h);
    process->setupStoreTracj(nIter);
    /*
    Simulation objects
    */
  
  std::mt19937 random_engine;
  std::normal_distribution<double> normal;
  std::default_random_engine gammagenerator;
  std::default_random_engine subsample_generator;
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
  subsampler* sampler = new subsampler;
  sampler->initFromList(in_list);
  

  int par_burnin = 0;
  //Start main loop
  for(int iter=0; iter < nIter; iter++){
    /*
    printing output
    */

    
    
      
    sampler->sample(iter,subsample_generator);
    
    //loop over patients/replicates
    for(int ilong = 0; ilong < sampler->nSubsample_i; ilong++ ){

      

      int i = sampler->longInd[ilong];

      int n_simulations = nSim;
      
      //Gibbs sampling
      for(int ii = 0; ii < n_simulations; ii ++)
      {
        Eigen::SparseMatrix<double, 0, int> K;
         K = Eigen::SparseMatrix<double,0,int>(Kobj->Q[i]);
        process->sample_V(i, rgig, K);
        Eigen::VectorXd  res0;
        Eigen::VectorXd  V0;
        //***************************************
        //  computing gradients
        //***************************************
        if(iter >= nBurnin)
            process->gradient_v2(i,K,A,res0,0,
                            V0,
                            0,
                            0,
                            sampler->weight[i]);




    }
    }
     if(iter >= nBurnin){
      double stepsize = step0 / pow(iter + 1, .5);
      process->step_theta(stepsize, 1., 0);
      process->clear_gradient();
      process->printIter();
      Rcpp::Rcout << "\n";
    }
    }
      Rcpp::List out_list;
   out_list["Vs"]               = process->Vs;
   Rcpp::List process_list      = process->toList();
   out_list["processes_list"]   = process_list;
   return(out_list);
}


