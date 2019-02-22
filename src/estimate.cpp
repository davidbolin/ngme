#include "estimate_util.h"
using namespace Rcpp;




// [[Rcpp::export]]
List estimateLong_cpp(Rcpp::List in_list)
{
  //**********************************
  //      basic parameter
  //**********************************

  int debug = 0;
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
  std::vector< Eigen::SparseMatrix<double,0,int> > As( nindv);
  std::vector< Eigen::VectorXd > Ys( nindv);
  std::vector< int > burnin_done(nindv);
  
  int count = 0;
  for( List::iterator it = obs_list.begin(); it != obs_list.end(); ++it ) {
    List obs_tmp = Rcpp::as<Rcpp::List>(*it);
    if(process_active)
      As[count] = Rcpp::as<Eigen::SparseMatrix<double,0,int> >(obs_tmp["A"]);
    Ys[count] = Rcpp::as<Eigen::VectorXd>(obs_tmp["Y"]);
    burnin_done[count] = 0;
    count++;
  }

  //*********************************
  //Init subsampler
  //********************************
  if(silent == 0){
    Rcpp::Rcout << " Setup subsampler\n";
  }
  subsampler* sampler = new subsampler;
  sampler->initFromList(in_list);
  
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
  std::vector<  solver* >  Solver; //solvers for Q
  Eigen::SparseMatrix<double, 0, int> Q, K;
  Eigen::VectorXd  z;

  if(process_active){
    Rcpp::List operator_list  = Rcpp::as<Rcpp::List> (in_list["operator_list"]);
    operator_list["nIter"] = nIter;
    std::string type_operator = Rcpp::as<std::string>(operator_list["type"]);

    operator_select(type_operator, &Kobj);
    Kobj->initFromList(operator_list, List::create(Rcpp::Named("use.chol") = 1));

    count = 0;
    for( List::iterator it = obs_list.begin(); it != obs_list.end(); ++it ) {
      List obs_tmp = Rcpp::as<Rcpp::List>( *it);
      Solver.push_back( new cholesky_solver );  
      Solver[count]->init(Kobj->d[count], 0, 0, 0);
      K = Eigen::SparseMatrix<double,0,int>(Kobj->Q[count]);

      Q = K.transpose();
      Q = Q * K;
      Q = Q + As[count].transpose()*As[count];
      Solver[count]->analyze(Q);
      Solver[count]->compute(Q);
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
  }else if(type_mixedEffect == "tdist") {
    mixobj   = new tdMixedEffect;
  }else if(type_mixedEffect == "NIG") {
    mixobj   = new NIGMixedEffect;
  } else {
    Rcpp::Rcout << "Wrong mixed effect distribution";
  }



  if(silent == 0){
    Rcpp::Rcout << " init";
  }
  mixobj->initFromList(mixedEffect_list);
  if(silent == 0){
    Rcpp::Rcout << " store";
  }
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
  else if(type_MeasurementError == "nsNormal")
    errObj = new nsGaussianMeasurementError;
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

  
  std::vector<Eigen::VectorXd> Vmean;
  Eigen::VectorXd count_vec;
  count_vec.setZero(nindv);
  Vmean.resize(nindv);
  if(process_active){
    for(int i = 0; i < nindv; i++ ){
      Vmean[i].setZero(Kobj->h[i].size());
    }
  }


  int par_burnin = 0;


  int npars =   mixobj->npars + errObj->npars;
  if(process_active)
    npars += process->npars + Kobj->npars;

  Eigen::MatrixXd Fisher_information(npars, npars); // If computing the fisher information
  Eigen::MatrixXd GradientVariance(npars, npars);
  Eigen::MatrixXd Fisher_information0(npars, npars);
  GradientVariance.setZero(npars, npars);
  Fisher_information.setZero(npars, npars);
  Fisher_information0.setZero(npars, npars);

  double burnin_rate = 1;
  
  //Start main loop
  for(int iter=0; iter < nIter; iter++){
    /*
    printing output
    */


    if(silent == 0){
      if(estimate_fisher == 0){
        if(debug || (burnin_rate>0) || (iter % (int)ceil(nIter/1000.) == 0)){
          Rcpp::Rcout << "--------------------------------\n";
          Rcpp::Rcout << "iter = " << iter << " (of " << nIter << "): \n";
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
      } else {
        if(nIter>1){
          Rcpp::Rcout << "iter = " << iter << " (of " << nIter << "): \n";
        } else {
          Rcpp::Rcout << "Estimated parameters:\n";
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
          Rcpp::Rcout << "Start estimation of Fisher information matrix\n";
        }

      }
    }
    if(iter < nPar_burnin)
      par_burnin = 1;


    if(debug){
      Rcpp::Rcout << "estimate::subsample \n";
      Rcpp::Rcout << "estimate::subsample_type = " << sampler->subsample_type << " \n";
    }
    int nSubsample_i = sampler->nSubsample;
    
      
    sampler->sample(iter,subsample_generator);
    if(debug)
      Rcpp::Rcout << "estimate::subsample done \n";
      
    burnin_rate = 0;
    if(process_active)
      Kobj->gradient_init(nSubsample_i,nSim);

    Eigen::VectorXd  grad_last(npars); // grad_last helper vector
    grad_last.setZero(npars);

    Eigen::MatrixXd grad_outer(sampler->nSubsample_i, npars); // outside person variation (subsampling)
    Eigen::MatrixXd grad_outer_unweighted(sampler->nSubsample_i, npars); // outside person variation (subsampling)

    Eigen::MatrixXd Vgrad_inner(npars, npars); // outside person variation (subsampling)
    Eigen::VectorXd Ebias_inner(npars);
    Ebias_inner.array() *= 0;
    Vgrad_inner.setZero(npars, npars);
    if(debug)
      Rcpp::Rcout << "estimate::start patient loop, number of patients:" << sampler->nSubsample_i << " \n";

    double pdone = 0.0;
    double next_disp = 0.00;

    //loop over patients/replicates
    for(int ilong = 0; ilong < sampler->nSubsample_i; ilong++ ){

      if(silent == 0 && estimate_fisher && nIter == 1){

        pdone = (double) ilong / (double) sampler->nSubsample_i;
        //Rcpp::Rcout << ilong << " " << pdone << " " << next_disp << "\n";
        if(pdone > next_disp){
          Rcpp::Rcout << round(100*pdone) << " % done\n";
          next_disp += 0.01;
        }

      }
      int i = sampler->longInd[ilong];
      if(debug)
        Rcpp::Rcout << "i = " << i << " "<< sampler->weight[i] << "\n";

      Eigen::SparseMatrix<double,0,int> A;
      if(process_active)
        A = As[i];
      Eigen::VectorXd  Y = Ys[i];

      /*

          If estimate_fisher is 1 then one simulate Y

      */
      if(estimate_fisher == 1){
        if(process_active){
          K = Eigen::SparseMatrix<double,0,int>(Kobj->Q[i]);
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

          process->simulate(i, z, A, K, Y,*Solver[i]);
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
      if(process_active){
        z.setZero(Kobj->h[i].size());
      } else {
        z.setZero(0);
      }
      //Gibbs sampling
      for(int ii = 0; ii < n_simulations; ii ++)
      {
        if(process_active){
          if(debug)
            Rcpp::Rcout << "simulate z\n";
          for(int j =0; j < Kobj->h[i].size(); j++)
            z[j] =  normal(random_engine);
        }

        if(debug)
          Rcpp::Rcout << "Enter Gibbs\n";
        Eigen::VectorXd res =  GibbsSampling(i,
                                             Y,
                                             A,
                                             sampleX,
                                             sampleV,
                                             process_active,
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
                           sampler->weight[i],
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
          Eigen::MatrixXd Fisher_temp  = 0.5 * grad_inner.col(count_inner)*grad_inner.col(count_inner).transpose() /  sampler->weight[i];
          Fisher_information -=  Fisher_temp + Fisher_temp.transpose(); //gives -sum g_i*g_i'
          GradientVariance += Fisher_temp + Fisher_temp.transpose();
          grad_last = grad_last_temp;
          count_inner++;
        }
      }

      // TODO::redo and move the above note i is the same always, and
      // we need both weighted unweighted
      // grad_inner.array() /= weight[i]; // dont use weight here(?)
      // do grad_outer_unweighted
      Eigen::VectorXd Mgrad_inner = grad_inner.rowwise().mean(); //gives E(g)
      grad_outer.row(ilong) = Mgrad_inner;
      grad_outer_unweighted.row(ilong) = Mgrad_inner;
      grad_outer_unweighted.row(ilong) /= sampler->weight[i];
      Eigen::MatrixXd Fisher_add  = 0.5 * nSim  * (Mgrad_inner/sampler->weight[i]) * Mgrad_inner.transpose();
      Fisher_information  +=  Fisher_add + Fisher_add.transpose() ; //add N*E(g)*E(g)'
      GradientVariance    -=  Fisher_add + Fisher_add.transpose();
      Fisher_information0 +=  Fisher_add + Fisher_add.transpose();

      Eigen::MatrixXd centered = grad_inner.colwise() - Mgrad_inner;
      Ebias_inner.array() += centered.col(nSim-1).array();
      Ebias_inner.array() -= centered.col(0).array();
      Eigen::MatrixXd cov = (centered * centered.transpose()) /  (sampler->weight[i]* double(grad_inner.cols() - 1));
      Vgrad_inner.array() += cov.array();
      if(process_active)
        Vmean[i] += process->Vs[i];
      count_vec[i] += 1;
    }
    // update weights given the gradient
    // change here to unweighted gradient!
    if(sampler->subsample_type == 3){
      Eigen::VectorXd W = gradientWeight(grad_outer_unweighted);
      for( int id = 0; id < sampler->nSubsample_i; id++)
        sampler->p_inv[sampler->longInd[id]] = W[id];
    }
   
    if((estimate_fisher == 0) && (silent == 0) && burnin_rate>0){
      Rcpp::Rcout << "Burnin percentage " << burnin_rate/sampler->nSubsample_i << "\n";
    }


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
    out_list["FisherMatrix0"]    = Fisher_information0 / (nIter * nSim);

    out_list["GradientVariance"] = GradientVariance / (nIter * nSim) ;

  }

  out_list["pSubsample"]       = sampler->pSubsample;
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
  
  //free solver
  for( std::vector<solver*>::iterator i = Solver.begin(); i != Solver.end(); ++i )
  {
    delete *i;
  }
  Solver.clear();
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
