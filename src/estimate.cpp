#include "estimate_util.h"
#include "MatrixAlgebra.h"
using namespace Rcpp;

/*
  Caculating the the variance of the gradient conditioning only on the Varaiance components
  For fixed and random effects
*/
Eigen::MatrixXd var_term_calc(
                                Eigen::VectorXd&  Y,
                                int i,
                                Eigen::SparseMatrix<double, 0, int>  A,
                                MixedEffect& mixobj,
                                operatorMatrix& Kobj,
                                MeasurementError& errObj,
                                Process& process,
                                int process_active,
                                int calc_mean,
                                double weight,
                                int debug
                              )

{
  if(debug)
    Rcpp::Rcout << "in var_term_calc\n";
  int n_r =  0;
  if(mixobj.Br.size() > 0)
    n_r =  mixobj.Br[i].cols();

  int n_f = 0;
  if(mixobj.Bf.size() > 0)
    n_f = mixobj.Bf[i].cols();
  int n_p = 0;
  if(process_active)
    n_p = A.cols();
  
  if(debug)
    Rcpp::Rcout << "n_r = " << n_r << ", n_f = " << n_f << ", n_p = " << n_p << "\n";

  //###############################################
  //#create joint observation matrix Ajoint = [B A]
  //################################################
  
  if(debug)
    Rcpp::Rcout << "create Ajoint \n";
  
  Eigen::SparseMatrix<double, 0, int> Ajoint;
  Ajoint.resize(Y.size(), n_p + n_r);

  if(n_r>0){
    Eigen::SparseMatrix<double, 0, int> B = full2sparse(mixobj.Br[i]);
    setSparseBlock(&Ajoint, 0, 0, B);
  }

  if(process_active)
    setSparseBlock(&Ajoint, 0, n_r, A);

  //###########################################################
  //#create joint operator matrix Kjoint = [sqrt(Sigma_u) K]  #
  //###########################################################
  if(debug)
    Rcpp::Rcout << "create Kjoint \n";
  
  Eigen::SparseMatrix<double, 0, int> Kjoint;
  Kjoint.resize(n_p + n_r,n_p + n_r);
  Eigen::VectorXd mu_prior(n_r + n_p);
  Eigen::MatrixXd iSigma;
  if(n_r > 0){
    iSigma = mixobj.Sigma.inverse();
    Eigen::MatrixXd L = iSigma.llt().matrixL().transpose();
    Eigen::SparseMatrix<double,0,int> Sigma_root = full2sparse(L);
    setSparseBlock(&Kjoint,0,0,Sigma_root);
  }
  
  if(n_p > 0){
    Eigen::SparseMatrix<double, 0, int>  K;
    if(Kobj.nop==1){
      K= Eigen::SparseMatrix<double,0,int>(Kobj.Q[0]);
    } else {
      K= Eigen::SparseMatrix<double,0,int>(Kobj.Q[i]);  
    }
    
    setSparseBlock(&Kjoint, n_r, n_r,K);
    mu_prior.tail(n_p)   = process.get_mean_prior(i, K);
  }

  //################################################
  //#Compute residual Y - Xbeta - mean
  //###############################################
  if(debug)
    Rcpp::Rcout << "compute residual \n";
  
  Eigen::VectorXd  res = Y;
  mixobj.remove_cov(i,  res);
  errObj.remove_asym(i, res);

  //################################################
  //#Build b and Q for the joint distribution
  //###############################################
  if(debug)
    Rcpp::Rcout << "build b and Q \n";
  
  Eigen::VectorXd Vinv_joint(n_r + n_p);
  
  if(n_r > 0){
    Eigen::VectorXd V_r;
    V_r.setOnes(n_r);
    V_r /=  mixobj.V(i);
    Vinv_joint.head(n_r) = V_r;
    mu_prior.head(n_r)   = mixobj.get_mean_prior(i);
  }

  if(n_p > 0)
    Vinv_joint.tail(n_p) = process.Vs[i].cwiseInverse();
    
  Eigen::VectorXd Q_e = errObj.getSigma(i, Y.size()).cwiseInverse();
  Eigen::SparseMatrix<double,0,int> Q_hat = Kjoint.transpose()*Vinv_joint.asDiagonal()*Kjoint;
  Eigen::VectorXd b_hat = Q_hat*mu_prior;
  Q_hat +=  Ajoint.transpose()*Q_e.asDiagonal()*Ajoint;
  b_hat += Ajoint.transpose()*Q_e.asDiagonal()*res;
  Eigen::SimplicialLDLT< Eigen::SparseMatrix<double,0,int> > Qsolver;
  Qsolver.compute(Q_hat);    
  Eigen::VectorXd mu_hat = Qsolver.solve(b_hat);
//########################################
//#computation of the variance
//########################################
if(debug)
  Rcpp::Rcout << "compute variance component \n";

 Eigen::MatrixXd X(Y.size(), mixobj.nfr);
 Eigen::MatrixXd Br;
 if(n_r > 0){
  if(mixobj.noise == "Normal"){
      Br = mixobj.Br[i];
    }else{
      Br.resize(mixobj.Br[i].rows(), 2*mixobj.Br[i].cols());
      Eigen::MatrixXd Br_mu = (mixobj.V(i)-1.) * mixobj.Br[i];
      Br << mixobj.Br[i], Br_mu;
    }
   
    //Rcpp::Rcout << "E = \n" << E << "\n";
 }

  if(n_r == 0){
    X = mixobj.Bf[i];
  }else if(n_f == 0){
      X = Br;
  }else{
    X << mixobj.Bf[i], Br;
  }

  Eigen::MatrixXd QX = Q_e.asDiagonal()*X;
  Eigen::MatrixXd Results = - X.transpose() * QX; 
  
  Eigen::MatrixXd  XtQ = QX.transpose();
  Eigen::MatrixXd AtQX = Ajoint.transpose()*QX;
  Eigen::MatrixXd XtQA = AtQX.transpose();
  Eigen::MatrixXd Res;

  //Result_full
  Eigen::MatrixXd ResultsFull;
  Eigen::VectorXd gradFull;
  gradFull.setZero(mixobj.npars);
  ResultsFull.setZero(mixobj.npars, mixobj.npars);

  Eigen::VectorXd grad  = XtQ * (res - Ajoint * mu_hat);
  gradFull.head(grad.size()) = grad;
  int n_sigma = 0;
  if(n_r > 0){
     /*
      Gradient of Sigma and Fisher
    */
      Eigen::VectorXd  dSigma_vech;
      Eigen::MatrixXd PreCond = 0.5 * mixobj.Dd.transpose() * mixobj.iSkroniS;
      Eigen::VectorXd dSigma_vech1, dSigma_vech2;
      Eigen::MatrixXd E;
      E.setIdentity(b_hat.size(), n_r);
      Eigen::MatrixXd Sigma_Random = Qsolver.solve(E);
      Eigen::MatrixXd  EUUt = NormalOuterExpectation(Sigma_Random.topLeftCorner(n_r,n_r),
                                                     mu_hat.head(n_r),
                                                     mixobj.get_mean_prior(i));
      EUUt.array() /= mixobj.V(i);
      dSigma_vech1 = vec(EUUt);
      Eigen::VectorXd vSigma = vec(mixobj.Sigma);
      dSigma_vech2 =  vSigma;
      dSigma_vech = PreCond * (dSigma_vech1 - dSigma_vech2);

      Eigen::MatrixXd EgradgradT, Fisher_Sigma_temp;
      EgradgradT =  dSigma_vech2 * dSigma_vech2.transpose();
      EgradgradT -= dSigma_vech1 * dSigma_vech2.transpose();
      EgradgradT -= dSigma_vech2 * dSigma_vech1.transpose();
      Eigen::VectorXd mean_temp = mu_hat.head(n_r) - mixobj.get_mean_prior(i);
      EgradgradT += NormalouterKron(Sigma_Random.topLeftCorner(n_r,n_r),
                                      mean_temp)/pow(mixobj.V(i),2);
      Fisher_Sigma_temp = PreCond * EgradgradT * PreCond.transpose();
      n_sigma = Fisher_Sigma_temp.cols();
      ResultsFull.block(mixobj.nfr, mixobj.nfr, n_sigma, n_sigma) += Fisher_Sigma_temp;
      ResultsFull.block(mixobj.nfr, mixobj.nfr, n_sigma, n_sigma) -= dSigma_vech * dSigma_vech.transpose();
      Eigen::MatrixXd Hessian = HessianSigma(EUUt,
                                             iSigma,
                                             mixobj.Sigma,
                                             1);
      ResultsFull.block(mixobj.nfr, mixobj.nfr, n_sigma, n_sigma) +=  mixobj.Dd.transpose() * Hessian * mixobj.Dd;
      gradFull.segment(grad.size(), dSigma_vech.size()) = dSigma_vech;

      /*
        Fisher [beta,mu] * [Sigma]^T
      */
      Eigen::VectorXd aTilde = XtQ * (res + Br * mixobj.get_mean_prior(i));
      Eigen::VectorXd mu_hat_tilde = mu_hat;
      mu_hat_tilde.head(n_r).array() -= mixobj.get_mean_prior(i).array();
      // (-a + AX_i) * vec(Sigma)^T * H^T
      ResultsFull.block(0, mixobj.nfr, mixobj.nfr, n_sigma) -= (grad * vSigma.transpose())*PreCond.transpose();

      ResultsFull.block(0, mixobj.nfr, mixobj.nfr, n_sigma) += aTilde * dSigma_vech1.transpose() * PreCond.transpose();
      Eigen::MatrixXd EUUxT = Normalxxty(Sigma_Random.topLeftCorner(n_r,n_r),
                                         Sigma_Random.transpose(),
                                         mu_hat_tilde.head(n_r), 
                                         mu_hat_tilde); 
      ResultsFull.block(0, mixobj.nfr, mixobj.nfr, n_sigma) -= (XtQ * (Ajoint * EUUxT.transpose() ))* PreCond.transpose();
      ResultsFull.block(0, mixobj.nfr, mixobj.nfr, n_sigma) -=  grad * dSigma_vech.transpose();
      ResultsFull.block(mixobj.nfr, 0, n_sigma, mixobj.nfr) += ResultsFull.block(0, mixobj.nfr, mixobj.nfr, n_sigma).transpose();
     
    Eigen::VectorXd grad2;
    if(mixobj.noise != "Normal"){
      Eigen::VectorXd grad2_temp = iSigma * mu_hat.head(n_r);
      grad2.resize(2*n_r);
      grad2_temp -= (mixobj.V(i)-1)* iSigma * mixobj.mu;
      grad2_temp /= mixobj.V(i);
      grad2.head(n_r) = grad2_temp;
      grad2.tail(n_r) = (mixobj.V(i) -1.)*grad2_temp;
    }else{
      grad2.resize(n_r);
      grad2 = iSigma * mu_hat.head(n_r);
    }
    grad2.array() *= weight;
      mixobj.add_gradient2(grad2);
    //return Results;
  }else{
    grad *= weight;  
      mixobj.add_gradient(grad);
  }

   
  
  mixobj.updateFisher(i, ResultsFull, gradFull);
  if(n_r >0 ){
    gradFull.array() *=  weight;
    mixobj.add_gradient(gradFull);
  }
  Eigen::MatrixXd tmp = Qsolver.solve(AtQX);
  
  Results += XtQA*tmp;
  ResultsFull.topLeftCorner(mixobj.nfr, mixobj.nfr) += Results;
  Res = weight * Results;
  mixobj.get_Hessian(Res);
  return ResultsFull;
}


// [[Rcpp::export]]
List estimateLong_cpp(Rcpp::List in_list)
{
  //**********************************
  //      basic parameter
  //**********************************

  int debug = 0;
  int calc_grad=0;
  int nIter      = Rcpp::as< int > (in_list["nIter"]);
  int nSim       = Rcpp::as< int > (in_list["nSim"]);
  int nBurnin    = Rcpp::as< int > (in_list["nBurnin"] );
  int nBurnin_base    = Rcpp::as< int > (in_list["nBurnin_base"] );
  int nBurnin_learningrate = nBurnin;
  int iter_start = 0;
  if(in_list.containsElementNamed("iter_start"))
    iter_start    = Rcpp::as< int > (in_list["iter_start"] );
  
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
      if(Kobj->nop == 1){
        Solver[count]->init(Kobj->d[0], 0, 0, 0);
        K = Eigen::SparseMatrix<double,0,int>(Kobj->Q[0]);  
      } else {
        Solver[count]->init(Kobj->d[count], 0, 0, 0);
        K = Eigen::SparseMatrix<double,0,int>(Kobj->Q[count]);
      }
      
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
  mixobj->set_calc_grad(calc_grad);


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
  int nfr = mixobj->get_nf() + mixobj->get_nr(); //number of covariates
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
    type_processes= Rcpp::as<std::string> (processes_list["noise"]);

    if(debug)
      Rcpp::Rcout << "process name = " << type_processes << "\n";
    if (type_processes == "Normal"){
       process  = new GaussianProcess;     
    }else if(type_processes == "MultiGH" ){
      process  = new MGHProcess;
    }else{
      process  = new GHProcess;
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
      if(Kobj->nop==1){
        Vmean[i].setZero(Kobj->h[0].size());
      } else {
        Vmean[i].setZero(Kobj->h[i].size());  
      }
      
    }
  }


  int par_burnin = 0;


  int npars =   mixobj->npars + errObj->npars;
  if(process_active)
    npars += process->npars + Kobj->npars;

  Eigen::MatrixXd Fisher_information(npars, npars); // If computing the fisher information
  Eigen::MatrixXd GradientVariance(npars, npars);
  Eigen::MatrixXd Fisher_information0(npars, npars);
  Eigen::VectorXd grad_mix(mixobj->npars );
  grad_mix.array() *= 0;
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


    Eigen::MatrixXd Vmf(0,0); // variance of fixed and mixed effect
    Eigen::VectorXd grad_mix_temp(0); // variance of fixed and mixed effect
    if(mixobj->get_nf() + mixobj->get_nr()>0){
      Vmf.setZero( mixobj->npars, mixobj->npars);
      grad_mix_temp.setZero(mixobj->npars);
      }
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
          if(Kobj->nop == 1){
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
        if(Kobj->nop==1){
          z.setZero(Kobj->h[0].size());
        } else {
          z.setZero(Kobj->h[i].size());  
        }
      } else {
        z.setZero(0);
      }
      //Gibbs sampling
      for(int ii = 0; ii < n_simulations; ii ++)
      {
        if(process_active){
          if(debug)
            Rcpp::Rcout << "simulate z\n";
          for(int j =0; j < z.size(); j++)
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
          if(calc_grad == 0 && estimate_fisher == 0){
            var_term_calc(Y,
                          i,
                          A,
                          *mixobj,
                          *Kobj,
                          *errObj,
                          *process,
                          process_active,
                          1,
                          sampler->weight[i],
                          debug);
          }
          if(estimate_fisher && nfr > 0){
            if(debug)
              Rcpp::Rcout << "enter fisher var_term_calc \n";
            
              Vmf +=          var_term_calc(Y,
                                            i,
                                            A,
                                            *mixobj,
                                            *Kobj,
                                            *errObj,
                                            *process,
                                            process_active,
                                            0,
                                            sampler->weight[i],
                                            debug);
              if(debug)
                Rcpp::Rcout << "var_term_calc done\n";
              //if(i==0)
              //Rcpp::Rcout << "Vmf = " << Vmf << "\n";
              grad_mix_temp += mixobj->get_gradient();

          }
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
          count_inner++;

        }
      }
      if(nfr> 0){
          Fisher_information.topLeftCorner(mixobj->npars, mixobj->npars) -= Vmf * (sampler->weight[i] );
          grad_mix += grad_mix_temp;
          Vmf *= 0; 
        }
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
      //double stepsize = step0/pow(iter + 1 + iter_start, alpha);
      double stepsize = step0/pow(iter + 1 + iter_start, alpha*(iter+iter_start)/(nIter+iter_start));
      
      double polyak_rate_temp = polyak_rate;
      double learning_rate_temp  =learning_rate;
      if(polyak_rate == 0)
        polyak_rate_temp = 1./ (iter + 1);
      if(iter < nBurnin_learningrate)
        learning_rate_temp = 0;
      if(debug)
        Rcpp::Rcout << "polyak_rate_temp = " << polyak_rate_temp <<"\n";
      mixobj->step_theta(stepsize,  learning_rate_temp, polyak_rate_temp, iter + iter_start);
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
  Fisher_information.array()  *= Ys.size();
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
    out_list["gradient"]         = grad_mix/(nIter * nSim);
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
