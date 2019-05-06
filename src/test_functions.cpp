#include <Rcpp.h>
#include <RcppEigen.h>
#include "sample.h"
#include "solver.h"
#include "GIG.h"
#include "GHmisc.h"
#include "MixedEffect.h"
#include "MatrixAlgebra.h"
#include "operatorMatrix.h"
#include "operator_helper.h"
#include "latentprocess.h"
#include "measError.h"
#include "estimate_util.h"
//[[Rcpp::export]]
Eigen::MatrixXi getDuplicateM(const int n)
{

  return(duplicatematrix(n));
}


// [[Rcpp::export]]
Rcpp::List test_d2_process(Rcpp::List Y,
                           Rcpp::List process_list,
                           Rcpp::List operator_list)
{
  operatorMatrix* Kobj;
  std::string type_operator = Rcpp::as<std::string>(operator_list["type"]);
  operator_select(type_operator, &Kobj);
  Kobj->initFromList(operator_list, Rcpp::List::create(Rcpp::Named("use.chol") = 1));
  Process *process = NULL;
  process  = new GHProcess;
  process->initFromList(process_list, Kobj->h);
  Eigen::SparseMatrix<double,0,int> A = Rcpp::as<Eigen::SparseMatrix<double,0,int> >(Y["A"]);
  Eigen::VectorXd Y_ =Rcpp::as<Eigen::VectorXd  >(Y["Y"]);
  Eigen::MatrixXd d2 = process->d2Given(0, Kobj->Q[0], A, Y_,
                                         1.,
                                         1.,
                                         1.);

  process_list["d2"] = d2;
  delete process;
  delete Kobj;
  return(process_list);
}

// [[Rcpp::export]]
Rcpp::List  test_sampling_NIG(Rcpp::List mixedEffect_list,
                       Rcpp::List meas_list,
                       int nsamples)
{
  NIGMixedEffect *mixobj = NULL;
  mixobj   = new NIGMixedEffect;
  mixobj->initFromList(mixedEffect_list);
  Rcpp::List Y =  Rcpp::as<  Rcpp::List    >(meas_list["Y"]);
  double sigma = Rcpp::as< double >(meas_list["sigma_eps"]);
  
  Eigen::MatrixXd U_MALA, U_Gibbs;
  U_MALA.setZero(nsamples, mixobj->U.rows());
  for(int i =0; i < nsamples; i++){
    //mixobj->sampleU_MALA(0, Rcpp::as< Eigen::VectorXd>(Y[0]), 2*log(sigma));
    U_MALA.row(i) = mixobj->U.col(0);  
  }
  
  U_Gibbs.setZero(nsamples, mixobj->U.rows());
  for(int i =0; i < nsamples; i++){
    mixobj->sampleU( 0, Rcpp::as< Eigen::VectorXd>(Y[0]),  2*log(sigma));
    U_Gibbs.row(i) = mixobj->U.col(0);  
  }
  
  Rcpp::List out;
  //out["acc_MALA"] = ((double)mixobj->accept_MALA ) /((double) mixobj->count_MALA );
  out["U_MALA"]   = U_MALA;
  out["U_Gibbs"]  = U_Gibbs;
  free(mixobj);
  
  return(out);
}

// [[Rcpp::export]]
double test_logf_NIG(const Eigen::VectorXd & U,
                     const Eigen::VectorXd & mu,
                     const Eigen::VectorXd & delta,
                     const Eigen::MatrixXd & iSigma,
                     const double nu )
{
  return(logNIG(U,
                mu,
                delta,
                iSigma,
                nu ));
  
  
}


// [[Rcpp::export]]
double test_logf_GH(const Eigen::VectorXd & U,
                     const Eigen::VectorXd & mu,
                     const Eigen::VectorXd & delta,
                     const Eigen::MatrixXd & iSigma,
                     const double p,
                     const double a,
                     const double b)
{
  return(logGH(U,
                mu,
                delta,
                iSigma,
                p,
                a,
                b));
  
  
}

// [[Rcpp::export]]
Rcpp::List test_dU_EiV(
                 const Eigen::VectorXd & U,
                 const Eigen::MatrixXd & Sigma,
                 const Eigen::VectorXd & delta,
                 const Eigen::VectorXd & mu,
                 const double p_GIG,
                 const double a_GIG,
                 const double b_GIG,
                 const Eigen::VectorXd & res,
                 const Eigen::MatrixXd & Q_noise,
                 const Eigen::MatrixXd & B){
  
  Eigen::VectorXd dU;
  Eigen::MatrixXd ddU;
  dU_ddU_GH(dU,
             ddU,
             U,
             Sigma.inverse(),
             delta,
              mu,
            -0.5,
             a_GIG,
             b_GIG,
             res,
             Q_noise,
             B);
  
  Rcpp::List out;
  out["dU"] = dU;
  out["ddU"] = ddU;
  return(out);
  
}

// [[Rcpp::export]]
double test_db_EiV_GIG(double p, double a, double b){
  return( db_EiV_GIG(p,a,b));
}

// [[Rcpp::export]]
double test_EiV_NGIG( Eigen::VectorXd & U,
                      Eigen::MatrixXd & Sigma,
                      Eigen::VectorXd & delta,
                      Eigen::VectorXd & mu,
                      double p,
                      double a,
                      double b){
  return EiV_NGIG(U, Sigma, delta, mu, p, a, b);
}

// [[Rcpp::export]]
double test_EiV_GIG(double p, double a, double b) {
	return EiV_GIG(p, a, b);
}
// [[Rcpp::export]]
double test_PreDiagsolver(Rcpp::List in_list)
{
	Eigen::SparseMatrix<double, 0, int> Q = Rcpp::as<Eigen::SparseMatrix<double,0,int> >(in_list["Q"]);
	Eigen::VectorXd z = Rcpp::as<Eigen::VectorXd > (in_list["z"]);
	Eigen::VectorXd b = Rcpp::as<Eigen::VectorXd > (in_list["b"]);
	cholesky_solver Solver;
	Solver.init(Q.rows(), 0, 0, 0);
	Solver.analyze(Q);
  	Solver.compute(Q);
  	
    Eigen::VectorXd res = Solver.rMVN(b, z);
    
    Eigen::SparseMatrix<double, 0, int> Q2(Q);
    Eigen::VectorXd D_12 = Q.diagonal().cwiseInverse().cwiseSqrt();
    Q2 = D_12.asDiagonal() * Q * D_12.asDiagonal();
    Solver.analyze(Q2);
    Solver.compute(Q2);
    Eigen::VectorXd b_ = b.cwiseProduct(D_12);
    Eigen::VectorXd res2 = Solver.rMVN(b_, z);
    res2 = res2.cwiseProduct(D_12);
    return((res - res2).sum());
}

// [[Rcpp::export]]
Eigen::VectorXi   sampleR(int n, Eigen::VectorXd w_in)
{
	Eigen::VectorXd w;
  	w = w_in;
  	w /= w.sum();
  	return(ProbSampleNoReplace(n, w));
}


// [[Rcpp::export]]
Rcpp::List  sample_internalR(int n, 
                                   Eigen::VectorXd p_in,
                                   Eigen::VectorXd selected_in,
                                   Eigen::VectorXd w_in)
{
  std::vector<int> ans;
  std::vector<int> selected(p_in.size());
  for(int i = 0; i < p_in.size(); i++)
    selected[i] = (int) (selected_in[i]);
  
  
  Eigen::VectorXd p;
  p = p_in;
  p /= p.sum();
  Rcpp::List out;
  poissonSampling_internal(n, p, w_in, ans, selected);
  
  out["ans"]   = ans;
  out["w_in"]  = w_in;
  out["selected"]  = selected;
  return(out);
}


/*
  simple optimazation of mixed effect object
  niter      - number of set in conjuagete gradient
  Y          - observtions
  mixed_list - conatins init ofr mixed_list

*/
// [[Rcpp::export]]
Rcpp::List test_mixed(int niter,
                      Rcpp::List Y,
                      Rcpp::List mixed_list,
                      Rcpp::List error_list)
{
  int nindv = Y.size(); 
  std::vector< Eigen::VectorXd > Ys( nindv);
  for(int i=0; i < nindv; i++)
    Ys[i] = Rcpp::as<Eigen::VectorXd>(Y[i]);
  
  
  MixedEffect *mixobj = NULL;
  std::string type_mixed = Rcpp::as <std::string> (mixed_list["name"]);
  if(type_mixed == "Normal"){
      mixobj = new NormalMixedEffect;
    }else if(type_mixed == "tdist") {
      mixobj   = new tdMixedEffect;
    }else if(type_mixed == "NIG") {
      mixobj   = new NIGMixedEffect;
    } else {
      Rcpp::Rcout << "Wrong mixed effect distribution";
    }
    mixobj->initFromList(mixed_list);



  MeasurementError *errObj;
  std::string type_me = Rcpp::as <std::string> (error_list["name"]);
  if(type_me == "Normal")
    errObj = new GaussianMeasurementError;
  else if(type_me == "nsNormal")
    errObj = new nsGaussianMeasurementError;
  errObj->initFromList(error_list);
  
  mixobj->setupStoreTracj(niter);
  errObj->setupStoreTracj(niter);
  for(int iter=0; iter < niter; iter++){

    for(int i=0; i < nindv; i++){
      Eigen::VectorXd  res = Y[i];
      mixobj->remove_cov(i, res);

      //mixobj->remove_inter(i, res);
      //mixobj->add_inter(i, res);
      Eigen::VectorXd sigmas; 

      if(errObj->nsSigma)
        sigmas = errObj->sigmas[i];
      if(errObj->nsSigma>0){
        mixobj->sampleU2( i, res, sigmas.array().pow(2).cwiseInverse(),2 * log(errObj->sigma));
      }else{
        mixobj->sampleU( i, res, 2 * log(errObj->sigma));
      }


      if(type_me == "Normal"){
        mixobj->gradient(i, res, 2 * log(errObj->sigma), 1., 1);
      }else{
      mixobj->gradient2(i,
                      res,
                      errObj->Vs[i].cwiseInverse(),
                      sigmas,
                      2 * log(errObj->sigma),
                      1.,
                      1., //w 
                      1.,
                      errObj->nsSigma);
      }

      mixobj->remove_inter(i, res);
      
      errObj->gradient(i, res, 1.);
    }

    mixobj->step_theta(0.5, 0, -1);

    errObj->step_theta(0.5, 0, -1);
    mixobj->clear_gradient();
    errObj->clear_gradient();
  }
  
  Rcpp::List out_list;
  Rcpp::List mixobj_list       = mixobj->toList();
  out_list["mixedEffect_list"] = mixobj_list;
  Rcpp::List errobj_list            = errObj->toList();
  out_list["measurementError_list"] = errobj_list;

  return(out_list);
}


/*
  simple optimazation of error effect object espically assymetric eerror
  niter      - number of set in conjuagete gradient
  Y          - observtions
  mixed_list - conatins init ofr mixed_list

*/
// [[Rcpp::export]]
Rcpp::List test_error(int niter,
                      Rcpp::List Y,
                      Rcpp::List error_list)
{
 int nindv = Y.size(); 
  std::vector< Eigen::VectorXd > Ys( nindv);
  for(int i=0; i < nindv; i++)
    Ys[i] = Rcpp::as<Eigen::VectorXd>(Y[i]);
  

  MeasurementError *errObj;
  std::string type_me = Rcpp::as <std::string> (error_list["name"]);
  if(type_me == "Normal")
    errObj = new GaussianMeasurementError;
  else if(type_me == "nsNormal")
    errObj = new nsGaussianMeasurementError;
  else if(type_me == "NIG")
    errObj = new NIGMeasurementError;
  
  errObj->initFromList(error_list);
  
  errObj->setupStoreTracj(niter);
  for(int iter=0; iter < niter; iter++){

    for(int i=0; i < nindv; i++){
      Eigen::VectorXd  res = Y[i];
     
      Eigen::VectorXd sigmas; 

      if(errObj->nsSigma)
        sigmas = errObj->sigmas[i];

      errObj->sampleV(i, res);
      errObj->remove_asym(i, res);
      errObj->gradient(i, res, 1.);
    }
    
    errObj->step_theta(0.5, 0, -1);
    errObj->clear_gradient();
  }
  
  Rcpp::List out_list;
  Rcpp::List errobj_list            = errObj->toList();
  out_list["measurementError_list"] = errobj_list;
  return(out_list);

}


/*
  simple optimazation of mixed effect object
  niter      - number of set in conjuagete gradient
  Y          - observtions
  mixed_list - conatins init ofr mixed_list

*/
// [[Rcpp::export]]
Rcpp::List test_mixed_Fisher(int niter,
                      Rcpp::List Y,
                      Rcpp::List mixed_list,
                      Rcpp::List error_list)
{
  int nindv = Y.size(); 
  std::vector< Eigen::VectorXd > Ys( nindv);
  for(int i=0; i < nindv; i++)
    Ys[i] = Rcpp::as<Eigen::VectorXd>(Y[i]);
  
  
  MixedEffect *mixobj = NULL;
  std::string type_mixed = Rcpp::as <std::string> (mixed_list["name"]);
  if(type_mixed == "Normal"){
      mixobj = new NormalMixedEffect;
    }else if(type_mixed == "tdist") {
      mixobj   = new tdMixedEffect;
    }else if(type_mixed == "NIG") {
      mixobj   = new NIGMixedEffect;
    } else {
      Rcpp::Rcout << "Wrong mixed effect distribution";
    }
    mixobj->initFromList(mixed_list);



  MeasurementError *errObj;
  std::string type_me = Rcpp::as <std::string> (error_list["name"]);
  if(type_me == "Normal")
    errObj = new GaussianMeasurementError;
  else if(type_me == "nsNormal")
    errObj = new nsGaussianMeasurementError;
  errObj->initFromList(error_list);
  
  mixobj->setupStoreTracj(niter);
  errObj->setupStoreTracj(niter);
  Eigen::MatrixXd d2Given;
  Eigen::MatrixXd VVt;
  Eigen::VectorXd V;
  Eigen::VectorXd Vtemp;
  Eigen::VectorXd Vtemp_old;
  d2Given.setZero(mixobj->npars, mixobj->npars);
  VVt.setZero(mixobj->npars, mixobj->npars);
  V.setZero(mixobj->npars);
  Vtemp.setZero(mixobj->npars);
  for(int iter=0; iter < niter; iter++){

    Vtemp_old= Vtemp;
    Vtemp.setZero(mixobj->npars);
    for(int i=0; i < nindv; i++){
      Eigen::VectorXd  res = Y[i];
      mixobj->remove_cov(i, res);

      //mixobj->remove_inter(i, res);
      //mixobj->add_inter(i, res);
      Eigen::VectorXd sigmas; 
      Eigen::MatrixXd d2Given_temp;
      if(errObj->nsSigma)
        sigmas = errObj->sigmas[i];
      if(errObj->nsSigma>0){
        mixobj->sampleU2( i, res, sigmas.array().pow(2).cwiseInverse(),2 * log(errObj->sigma));
      }else{
        mixobj->sampleU( i, res, 2 * log(errObj->sigma));
      }

      if(type_me == "Normal"){
        mixobj->gradient(i, res, 2 * log(errObj->sigma), 1., 0);
        d2Given_temp = mixobj->d2Given(i,
                                  res,
                                  2 * log(errObj->sigma),
                                  1.);
      }else{
      mixobj->gradient2(i,
                      res,
                      errObj->Vs[i].cwiseInverse(),
                      sigmas,
                      2 * log(errObj->sigma),
                      1.,
                      1., //w 
                      0,
                      errObj->nsSigma);

        d2Given_temp = mixobj->d2Given2(i,
                                   res,
                                   errObj->Vs[i].cwiseInverse(),
                                   2 * log(errObj->sigma),
                                   1.,
                                   0);

      }
      mixobj->remove_inter(i, res);
      d2Given += d2Given_temp.topLeftCorner(mixobj->npars,mixobj->npars);
      Vtemp     += mixobj->gradientVec;
      errObj->gradient(i, res, 1.);
    }
    VVt += Vtemp * Vtemp.transpose();
    V   += Vtemp;
    //mixobj->step_theta(0.5, 0, -1);
    //errObj->step_theta(0.5, 0, -1);
    mixobj->clear_gradient();
    errObj->clear_gradient();
  }
  Rcpp::List out_list;
  Rcpp::List mixobj_list       = mixobj->toList();
  out_list["mixedEffect_list"] = mixobj_list;
  Rcpp::List errobj_list            = errObj->toList();
  out_list["measurementError_list"] = errobj_list;
  out_list["V"] = V/niter;
  out_list["VVt"] = VVt/niter;
  out_list["d2Given"] = d2Given/niter;

  return(out_list);
}


// [[Rcpp::export]]
Rcpp::List test_process(int niter,
                      Rcpp::List U,
                      Rcpp::List process_list)
{

  std::mt19937 random_engine;
  gig rgig;
  rgig.seed(random_engine());
  gig *rgig_pointer = &rgig;
  int sampleV = 1;
  if(process_list.containsElementNamed("sampleV" ))
    sampleV = Rcpp::as<int>(process_list["sampleV"]);
  
  Rcpp::List h_list  = Rcpp::as<Rcpp::List> (process_list["h"]);
  
  std::vector<Eigen::VectorXd > h;
  h.resize(h_list.length());
  for(int i = 0; i < h_list.length(); i++)
    h[i] = Rcpp::as<Eigen::VectorXd> (h_list[i]);


  Process *process = NULL;
  process = new GHProcess;
  process->initFromList(process_list, h);
  process->setupStoreTracj(niter);
  Eigen::MatrixXd I = MatrixXd::Identity(h[0].size(), h[0].size());

  SparseMatrix<double> spI = I.sparseView();
  for(int j = 0; j< niter; j++){
    for(int i = 0; i < h_list.length(); i++){
      if(sampleV)
        process->sample_V(i ,
                          rgig,
                          spI);
      process->gradient(i, spI, spI, h[0], 1., 1., 1.);
    

      //process->gradient_v2(i , spI, spI, h[0], 1., h[0], 1., 1., 1.);
    }
    process->step_theta( 0.95, 0., -1, 0);
    process->clear_gradient();
  }

  Rcpp::List out_list;
  Rcpp::List process_list_out  = process->toList();
  out_list["process_list"]     = process_list_out;
  return(out_list);

}
// [[Rcpp::export]]
Rcpp::List test_Mprocess(int niter,
                      Rcpp::List U,
                      Rcpp::List process_list)
{

  std::mt19937 random_engine;
  gig rgig;
  rgig.seed(random_engine());
  gig *rgig_pointer = &rgig;
  int sampleV = 1;
  if(process_list.containsElementNamed("sampleV" ))
    sampleV = Rcpp::as<int>(process_list["sampleV"]);
 
  Rcpp::List h_list  = Rcpp::as<Rcpp::List> (process_list["h"]);
  
  std::vector<Eigen::VectorXd > h;
  h.resize(h_list.length());
  for(int i = 0; i < h_list.length(); i++)
    h[i] = Rcpp::as<Eigen::VectorXd> (h_list[i]);

  Process *process = NULL;
  process = new MGHProcess;
  process->initFromList(process_list, h);
  process->setupStoreTracj(niter);
  Eigen::MatrixXd I = MatrixXd::Identity(h[0].size(), h[0].size());

  SparseMatrix<double> spI = I.sparseView();
  for(int j = 0; j< niter; j++){
    for(int i = 0; i < h_list.length(); i++){
      if(sampleV)
        process->sample_V(i ,
                          *rgig_pointer,
                          spI);
      process->gradient(i, spI, spI, h[0], 1., 1., 1.);
      //process->gradient_v2(i , spI, spI, h[0], 1., h[0], 1., 1., 1.);
    }

    process->step_theta( 0.95, 0., -1, 0);
    process->clear_gradient();
  }

  Rcpp::List out_list;
  Rcpp::List process_list_out  = process->toList();
  out_list["process_list"]     = process_list_out;
  return(out_list);

}



// [[Rcpp::export]]
Rcpp::List test_2Doperator(Rcpp::List process_list,
                           Rcpp::List operator_list)
{
  
  operator_list["nIter"] = 1;
  operator_list["estimate_theta"] = 1;
  MaternOperator2D* Kobj = new MaternOperator2D;
  Kobj->initFromList(operator_list, Rcpp::List::create(Rcpp::Named("use.chol") = 1));
  
  Process *process = NULL;
  process  = new GaussianProcess;
  process->initFromList(process_list, Kobj->h);
  Rcpp::List V_list  = Rcpp::as<Rcpp::List>  (process_list["V"]);
  Eigen::VectorXd V = Rcpp::as<Eigen::VectorXd>( V_list[0]);
  Kobj->gradient_init(0, 0);
  
  Kobj->gradient_add( process->Xs[0],
                      V.cwiseInverse(),
                      process->mean_X(0),
                      0,
                      1);
  
  Eigen::VectorXd g =  Kobj->get_gradient();
  Rcpp::List olist;

  olist["dtau1"] = g[0];
  olist["dtau2"] = g[1];
  olist["dkappa1"] = g[2];
  olist["dkappa2"] = g[3];
  olist["drho"] = g[4];
  olist["dtheta"] = g[5];
  
  olist["d2tau1"] = Kobj->ddtau1;
  olist["d2tau2"] = Kobj->ddtau2;
  olist["d2kappa1"] = Kobj->ddkappa1;
  olist["d2kappa2"] = Kobj->ddkappa2;
  olist["d2rho"] = Kobj->ddrho;
  olist["d2theta"] = Kobj->ddtheta;
  
  delete process;
  delete Kobj;
  return(olist);
}
