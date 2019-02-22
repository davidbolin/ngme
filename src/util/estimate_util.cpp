#include "estimate_util.h"
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
                              MixedEffect       & mixobj,
                              operatorMatrix    & Kobj,
                              MeasurementError  & errObj,
                              Process           & process,
                              int debug,
                              Eigen::VectorXd & z,
                              gig & rgig,
                              std::vector<  solver* >  & Solver)
{
  Eigen::VectorXd b;
  if(process_active){
    b.setZero(Kobj.d[i]);
  }

  Eigen::VectorXd  res = Y;
  //***************************************
  //***************************************
  //   building the residuals and sampling
  //***************************************
  //***************************************

  // removing fixed effect from Y
  // remove assymetric from errobj
  mixobj.remove_cov(i, res);
  errObj.remove_asym(i, res);
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

    K = Eigen::SparseMatrix<double,0,int>(Kobj.Q[i]);

    Q = K.transpose();
    Q =  Q * iV.asDiagonal();
    Q =  Q * K;

    res += A * process.Xs[i]; //res is now Y - Bfix*beta_fix - Brand*U = A*X + E
    if(debug)
      Rcpp::Rcout << "estimate::sample X\n";

    //Sample X|Y, V, sigma
    if(sampleX){
      if(errObj.noise == "Normal")
        process.sample_X(i, z, res, Q, K, A, errObj.sigma, *Solver[i]);
      else
        process.sample_Xv2( i, z, res, Q, K, A, errObj.sigma, *Solver[i], errObj.Vs[i].cwiseInverse());
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

  // add assym
  errObj.add_asym(i, res);
  if(errObj.noise != "Normal")
    errObj.sampleV(i, res);
  errObj.remove_asym(i, res);


  // rem assym
  
  return(res);
}

/*
internal function handling all gradient caculations

*/

void grad_caculations(int i,
                      Eigen::VectorXd&  res,
                      Eigen::SparseMatrix<double,0,int>& A,
                      double w,
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
  
  Eigen::VectorXd sigmas; 

  if(errObj.nsSigma)
    sigmas = errObj.sigmas[i];
  
  if(errObj.noise != "Normal"){
    mixobj.gradient2(i,
                     res,
                     errObj.Vs[i].cwiseInverse(),
                     sigmas,
                     2 * log(errObj.sigma),
                     errObj.EiV,
                     w, 
                     use_EU,
                     errObj.nsSigma);
    if(estimate_fisher)
      Fisher_information.block(0, 0, mixobj.npars + 1, mixobj.npars + 1) += mixobj.d2Given2(i,res,errObj.Vs[i].cwiseInverse(), 2 * log(errObj.sigma),errObj.EiV,w);
  }else{
    mixobj.gradient(i,
                    res,
                    2 * log(errObj.sigma),
                    w, 
                    use_EU);
    if(estimate_fisher){
      Fisher_information.block(0, 0, mixobj.npars + 1, mixobj.npars + 1) += mixobj.d2Given(i,res,2 * log(errObj.sigma),w);
    }
  }
  

  mixobj.remove_inter(i, res);

  if(debug)
    Rcpp::Rcout << "estimate::gradcalc::errObj \n";
  // measurent error  gradient
  //TODO:: ADDD SCALING WITH W FOR ERROR GRADIENT
  //errObj.add_asym(i, res);
  errObj.gradient(i, res, w);
  if( (errObj.noise != "Normal") && (estimate_fisher > 0) && (errObj.npars > 1) )
    Fisher_information.block(mixobj.npars + 1, mixobj.npars + 1, errObj.npars - 1, errObj.npars - 1) += errObj.d2Given(i, res, w);
  
  // add ass
  //errObj.remove_asym(i, res);
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

    K = Eigen::SparseMatrix<double,0,int>(Kobj.Q[i]);

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