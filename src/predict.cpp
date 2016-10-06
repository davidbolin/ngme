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
#ifdef _OPENMP
#include<omp.h>
#endif
using namespace Rcpp;
#define max(a,b) (((a)>(b))?(a):(b))
#define min(a,b) (((a)<(b))?(a):(b))



// [[Rcpp::export]]
List predictLong_cpp(Rcpp::List in_list)
{


  //**********************************
  //      basic parameter
  //**********************************

  int nSim       = Rcpp::as< int > (in_list["nSim"]);
  int nBurnin    = Rcpp::as< int > (in_list["nBurnin"] );
  int silent     = Rcpp::as< int > (in_list["silent"]);
  int n_threads  = Rcpp::as< int > (in_list["n_threads"]);

  //**********************************
  //     setting up the main data
  //**********************************
  if(silent == 0){
    Rcpp::Rcout << " Setup data\n";
  }
  Rcpp::List obs_list  = Rcpp::as<Rcpp::List> (in_list["obs_list"]);
  int nindv = obs_list.length(); //count number of patients
  std::vector< Eigen::SparseMatrix<double,0,int> > As( nindv);
  std::vector< Eigen::SparseMatrix<double,0,int> > As_pred( nindv);
  std::vector< Eigen::VectorXd > Ys( nindv);
  std::vector< Eigen::MatrixXd > pred_ind(nindv);
  std::vector< Eigen::MatrixXd > obs_ind(nindv);
  std::vector< Eigen::MatrixXd > Bfixed_pred(nindv);
  std::vector< Eigen::MatrixXd > Brandom_pred(nindv);
  int count;
  count = 0;
  for( List::iterator it = obs_list.begin(); it != obs_list.end(); ++it ) {
    List obs_tmp = Rcpp::as<Rcpp::List>(*it);
    As[count]            = Rcpp::as<Eigen::SparseMatrix<double,0,int> >(obs_tmp["A"]);
    As_pred[count]       = Rcpp::as<Eigen::SparseMatrix<double,0,int> >(obs_tmp["Apred"]);
    pred_ind[count]      = Rcpp::as<Eigen::MatrixXd>(obs_tmp["pred_ind"]);
    obs_ind[count]       = Rcpp::as<Eigen::MatrixXd>(obs_tmp["obs_ind"]);
    Ys[count]            = Rcpp::as<Eigen::VectorXd>(obs_tmp["Y"]);
    Brandom_pred[count]  = Rcpp::as<Eigen::MatrixXd>(obs_tmp["Brandom_pred"]);
    Bfixed_pred[count]   = Rcpp::as<Eigen::MatrixXd>(obs_tmp["Bfixed_pred"]);
    count++;
  }


  //**********************************
  //operator setup
  //***********************************
  if(silent == 0){
    Rcpp::Rcout << " Setup operator\n";
  }
  Rcpp::List operator_list  = Rcpp::as<Rcpp::List> (in_list["operator_list"]);
  std::string type_operator = Rcpp::as<std::string>(operator_list["type"]);
  operator_list["nIter"] = 1;
  operatorMatrix* Kobj;
  //Kobj = new MaternMatrixOperator;
  operator_select(type_operator, &Kobj);
  Kobj->initFromList(operator_list, List::create(Rcpp::Named("use.chol") = 1));


  int common_grid = 1;
  if(Kobj->nop>1){
    common_grid = 0;
  }

  std::vector<  cholesky_solver >  Solver( nindv);
  Eigen::SparseMatrix<double, 0, int> Q, K;
//  Rcpp::Rcout << "d = " << Kobj->d[0] << "\n";
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
  //**********************************
  // mixed effect setup
  //***********************************
  if(silent == 0){
    Rcpp::Rcout << " Setup mixed effect\n";
  }
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




  //**********************************
  // stochastic processes setup
  //***********************************
  if(silent == 0){
    Rcpp::Rcout << " Setup process\n";
  }
  Rcpp::List processes_list   = Rcpp::as<Rcpp::List>  (in_list["processes_list"]);
  std::string type_processes  = Rcpp::as<std::string> (processes_list["noise"]);


  Process *process;

  if (type_processes != "Normal"){
    process  = new GHProcess;
  }else{ process  = new GaussianProcess;}

  process->initFromList(processes_list, Kobj->h);


  //**********************************
  // OpenMP setup
  //***********************************
  if(silent == 0){
    Rcpp::Rcout << " Setup OMP: " << n_threads << "\n";
  }
#ifdef _OPENMP
  Eigen::initParallel();
  const int max_nP = omp_get_num_procs();
  int nPtmp;
  if(n_threads == 0){
    nPtmp = max_nP;
  } else {
    nPtmp = min(max_nP, max(n_threads,1));
  }
  const int nP = nPtmp;
  omp_set_num_threads(nP);
#else
  const int nP = 1;
#endif


  /*
  Simulation objects
  */

  std::normal_distribution<double> normal;

  std::vector<std::mt19937> random_engine(nP);
  std::vector<gig> rgig(nP);
  for (int i = 0; i < nP; ++i) {
    random_engine[i].seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    rgig[i].seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());
  }


  Eigen::VectorXd  z;
  z.setZero(Kobj->d[0]);

  Eigen::VectorXd b, Ysim;
  b.setZero(Kobj->d[0]);

  std::vector<int> longInd;
  for (int i=0; i< nindv; i++) longInd.push_back(i);



  std::vector< Eigen::MatrixXd > WVec(nindv);
  std::vector< Eigen::MatrixXd > VVec(nindv);
  std::vector< Eigen::MatrixXd > XVec(nindv);
  std::vector< Eigen::MatrixXd > YVec(nindv);

  int rank = 1;
  double percent_done = 0;
  int  d;
  #pragma omp parallel for private(K,Q)
  for(int i = 0; i < nindv; i++ ) {
    #ifdef _OPENMP
      rank = omp_get_thread_num();
    #endif

    percent_done++;
    if(silent == 0){
      std::stringstream stream;
      stream << " Prediction " << 100*percent_done/nindv << " % done (" << rank<< ")\n";
      Rcpp::Rcout << stream.str();
    }

    XVec[i].resize(As_pred[i].rows(), nSim);
    WVec[i].resize(As_pred[i].rows(), nSim);
    YVec[i].resize(As_pred[i].rows(), nSim);
    VVec[i].resize(As_pred[i].rows(), nSim);
    Eigen::MatrixXd random_effect = mixobj->Br[i];
    Eigen::MatrixXd fixed_effect = mixobj->Bf[i];
    for(int ipred = 0; ipred < pred_ind[i].rows(); ipred++){
      //extract data to use for prediction:
      Eigen::SparseMatrix<double,0,int> A = As[i].middleRows(obs_ind[i](ipred,0),obs_ind[i](ipred,1));
      Eigen::VectorXd Y = Ys[i].segment(obs_ind[i](ipred,0),obs_ind[i](ipred,1));
      mixobj->Br[i] = random_effect.middleRows(obs_ind[i](ipred,0),obs_ind[i](ipred,1));
      mixobj->Bf[i] = fixed_effect.middleRows(obs_ind[i](ipred,0),obs_ind[i](ipred,1));

      for(int ii = 0; ii < nSim + nBurnin; ii ++){
        //Rcpp::Rcout << "iter = " << ii << "\n";
        Eigen::VectorXd res = Y;
        //   building the residuals and sampling

        // removing fixed effect from Y
        mixobj->remove_cov(i, res);
        res -= A * process->Xs[i];


        //***********************************
        // mixobj sampling
        //***********************************
        //Rcpp::Rcout << "here2\n";
        if(type_MeasurementError == "Normal"){
          //mixobj->sampleU( i, res, 2 * log(errObj->sigma),random_engine[rank]);
          mixobj->sampleU_par( i, res, 2 * log(errObj->sigma),random_engine[rank]);

        } else {
          //errObj->Vs[i] = errObj->Vs[i].head(ipred+1);
          // mixobj->sampleU2( i,
          //                  res,
          //                  errObj->Vs[i].segment(obs_ind[i](ipred,0),obs_ind[i](ipred,1)).cwiseInverse(),
          //                  2 * log(errObj->sigma));
          mixobj->sampleU2_par( i,
                            res,
                            errObj->Vs[i].segment(obs_ind[i](ipred,0),obs_ind[i](ipred,1)).cwiseInverse(),
                            random_engine[rank],
                            2 * log(errObj->sigma));
        }


        mixobj->remove_inter(i, res);

        //***********************************
        // sampling processes
        //***********************************
        //Rcpp::Rcout << "here3\n";
        Eigen::SparseMatrix<double, 0, int> Ki, Qi;
        if(common_grid){
          Ki = Eigen::SparseMatrix<double,0,int>(Kobj->Q[0]);
        } else {
          Ki = Eigen::SparseMatrix<double,0,int>(Kobj->Q[i]);
        }
        Eigen::VectorXd iV(process->Vs[i].size());
        iV.array() = process->Vs[i].array().inverse();

        Qi = Ki.transpose();
        Qi =  Qi * iV.asDiagonal();
        Qi =  Qi * Ki;
        //before here

        // Rcpp::Rcout << "here4\n
        int d = 1;
        if(common_grid == 1){
          d = Kobj->d[0];
        } else {
          d = Kobj->d[i];
        }
        // Rcpp::Rcout << "here5\n";
        Eigen::VectorXd zi;
        zi.setZero(d);
        for(int j =0; j < d; j++)
          zi[j] =  normal(random_engine[rank]);
        //Rcpp::Rcout << "here5.1\n";
        res += A * process->Xs[i];
        //Sample X|Y, V, sigma

        //Rcpp::Rcout << "here6\n";
        //Rcpp::Rcout << d << "\n";
        //Rcpp::Rcout << zi.size() << "\n";
        if(type_MeasurementError == "Normal"){
          process->sample_X(i,
                            zi,
                            res,
                            Qi,
                            Ki,
                            A,
                            errObj->sigma,
                            Solver[i]);
        } else {
          process->sample_Xv2( i,
                               zi,
                               res,
                               Qi,
                               Ki,
                               A,
                               errObj->sigma,
                               Solver[i],
                               errObj->Vs[i].segment(obs_ind[i](ipred,0),obs_ind[i](ipred,1)).cwiseInverse());
        }
        res -= A * process->Xs[i];

      //after here



        //if(res.cwiseAbs().sum() > 1e16){
        //  throw("res outof bound\n");
        //}
        // sample V| X
        process->sample_V(i, rgig[rank], Ki);
        //***********************************
        // random variance noise sampling
        //***********************************
        //Rcpp::Rcout << "here8\n";
        if(type_MeasurementError != "Normal"){
          errObj->sampleV(i, res,obs_ind[i](ipred,0)+obs_ind[i](ipred,1));
        }

        // save samples
        if(ii >= nBurnin){
          Eigen::SparseMatrix<double,0,int> Ai = As_pred[i];
          Ai = Ai.middleRows(pred_ind[i](ipred,0),pred_ind[i](ipred,1));
          Eigen::VectorXd random_effect = Bfixed_pred[i]*mixobj->beta_fixed + Brandom_pred[i]*(mixobj->U.col(i)+mixobj->beta_random);
          Eigen::VectorXd random_effect_c = random_effect.segment(pred_ind[i](ipred,0),pred_ind[i](ipred,1));
          Eigen::VectorXd AX = Ai * process->Xs[i];
          Eigen::VectorXd mNoise = errObj->simulate_par(AX,random_engine[rank]);

          WVec[i].block(pred_ind[i](ipred,0), ii - nBurnin, pred_ind[i](ipred,1), 1) = AX;
          XVec[i].block(pred_ind[i](ipred,0), ii - nBurnin, pred_ind[i](ipred,1), 1) = random_effect_c + AX;
          YVec[i].block(pred_ind[i](ipred,0), ii - nBurnin, pred_ind[i](ipred,1), 1) = random_effect_c + AX + mNoise;
          //V process
          VVec[i].block(pred_ind[i](ipred,0), ii - nBurnin, pred_ind[i](ipred,1), 1) = Ai * process->Vs[i];
        }
      }
    }
  }

  //Rcpp::Rcout << "store results\n";
  // storing the results
  Rcpp::List out_list;
  out_list["YVec"] = YVec;
  out_list["XVec"] = XVec;
  out_list["WVec"] = WVec;
  out_list["VVec"] = VVec;
  return(out_list);
}
