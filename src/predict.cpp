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


// [[Rcpp::export]]
List predictLong_cpp(Rcpp::List in_list)
{


  //**********************************
  //      basic parameter
  //**********************************

  int nSim       = Rcpp::as< int > (in_list["nSim"]);
  int nBurnin    = Rcpp::as< int > (in_list["nBurnin"] );
  int silent     = Rcpp::as< int > (in_list["silent"]);

  //**********************************
  //     setting up the main data
  //**********************************
  //if(silent == 0){
  //  Rcpp::Rcout << " Setup data\n";
  //}
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
  //if(silent == 0){
  //  Rcpp::Rcout << " Setup operator\n";
  //}
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
  //if(silent == 0){
  //  Rcpp::Rcout << " Setup mixed effect\n";
  //}
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
  //if(silent == 0){
  //  Rcpp::Rcout << " Setup noise\n";
  //}
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
  //if(silent == 0){
  //  Rcpp::Rcout << " Setup process\n";
  //}
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
  random_engine.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());
  gig rgig;
  rgig.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());
  Eigen::VectorXd  z;
  z.setZero(Kobj->d[0]);

  Eigen::VectorXd b, Ysim;
  b.setZero(Kobj->d[0]);

  std::vector<int> longInd;
  for (int i=0; i< nindv; i++) longInd.push_back(i);



  std::vector< Eigen::MatrixXd > WVec(nindv);
  std::vector< Eigen::MatrixXd > XVec(nindv);

  for(int i = 0; i < nindv; i++ ) {
    if(silent == 0){
      Rcpp::Rcout << " Compute prediction for patient " << i << "\n";
    }

    XVec[i].resize(As_pred[i].rows(), nSim);
    WVec[i].resize(As_pred[i].rows(), nSim);
    Eigen::MatrixXd random_effect = mixobj->Br[i];
    Eigen::MatrixXd fixed_effect = mixobj->Bf[i];
    for(int ipred = 0; ipred < pred_ind[i].rows(); ipred++){
      //if(silent == 0){
      //  Rcpp::Rcout << " location = " << ipred << ": "<< obs_ind[i](ipred,0)<< " " << obs_ind[i](ipred,1) << "\n";
      //  Rcpp::Rcout << " A size = " <<  As[i].rows() << " " << As[i].cols() << "\n";
      //  Rcpp::Rcout << " Y size = " <<  Ys[i].size() << ", re = " << random_effect.rows() << ", fe = "<< fixed_effect.size() << "\n";
      //}
      //extract data to use for prediction:
      Eigen::SparseMatrix<double,0,int> A = As[i].middleRows(obs_ind[i](ipred,0),obs_ind[i](ipred,1));
      Eigen::VectorXd  Y = Ys[i].segment(obs_ind[i](ipred,0),obs_ind[i](ipred,1));
      mixobj->Br[i] = random_effect.middleRows(obs_ind[i](ipred,0),obs_ind[i](ipred,1));
      mixobj->Bf[i] = fixed_effect.middleRows(obs_ind[i](ipred,0),obs_ind[i](ipred,1));
      //Rcpp::Rcout  << "here1111\n";
      for(int ii = 0; ii < nSim + nBurnin; ii ++){
        //Rcpp::Rcout << "iter = " << ii << "\n";
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
        //Rcpp::Rcout << "here2\n";
        if(type_MeasurementError == "Normal")
          mixobj->sampleU( i, res, 2 * log(errObj->sigma));
        else
          //errObj->Vs[i] = errObj->Vs[i].head(ipred+1);
          mixobj->sampleU2( i,
                            res,
                            errObj->Vs[i].segment(obs_ind[i](ipred,0),obs_ind[i](ipred,1)).cwiseInverse(),
                            2 * log(errObj->sigma));

          mixobj->remove_inter(i, res);

        //***********************************
        // sampling processes
        //***********************************
        //Rcpp::Rcout << "here3\n";
        if(common_grid){
          K = Eigen::SparseMatrix<double,0,int>(Kobj->Q[0]);
        } else {
          K = Eigen::SparseMatrix<double,0,int>(Kobj->Q[i]);
        }
        Eigen::VectorXd iV(process->Vs[i].size());
        iV.array() = process->Vs[i].array().inverse();

        Q = K.transpose();
        Q =  Q * iV.asDiagonal();
        Q =  Q * K;

        int d;
        if(common_grid == 1){
          d = Kobj->d[0];
        } else {
          d = Kobj->d[i];
        }
        for(int j =0; j < d; j++)
          z[j] =  normal(random_engine);

        res += A * process->Xs[i];
        //Sample X|Y, V, sigma

        if(type_MeasurementError == "Normal")
          process->sample_X(i,
                            z,
                            res,
                            Q,
                            K,
                            A,
                            errObj->sigma,
                            Solver[i]);
        else
          process->sample_Xv2( i,
                               z,
                               res,
                               Q,
                               K,
                               A,
                               errObj->sigma,
                               Solver[i],
                               errObj->Vs[i].segment(obs_ind[i](ipred,0),obs_ind[i](ipred,1)).cwiseInverse());
        res -= A * process->Xs[i];

        //if(res.cwiseAbs().sum() > 1e16){
        //  throw("res outof bound\n");
        //}

        // sample V| X
        process->sample_V(i, rgig, K);

        //***********************************
        // random variance noise sampling
        //***********************************
        //Rcpp::Rcout << "here4\n";
        if(type_MeasurementError != "Normal"){
          errObj->sampleV(i, res,obs_ind[i](ipred,0)+obs_ind[i](ipred,1));
        }
        // save samples
        if(ii >= nBurnin){
          //if(silent == 0){
          //  Rcpp::Rcout << " save samples << " << i << ", " << ipred << ", " << ii << "\n";
          //}
          //Rcpp::Rcout << "here 1\n";
          Eigen::SparseMatrix<double,0,int> Ai = As_pred[i];
          Ai = Ai.middleRows(pred_ind[i](ipred,0),pred_ind[i](ipred,1));
          //Rcpp::Rcout << "here 2\n";
          Eigen::VectorXd random_effect = Bfixed_pred[i]*mixobj->beta_fixed;
          //Rcpp::Rcout << "here 3\n";
          random_effect += Brandom_pred[i]*(mixobj->U.col(i)+mixobj->beta_random);
          //Rcpp::Rcout << "here 4\n";
          Eigen::VectorXd random_effect_c = random_effect.segment(pred_ind[i](ipred,0),pred_ind[i](ipred,1));
          Eigen::VectorXd AX = Ai * process->Xs[i];
          //Rcpp::Rcout << "here 5\n";
          WVec[i].block(pred_ind[i](ipred,0), ii - nBurnin, pred_ind[i](ipred,1), 1) = AX;
          XVec[i].block(pred_ind[i](ipred,0), ii - nBurnin, pred_ind[i](ipred,1), 1) = random_effect_c + AX;
          //Rcpp::Rcout << "here 6\n";
        }
      }
    }
  }
  Rcpp::Rcout << "store results\n";
  // storing the results
  Rcpp::List out_list;
  out_list["XVec"] = XVec;
  out_list["WVec"] = WVec;
  return(out_list);
}
