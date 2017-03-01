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

  int debug = 0;
  //**********************************
  //      basic parameter
  //**********************************

  int nSim       = Rcpp::as< int > (in_list["nSim"]);
  int nBurnin    = Rcpp::as< int > (in_list["nBurnin"] );
  int silent     = Rcpp::as< int > (in_list["silent"]);
  int n_threads  = Rcpp::as< int > (in_list["n_threads"]);
  int mix_samp  = Rcpp::as< int > (in_list["mix_samp"]);
  int use_random_effect = Rcpp::as< int > (in_list["use_random_effect"]);
  int use_process = 0;
  if(in_list.containsElementNamed("processes_list"))
    use_process = 1;

  int predict_derivative = Rcpp::as< int > (in_list["predict_derivative"]);
  double deriv_scale = 1.0;
  if(predict_derivative == 1){
    deriv_scale = Rcpp::as< double > (in_list["derivative_scaling"]);
  }


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
  std::vector< Eigen::SparseMatrix<double,0,int> > As_pred_1( nindv);

  std::vector< Eigen::VectorXd > Ys( nindv);
  std::vector< Eigen::MatrixXd > pred_ind(nindv);
  std::vector< Eigen::MatrixXd > obs_ind(nindv);
  std::vector< Eigen::MatrixXd > Bfixed_pred(nindv);
  std::vector< Eigen::MatrixXd > Brandom_pred(nindv);

  std::vector< Eigen::MatrixXd > Bfixed_pred_1(nindv);
  std::vector< Eigen::MatrixXd > Brandom_pred_1(nindv);


  int count;
  count = 0;
  for( List::iterator it = obs_list.begin(); it != obs_list.end(); ++it ) {
    List obs_tmp = Rcpp::as<Rcpp::List>(*it);
    if(use_process == 1){
      As[count]            = Rcpp::as<Eigen::SparseMatrix<double,0,int> >(obs_tmp["A"]);
      As_pred[count]       = Rcpp::as<Eigen::SparseMatrix<double,0,int> >(obs_tmp["Apred"]);
    }
    if(predict_derivative == 1){
      if(use_process == 1){
        As_pred_1[count]       = Rcpp::as<Eigen::SparseMatrix<double,0,int> >(obs_tmp["Apred1"]);
      }
      Bfixed_pred_1[count]   = Rcpp::as<Eigen::MatrixXd>(obs_tmp["Bfixed_pred1"]);
      if(use_random_effect == 1){
        Brandom_pred_1[count]  = Rcpp::as<Eigen::MatrixXd>(obs_tmp["Brandom_pred1"]);
      }
    }
    pred_ind[count]      = Rcpp::as<Eigen::MatrixXd>(obs_tmp["pred_ind"]);
    obs_ind[count]       = Rcpp::as<Eigen::MatrixXd>(obs_tmp["obs_ind"]);
    Ys[count]            = Rcpp::as<Eigen::VectorXd>(obs_tmp["Y"]);
    if(use_random_effect == 1){
      Brandom_pred[count]  = Rcpp::as<Eigen::MatrixXd>(obs_tmp["Brandom_pred"]);
    }
    Bfixed_pred[count]   = Rcpp::as<Eigen::MatrixXd>(obs_tmp["Bfixed_pred"]);
    count++;
  }


  //**********************************
  //operator setup
  //***********************************
  operatorMatrix* Kobj;
  int common_grid = 1;
  std::vector<  cholesky_solver >  Solver( nindv);
  if(use_process == 1){
    if(silent == 0){
      Rcpp::Rcout << " Setup operator\n";
    }
    Rcpp::List operator_list  = Rcpp::as<Rcpp::List> (in_list["operator_list"]);
    std::string type_operator = Rcpp::as<std::string>(operator_list["type"]);
    operator_list["nIter"] = 1;

    operator_select(type_operator, &Kobj);
    Kobj->initFromList(operator_list, List::create(Rcpp::Named("use.chol") = 1));

    if(Kobj->nop>1){
      common_grid = 0;
    }


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
  }
  //**********************************
  // mixed effect setup
  //***********************************
  if(silent == 0){
    Rcpp::Rcout << " Setup mixed effect: RE = " << use_random_effect << " \n";
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
  Process *process;
  std::string type_processes;
  if(use_process == 1){

    if(silent == 0){
      Rcpp::Rcout << " Setup process\n";
    }
    Rcpp::List processes_list   = Rcpp::as<Rcpp::List>  (in_list["processes_list"]);
    type_processes  = Rcpp::as<std::string> (processes_list["noise"]);


    if (type_processes != "Normal"){
      process  = new GHProcess;
    }else{ process  = new GaussianProcess;}

    process->initFromList(processes_list, Kobj->h);
  }

  //**********************************
  // OpenMP setup
  //***********************************
  if(silent == 0){
    Rcpp::Rcout << " Setup OMP: " << n_threads << "\n";
  }
  #ifdef _OPENMP
    Eigen::initParallel();
    //Eigen::setNbThreads(0);
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

  Eigen::VectorXd b, Ysim;
  if(use_process == 1){
    Eigen::VectorXd  z;
    z.setZero(Kobj->d[0]);
    b.setZero(Kobj->d[0]);
  }


  std::vector<int> longInd;
  for (int i=0; i< nindv; i++){
    longInd.push_back(i);
  }



  std::vector< Eigen::MatrixXd > WVec(nindv);
  std::vector< Eigen::MatrixXd > VVec(nindv);
  std::vector< Eigen::MatrixXd > XVec(nindv);
  std::vector< Eigen::MatrixXd > XVec_deriv(nindv);
  std::vector< Eigen::MatrixXd > WVec_deriv(nindv);
  std::vector< Eigen::MatrixXd > YVec(nindv);

  /**
   * Timing setup
   */
  const double ticks_per_ms = static_cast<double>(CLOCKS_PER_SEC);


  int rank = 0;
  double percent_done = 0;

  #pragma omp parallel for private(K,Q,rank,debug)
  for(int i = 0; i < nindv; i++ ) {
    #ifdef _OPENMP
      rank = omp_get_thread_num();
    #endif
    clock_t start = clock();
    percent_done++;
    if(silent == 0){
      std::stringstream stream;
      if(nP>1){
        stream << " Prediction " << 100*percent_done/nindv << " % done (" << rank<< ")\n";
      } else {
        stream << " Prediction " << 100*percent_done/nindv << " % done\n";
      }
      Rcpp::Rcout << stream.str();
    }
    XVec[i].setZero(Bfixed_pred[i].rows(), nSim);
    WVec[i].setZero(Bfixed_pred[i].rows(), nSim);
    if(predict_derivative == 1){
      XVec_deriv[i].setZero(Bfixed_pred[i].rows(), nSim);
      WVec_deriv[i].setZero(Bfixed_pred[i].rows(), nSim);
    }
    YVec[i].setZero(Bfixed_pred[i].rows(), nSim);
    VVec[i].setZero(Bfixed_pred[i].rows(), nSim);
    Eigen::MatrixXd random_effect;
    if(use_random_effect == 1){
      random_effect = mixobj->Br[i];
    }
    Eigen::MatrixXd fixed_effect = mixobj->Bf[i];

    for(int ipred = 0; ipred < pred_ind[i].rows(); ipred++){
      int n_obs = obs_ind[i](ipred,1) - obs_ind[i](ipred,0);
      Eigen::VectorXd Y = Ys[i];
      Eigen::SparseMatrix<double,0,int> A = As[i];
      if(n_obs > 0){
        if(use_process == 1){
          A = As[i].middleRows(obs_ind[i](ipred,0),obs_ind[i](ipred,1));
        }
        Y = Ys[i].segment(obs_ind[i](ipred,0),obs_ind[i](ipred,1));
        if(use_random_effect == 1){
          mixobj->Br[i] = random_effect.middleRows(obs_ind[i](ipred,0),obs_ind[i](ipred,1));
        }
        mixobj->Bf[i] = fixed_effect.middleRows(obs_ind[i](ipred,0),obs_ind[i](ipred,1));
      }
      for(int ii = 0; ii < nSim + nBurnin; ii ++){
        Eigen::VectorXd res = Y;
        if(n_obs>0){
          mixobj->remove_cov(i, res);
          if(use_process == 1){
            if(n_obs>0)
              res -= A * process->Xs[i];
          }
          //***********************************
          // mixobj sampling
          //***********************************
          if(debug){
            Rcpp::Rcout << "Sample mixed effect (" << rank  << ")\n";
          }

          for(int kkk=0;kkk<mix_samp;kkk++){
            if(type_MeasurementError == "Normal"){
              mixobj->sampleU_par( i, res, 2 * log(errObj->sigma),random_engine[rank]);
            } else {
              //mixobj->sampleU2( i, res, errObj->Vs[i].cwiseInverse(), 2 * log(errObj->sigma));
              mixobj->sampleU2_par( i, res, errObj->Vs[i].segment(obs_ind[i](ipred,0),obs_ind[i](ipred,1)).cwiseInverse(),
                                    random_engine[rank], 2 * log(errObj->sigma));
            }
          }


          mixobj->remove_inter(i, res);

        } else {
          for(int kkk=0;kkk<mix_samp;kkk++){
            mixobj->simulate(Y, i);
          }
        }

        //***********************************
        // sampling processes
        //***********************************
        Eigen::SparseMatrix<double, 0, int> Ki;
        Eigen::SparseMatrix<double,0,int> Ai;
        if(use_process == 1){
          if(debug){
            Rcpp::Rcout << "Compute operator (" << rank  << ")\n";
          }

          Eigen::SparseMatrix<double, 0, int> Qi;
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

          if(debug){
            Rcpp::Rcout << "Sample normals (" << rank  << ")\n";
          }
          int d = 1;
          if(common_grid == 1){
            d = Kobj->d[0];
          } else {
            d = Kobj->d[i];
          }

          Eigen::VectorXd zi;
          zi.setZero(d);
          for(int j =0; j < d; j++)
            zi[j] =  normal(random_engine[rank]);
          if(n_obs>0)
            res += A * process->Xs[i];

          if(debug){
            Rcpp::Rcout << "Sample process (" << rank  << ")\n";
          }
          if(n_obs>0){
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
          } else {
            process->simulate(i,zi,A,Ki,Y,Solver[i]);
          }
        }

        if(debug){
          Rcpp::Rcout << "Sample variances (" << rank  << ")\n";
        }
        if(use_process == 1){
          if(n_obs>0){
              process->sample_V(i, rgig[rank], Ki);
            } else {
              process->simulate_V(i, rgig[rank]);
          }
        }
        //***********************************
        // random variance noise sampling
        //***********************************

        if(type_MeasurementError != "Normal"){
          if(n_obs>0){
            errObj->sampleV(i, res,obs_ind[i](ipred,0)+obs_ind[i](ipred,1));
          }
        }

          // save samples
          if(ii >= nBurnin){
            if(debug){
              Rcpp::Rcout << "Save samples\n";
            }
            int n_pred = pred_ind[i](ipred,1) - pred_ind[i](ipred,0);
          if(use_process == 1){
            Ai = As_pred[i];
            Ai = Ai.middleRows(pred_ind[i](ipred,0),pred_ind[i](ipred,1));
          }

          Eigen::VectorXd random_effect;
          if(use_random_effect == 1){
            random_effect= Bfixed_pred[i]*mixobj->beta_fixed + Brandom_pred[i]*(mixobj->U.col(i)+mixobj->beta_random);
          } else {
            random_effect = Bfixed_pred[i]*mixobj->beta_fixed;
          }

          Eigen::VectorXd random_effect_c = random_effect.segment(pred_ind[i](ipred,0),pred_ind[i](ipred,1));
          Eigen::VectorXd mNoise = errObj->simulate_par(random_effect_c,random_engine[rank]);
          Eigen::VectorXd AX;

          if(use_process == 1){
            AX = Ai * process->Xs[i];
            WVec[i].block(pred_ind[i](ipred,0), ii - nBurnin, pred_ind[i](ipred,1), 1) = AX;
            XVec[i].block(pred_ind[i](ipred,0), ii - nBurnin, pred_ind[i](ipred,1), 1) = random_effect_c + AX;
            YVec[i].block(pred_ind[i](ipred,0), ii - nBurnin, pred_ind[i](ipred,1), 1) = random_effect_c + AX + mNoise;
            VVec[i].block(pred_ind[i](ipred,0), ii - nBurnin, pred_ind[i](ipred,1), 1) = Ai * process->Vs[i];
          }  else {
            XVec[i].block(pred_ind[i](ipred,0), ii - nBurnin, pred_ind[i](ipred,1), 1) = random_effect_c;
            YVec[i].block(pred_ind[i](ipred,0), ii - nBurnin, pred_ind[i](ipred,1), 1) = random_effect_c + mNoise;
          }

          if(predict_derivative == 1){
            Eigen::VectorXd random_effect_1;
            if(use_random_effect == 1){
              random_effect_1= Bfixed_pred_1[i]*mixobj->beta_fixed + Brandom_pred_1[i]*(mixobj->U.col(i)+mixobj->beta_random);
            } else {
              random_effect_1 = Bfixed_pred_1[i]*mixobj->beta_fixed;
            }

            Eigen::VectorXd random_effect_c_1 = random_effect_1.segment(pred_ind[i](ipred,0),pred_ind[i](ipred,1));
            if(use_process == 1){
              Eigen::SparseMatrix<double,0,int> Ai_1 = As_pred_1[i];
              Ai_1 = Ai_1.middleRows(pred_ind[i](ipred,0),pred_ind[i](ipred,1));
              Eigen::VectorXd AX_1 = Ai_1 * process->Xs[i];
              WVec_deriv[i].block(pred_ind[i](ipred,0), ii - nBurnin, pred_ind[i](ipred,1), 1) = AX_1 - AX;
              XVec_deriv[i].block(pred_ind[i](ipred,0), ii - nBurnin, pred_ind[i](ipred,1), 1) = random_effect_c_1 + AX_1 - random_effect_c - AX;
            } else {
              XVec_deriv[i].block(pred_ind[i](ipred,0), ii - nBurnin, pred_ind[i](ipred,1), 1) = random_effect_c_1 - random_effect_c;
            }
          }
        }
      }

    }
    if(predict_derivative == 1){
      WVec_deriv[i] /= deriv_scale;
      XVec_deriv[i] /= deriv_scale;
    }

    double time_Ma = static_cast<double>(clock()-start)  / ticks_per_ms;
    if(silent == 0){
      std::stringstream stream;
      stream << ", time = " << time_Ma;
      Rcpp::Rcout << stream.str();
      Rcpp::Rcout << "\n";
    }
  }

  //Rcpp::Rcout << "store results\n";
  // storing the results
  Rcpp::List out_list;
  out_list["YVec"] = YVec;
  out_list["XVec"] = XVec;
  out_list["WVec"] = WVec;
  out_list["VVec"] = VVec;
  if(predict_derivative == 1)
  {
    out_list["XVec_deriv"] = XVec_deriv;
    out_list["WVec_deriv"] = WVec_deriv;
  }
  return(out_list);
}
