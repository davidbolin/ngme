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
#include <Rcpp/Benchmark/Timer.h>
using namespace Rcpp;
#define max(a,b) (((a)>(b))?(a):(b))
#define min(a,b) (((a)<(b))?(a):(b))



// [[Rcpp::export]]
List predictLong_cpp(Rcpp::List in_list)
{

  int debug = 0;
  int timeit = 1;

  double time_setup = 0;
  double time_sample = 0;
  double time_post = 0;
  double time_process_sample = 0;
  double time_process_simulate = 0;
  double time_v_sample = 0;
  double time_operator = 0;
  //**********************************
  //      basic parameter
  //**********************************
  int nSim       = Rcpp::as< int > (in_list["nSim"]);
  int nBurnin    = Rcpp::as< int > (in_list["nBurnin"] );
  int silent     = Rcpp::as< int > (in_list["silent"]);
  int mix_samp  = Rcpp::as< int > (in_list["mix_samp"]);
  int ind_general  = Rcpp::as< int > (in_list["ind_general"]);
  int use_random_effect = Rcpp::as< int > (in_list["use_random_effect"]);
  int use_process = 0;
  if(in_list.containsElementNamed("processes_list"))
    use_process = 1;

  int predict_derivative = Rcpp::as< int > (in_list["predict_derivative"]);
  double deriv_scale = 1.0;
  if(predict_derivative == 1){
    deriv_scale = Rcpp::as< double > (in_list["derivative_scaling"]);
  }


  unsigned long seed = 0;
  if(in_list.containsElementNamed("seed"))
    seed = Rcpp::as< unsigned long > (in_list["seed"]);


  Rcpp::List obs_list  = Rcpp::as<Rcpp::List> (in_list["obs_list"]);

  //**********************************
  // mixed effect setup
  //***********************************
  if(silent == 0){
    Rcpp::Rcout << " Setup mixed effect: RE = " << use_random_effect << " \n";
  }
  Rcpp::List mixedEffect_list  = Rcpp::as<Rcpp::List> (in_list["mixedEffect_list"]);
  std::string type_mixedEffect = Rcpp::as<std::string> (mixedEffect_list["noise"]);
  MixedEffect *mixobj;

  if(type_mixedEffect == "Normal"){
    mixobj = new NormalMixedEffect;
  }else if(type_mixedEffect == "tdist") {
    mixobj   = new tdMixedEffect;
  }else if(type_mixedEffect == "NIG") {
    mixobj   = new NIGMixedEffect;
  } else {
    Rcpp::Rcout << "Wrong mixed effect distribution\n";
    return(0);
  }


  mixobj->initFromList(mixedEffect_list);


  //**********************************
  //     setting up the prediction data
  //**********************************
  if(silent == 0){
    Rcpp::Rcout << " Setup predict data\n";
  }


  int nindv = obs_list.length(); //count number of patients
  std::vector< Eigen::SparseMatrix<double,0,int> > As( nindv);
  std::vector< Eigen::SparseMatrix<double,0,int> > As_pred( nindv);
  std::vector< Eigen::SparseMatrix<double,0,int> > As_pred_1( nindv);

  std::vector< Eigen::VectorXd > Ys( nindv);
  std::vector< Eigen::MatrixXi > pred_ind(nindv);
  std::vector< Eigen::MatrixXi > obs_ind(nindv);
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
    if(predict_derivative == 1)
    {
      if(use_process == 1)
      {
        As_pred_1[count]       = Rcpp::as<Eigen::SparseMatrix<double,0,int> >(obs_tmp["Apred1"]);
      }
      Bfixed_pred_1[count]   = Rcpp::as<Eigen::MatrixXd>(obs_tmp["Bfixed_pred1"]);
      if(use_random_effect == 1)
      {
        Brandom_pred_1[count]  = Rcpp::as<Eigen::MatrixXd>(obs_tmp["Brandom_pred1"]);

        if(Brandom_pred_1[count].cols() != mixobj->beta_random.size())
        {
          Rcpp::Rcout << "predict: error dimension missmatch :\n Brandom_1[" << count << "].cols() = " << Brandom_pred_1[count].cols() << " != "
                      << "mixobj->beta_random.size() = "<<  mixobj->beta_random.size()  << "\n";
            throw("input error\n");
        }

      }
    if(Bfixed_pred_1[count].cols() != mixobj->beta_fixed.size())
    {
      Rcpp::Rcout << "predict: error dimension missmatch :\n Bfixed_1[" << count << "].cols() = " << Bfixed_pred_1[count].cols() << " != "
                  << "mixobj->beta_fixed.size() = "<<  mixobj->beta_fixed.size()  << "\n";
      throw("input error\n");

    }

    }

    pred_ind[count]      = Rcpp::as<Eigen::MatrixXi>(obs_tmp["pred_ind"]);
    obs_ind[count]       = Rcpp::as<Eigen::MatrixXi>(obs_tmp["obs_ind"]);
    Ys[count]            = Rcpp::as<Eigen::VectorXd>(obs_tmp["Y"]);
    
    if(use_random_effect == 1)
    {
      Brandom_pred[count]  = Rcpp::as<Eigen::MatrixXd>(obs_tmp["Brandom_pred"]);

      if(Brandom_pred[count].cols() != mixobj->beta_random.size())
      {
        Rcpp::Rcout << "predict: error dimension missmatch :\n Brandom[" << count << "].cols() = " << Brandom_pred[count].cols() << " != "
                   << "mixobj->beta_random.size() = "<<  mixobj->beta_random.size()  << "\n";
        throw("input error\n");
      }
    }
    Bfixed_pred[count]   = Rcpp::as<Eigen::MatrixXd>(obs_tmp["Bfixed_pred"]);

    if(Bfixed_pred[count].cols() != mixobj->beta_fixed.size())
    {
      Rcpp::Rcout << "predict: error dimension missmatch :\n Bfixed[" << count << "].cols() = " << Bfixed_pred[count].cols() << " != "
                  << "mixobj->beta_fixed.size() = "<<  mixobj->beta_fixed.size()  << "\n";
      throw("input error\n");

    }
    count++;
  }

  //**********************************
  //operator setup
  //***********************************
  Rcpp::List operator_list;
  std::string type_operator;

  std::vector<  solver* >  Solver;
  operatorMatrix *Kobj;

  if(use_process == 1){

    if(silent == 0){
      Rcpp::Rcout << " Setup operator\n";
    }

    operator_list  = Rcpp::as<Rcpp::List> (in_list["operator_list"]);
    type_operator = Rcpp::as<std::string>(operator_list["type"]);
    operator_list["nIter"] = 1;

    operator_select(type_operator, &Kobj);

    Kobj->initFromList(operator_list, List::create(Rcpp::Named("use.chol") = 1));

    Eigen::SparseMatrix<double, 0, int> Q, K;

    count = 0;
    for( List::iterator it = obs_list.begin(); it != obs_list.end(); ++it ) {
      List obs_tmp = Rcpp::as<Rcpp::List>( *it);
      //Do not neet LU solver here since we only work with Q. 
      //if(type_operator == "bivariate matern"){
      //  Solver.push_back( new lu_solver );  
      //} else {
      //  Solver.push_back( new cholesky_solver );  
      //}
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
  }





  /*
  Simulation objects
  */
  const int nP = 1;
  std::normal_distribution<double> normal;

  std::vector<std::mt19937> random_engine(nP);
  std::vector<gig> rgig(nP);
  for (int i = 0; i < nP; ++i) {
    if(in_list.containsElementNamed("seed")){
      random_engine[i].seed(seed);
      rgig[i].seed(seed);
    }else{
      random_engine[i].seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());
      rgig[i].seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    }
  }

  std::vector<int> longInd;
  for (int i=0; i< nindv; i++){
    longInd.push_back(i);
  }


  std::vector< Eigen::MatrixXd > WVec(nindv); //process
  std::vector< Eigen::MatrixXd > VVec(nindv); //variance of process
  std::vector< Eigen::MatrixXd > XVec(nindv); //latent field (process + effects)
  std::vector< Eigen::MatrixXd > WnoiseVec(nindv); //process noise
  std::vector< Eigen::MatrixXd > VnoiseVec(nindv); //process variances
  std::vector< Eigen::MatrixXd > UVec(nindv); //random effects
  std::vector< Eigen::MatrixXd > XVec_deriv(nindv);
  std::vector< Eigen::MatrixXd > WVec_deriv(nindv);
  std::vector< Eigen::MatrixXd > YVec(nindv);

  /**
   * Timing setup
   */
  const double ticks_per_ms = static_cast<double>(CLOCKS_PER_SEC);


  int rank = 0;
  double percent_done = 0;

  for(int i = 0; i < nindv; i++ ) {

    clock_t start = clock();
    percent_done++;
    double start_time;
    if(timeit)
      start_time = static_cast<double>(clock())/ ticks_per_ms;
    if(silent == 0){
      std::stringstream stream;
      stream << " Prediction " << 100*percent_done/nindv << " % done\n";
      Rcpp::Rcout << stream.str();
    }
    XVec[i].setZero(Bfixed_pred[i].rows(), nSim);
    WVec[i].setZero(Bfixed_pred[i].rows(), nSim);
    if(use_process){
      if(Kobj->nop==1){
        WnoiseVec[i].setZero(Kobj->d[0], nSim);  
        VnoiseVec[i].setZero(Kobj->d[0], nSim);  
      } else {
        WnoiseVec[i].setZero(Kobj->d[i], nSim);  
        VnoiseVec[i].setZero(Kobj->d[i], nSim);   
      }
    }
    
    
    if(predict_derivative == 1){
      XVec_deriv[i].setZero(Bfixed_pred[i].rows(), nSim);
      WVec_deriv[i].setZero(Bfixed_pred[i].rows(), nSim);
    }
    YVec[i].setZero(Bfixed_pred[i].rows(), nSim);
    VVec[i].setZero(Bfixed_pred[i].rows(), nSim);
    Eigen::MatrixXd random_effect;
    if(use_random_effect == 1){
      random_effect = mixobj->Br[i];
      UVec[i].setZero(Brandom_pred[i].cols(), nSim);
    }
    Eigen::MatrixXd fixed_effect = mixobj->Bf[i];

    for(int ipred = 0; ipred < pred_ind[i].rows(); ipred++){
      int n_obs;
      if(ind_general){
        n_obs= obs_ind[i].row(ipred).sum();
      } else {
        n_obs= obs_ind[i](ipred,1) - obs_ind[i](ipred,0);  
      }
      
      Eigen::VectorXd Y = Ys[i];
      Eigen::SparseMatrix<double,0,int> A = As[i];
      if(n_obs > 0){
        if(use_process == 1){
          if(ind_general){
            VectorXi oi = obs_ind[i].row(ipred);
            get_submatrix(As[i],oi, &A);
          } else {
            A = As[i].middleRows(obs_ind[i](ipred,0),obs_ind[i](ipred,1));
          }
        }
        if(ind_general){
          get_subvector(Ys[i],obs_ind[i].row(ipred), Y);
        } else {
          Y = Ys[i].segment(obs_ind[i](ipred,0),obs_ind[i](ipred,1));
        }
        if(use_random_effect == 1){
          if(ind_general){
             get_submatrix(random_effect,obs_ind[i].row(ipred), mixobj->Br[i]);
          } else {
            mixobj->Br[i] = random_effect.middleRows(obs_ind[i](ipred,0),obs_ind[i](ipred,1));  
          }
        }
        if(ind_general){
          get_submatrix(fixed_effect,obs_ind[i].row(ipred), mixobj->Bf[i]);
        } else {
          mixobj->Bf[i] = fixed_effect.middleRows(obs_ind[i](ipred,0),obs_ind[i](ipred,1));
        }
        
      }

    if(timeit){
      time_setup += static_cast<double>(clock())/ ticks_per_ms - start_time;
      start_time = static_cast<double>(clock())/ ticks_per_ms;
    }
      for(int ii = 0; ii < nSim + nBurnin; ii ++){
        if(debug){
          Rcpp::Rcout << "Prediction number = " << ipred << ", sample = "<< ii<< "\n";
        }
        Eigen::VectorXd res = Y;
        if(n_obs>0){
          mixobj->remove_cov(i, res);
          errObj->remove_asym(i, res);
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
              VectorXd Vi;
              if(ind_general){
                get_subvector(errObj->Vs[i],obs_ind[i].row(ipred), Vi);
              } else {
                Vi = errObj->Vs[i].segment(obs_ind[i](ipred,0),obs_ind[i](ipred,1));  
              }
              mixobj->sampleU2_par( i, res, Vi.cwiseInverse(),
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
            Rcpp::Rcout << "Compute operator \n";
          }
          double time_operator_temp = static_cast<double>(clock())/ ticks_per_ms;
          Eigen::SparseMatrix<double, 0, int> Qi;
          int d = 1;
          if(Kobj->nop == 1){
            Ki = Eigen::SparseMatrix<double,0,int>(Kobj->Q[0]);  
            d = Kobj->d[0];
          } else {
            Ki = Eigen::SparseMatrix<double,0,int>(Kobj->Q[i]);
            d = Kobj->d[i];
          }

          Eigen::VectorXd iV(process->Vs[i].size());
          iV.array() = process->Vs[i].array().inverse();
          //Rcpp::Rcout << iV.maxCoeff()  << "\n";
          Qi = Ki.transpose();

          if(Ki.cols() != iV.size()){
              Rcpp::Rcout << "prediction::the columns of Ki ( " << Ki.cols() << ") does not match iV length ( " <<  iV.size() << " )\n";
              throw("input error\n");
          }

          Qi =  Qi * iV.asDiagonal();
          Qi =  Qi * Ki;
          if(debug){
            Rcpp::Rcout << "Sample normals \n";
          }
          time_operator += static_cast<double>(clock())/ ticks_per_ms  - time_operator_temp;
    
          Eigen::VectorXd zi;
          zi.setZero(d);
          for(int j =0; j < d; j++)
            zi[j] =  normal(random_engine[rank]);
          if(n_obs>0)
            res += A * process->Xs[i];

          if(debug){
            Rcpp::Rcout << "Sample process \n";
          }
          double sample_p_temp = static_cast<double>(clock())/ ticks_per_ms;
          if(n_obs>0){
            if(type_MeasurementError == "Normal"){
              process->sample_X(i,
                                zi,
                                res,
                                Qi,
                                Ki,
                                A,
                                errObj->sigma,
                                *Solver[i]);
            } else {
              VectorXd Vi;
              if(ind_general){
                get_subvector(errObj->Vs[i],obs_ind[i].row(ipred), Vi);
              } else {
                Vi = errObj->Vs[i].segment(obs_ind[i](ipred,0),obs_ind[i](ipred,1));  
              }
              process->sample_Xv2( i,
                                   zi,
                                   res,
                                   Qi,
                                   Ki,
                                   A,
                                   errObj->sigma,
                                   *Solver[i],
                                   Vi.cwiseInverse());
            }
            res -= A * process->Xs[i];
            time_process_sample += static_cast<double>(clock())/ ticks_per_ms - sample_p_temp;
          } else {
            process->simulate(i,zi,A,Ki,Y,*Solver[i]);
            time_process_simulate += static_cast<double>(clock())/ ticks_per_ms - sample_p_temp;
          }
          
        }

        if(debug){
          Rcpp::Rcout << "Sample variances \n";
        }
        double sample_v_temp = static_cast<double>(clock())/ ticks_per_ms;
        if(use_process == 1){
          if(debug){
            Rcpp::Rcout << "nobs = " << n_obs << "\n";
          }
          if(n_obs>0){
              process->sample_V(i, rgig[rank], Ki);
            } else {
              process->simulate_V(i, rgig[rank]);
          }
        }
        time_v_sample += static_cast<double>(clock())/ ticks_per_ms - sample_v_temp;
        //***********************************
        // random variance noise sampling
        //***********************************
        errObj->add_asym(i, res);
        if(type_MeasurementError != "Normal" || type_MeasurementError != "nsNormal"){
          if(n_obs>0){
            if(ind_general){
              errObj->sampleV(i, res,-1);
            } else {
              errObj->sampleV(i, res,obs_ind[i](ipred,0)+obs_ind[i](ipred,1));  
            }
            
          }
        }
        errObj->remove_asym(i, res);
        if(timeit){
          time_sample += static_cast<double>(clock())/ ticks_per_ms - start_time;
          start_time   = static_cast<double>(clock())/ ticks_per_ms;
        }

          // save samples
          if(ii >= nBurnin){
            if(debug){
              Rcpp::Rcout << "Save samples\n";
            }
          if(use_process == 1){
            if(ind_general){
              get_submatrix(As_pred[i],pred_ind[i].row(ipred), &Ai);
            } else {
              Ai = As_pred[i].middleRows(pred_ind[i](ipred,0),pred_ind[i](ipred,1));
            }
          }
          Eigen::VectorXd random_effect;
          if(use_random_effect == 1){
            random_effect= Bfixed_pred[i]*mixobj->beta_fixed + Brandom_pred[i]*(mixobj->U.col(i)+mixobj->beta_random);
          } else {
            random_effect = Bfixed_pred[i]*mixobj->beta_fixed;
          }
          Eigen::VectorXd random_effect_c;
          Eigen::VectorXd mNoise; 
          Eigen::VectorXd Noise = errObj->simulate_par(i,random_engine[rank],random_effect.size());
          if(ind_general){
            get_subvector(random_effect,pred_ind[i].row(ipred), random_effect_c);
            get_subvector(Noise,pred_ind[i].row(ipred), mNoise);
          }else{
            random_effect_c = random_effect.segment(pred_ind[i](ipred,0),pred_ind[i](ipred,1));  
            mNoise = Noise.segment(pred_ind[i](ipred,0),pred_ind[i](ipred,1));  
          }
          
          Eigen::VectorXd AX;
          if(use_random_effect == 1){
            UVec[i].col(ii-nBurnin) = mixobj->U.col(i); // this only makes sense for smoothing
          }
          
          if(use_process == 1){
            AX = Ai * process->Xs[i];

            WnoiseVec[i].col(ii-nBurnin) = process->Ws[i]; // this only makes sense for smoothing
            VnoiseVec[i].col(ii-nBurnin) = process->Vs[i]; // this only makes sense for smoothing

            if(ind_general){
              //Rcpp::Rcout << "WVec pre\n" << WVec[i] << "\n AX \n" << AX << "\n pred_ind = " << pred_ind[i].row(ipred) << "\n";
              set_subcol(WVec[i], ii-nBurnin, pred_ind[i].row(ipred), AX);
              //Rcpp::Rcout << "WVec post\n" << WVec[i] << "\n";
              set_subcol(XVec[i], ii-nBurnin, pred_ind[i].row(ipred), random_effect_c + AX);
              set_subcol(YVec[i], ii-nBurnin, pred_ind[i].row(ipred), random_effect_c + AX + mNoise);
              set_subcol(VVec[i], ii-nBurnin, pred_ind[i].row(ipred), Ai * process->Vs[i]);
            } else {
              WVec[i].block(pred_ind[i](ipred,0), ii - nBurnin, pred_ind[i](ipred,1), 1) = AX;
              XVec[i].block(pred_ind[i](ipred,0), ii - nBurnin, pred_ind[i](ipred,1), 1) = random_effect_c + AX;
              YVec[i].block(pred_ind[i](ipred,0), ii - nBurnin, pred_ind[i](ipred,1), 1) = random_effect_c + AX + mNoise;
              VVec[i].block(pred_ind[i](ipred,0), ii - nBurnin, pred_ind[i](ipred,1), 1) = Ai * process->Vs[i]; 
            }
            
          }  else {
            if(ind_general){
              set_subcol(XVec[i], ii-nBurnin, pred_ind[i].row(ipred), random_effect_c);
              set_subcol(YVec[i], ii-nBurnin, pred_ind[i].row(ipred), random_effect_c + mNoise);
            } else {
              XVec[i].block(pred_ind[i](ipred,0), ii - nBurnin, pred_ind[i](ipred,1), 1) = random_effect_c;
              YVec[i].block(pred_ind[i](ipred,0), ii - nBurnin, pred_ind[i](ipred,1), 1) = random_effect_c + mNoise;  
            }
            
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
              Eigen::SparseMatrix<double,0,int> Ai_1;
              if(ind_general){
                get_submatrix(As_pred_1[i],pred_ind[i].row(ipred), &Ai_1);
              } else {
                Ai_1 = As_pred_1[i].middleRows(pred_ind[i](ipred,0),pred_ind[i](ipred,1));  
              }
              
              Eigen::VectorXd AX_1 = Ai_1 * process->Xs[i];
              
              if(ind_general){
                set_subcol(WVec_deriv[i], ii-nBurnin, pred_ind[i].row(ipred), AX_1 - AX);
                set_subcol(XVec_deriv[i], ii-nBurnin, pred_ind[i].row(ipred), random_effect_c_1 + AX_1 - random_effect_c - AX);
              } else {
                WVec_deriv[i].block(pred_ind[i](ipred,0), ii - nBurnin, pred_ind[i](ipred,1), 1) = AX_1 - AX;
                XVec_deriv[i].block(pred_ind[i](ipred,0), ii - nBurnin, pred_ind[i](ipred,1), 1) = random_effect_c_1 + AX_1 - random_effect_c - AX;  
              }
            } else {
              if(ind_general){
                set_subcol(XVec_deriv[i], ii-nBurnin, pred_ind[i].row(ipred), random_effect_c_1 - random_effect_c);
              } else {
                XVec_deriv[i].block(pred_ind[i](ipred,0), ii - nBurnin, pred_ind[i](ipred,1), 1) = random_effect_c_1 - random_effect_c;  
              }
              
            }
          }
        }
        if(timeit){
            time_post += static_cast<double>(clock())/ ticks_per_ms - start_time;
            start_time   = static_cast<double>(clock())/ ticks_per_ms;
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
      stream << "time = " << time_Ma;
      Rcpp::Rcout << stream.str();
      Rcpp::Rcout << "\n";
    }
    if(timeit){
      Rcpp::Rcout << "setup = " << time_setup << "\n";
      Rcpp::Rcout << "sample = " << time_sample << "\n";
      if(use_process == 1){
        Rcpp::Rcout << "\tsample_process = " << time_process_sample << "\n";
        Rcpp::Rcout << "\tsimulate_process = " << time_process_simulate << "\n";
        Rcpp::Rcout << "\ttime_v_sample = " << time_v_sample << "\n";
        Rcpp::Rcout << "\ttime_operator = " << time_operator << "\n";
        time_v_sample = 0;
        time_process_simulate = 0;
        time_process_sample =0;
        time_operator = 0;
      }
      Rcpp::Rcout << "post = " << time_post << "\n";
      time_setup = 0;
      time_sample = 0;
      time_post = 0;

       
    }
  }
  
  if(debug){
    Rcpp::Rcout << "store results\n";  
  }
  
  // storing the results
  Rcpp::List out_list;
  out_list["YVec"] = YVec;
  out_list["XVec"] = XVec;
  out_list["WVec"] = WVec;
  if(use_process){
    out_list["WnoiseVec"] = WnoiseVec;
    out_list["VnoiseVec"] = VnoiseVec;
  }
  out_list["VVec"] = VVec;
  out_list["UVec"] = UVec;
  if(predict_derivative == 1)
  {
    out_list["XVec_deriv"] = XVec_deriv;
    out_list["WVec_deriv"] = WVec_deriv;
  }
  //Free memory
  delete mixobj;
  delete errObj;
  if(use_process == 1){
    delete process;
    delete Kobj;
  }
  return(out_list);
}
