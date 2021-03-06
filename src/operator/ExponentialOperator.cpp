#include "operatorMatrix.h"
#include "error_check.h"
#include "eigen_add_on.h"
#include "operator_helper.h"

ExponentialOperator::~ExponentialOperator()
{
  if(is_initialized == 1){
    delete[] Q;
    delete[] G;
    delete[] C;
    delete[] d2tauQ;
    delete[] dtauQ;
    delete[] dkappaQ;
    delete[] d2kappaQ;
    for(int i=0;i<nop;i++){
      delete Qsolver[i];
      delete Qepssolver[i];
    }
    delete[] Qsolver;
    delete[] Qepssolver;
  }
}

void ExponentialOperator::initFromList(Rcpp::List const & init_list, Rcpp::List const & solver_list)
{
  is_initialized = 1;
  npars = 2;
  std::vector<std::string> check_names =  {"C", "G", "kappa", "tau","h","beta"};
  check_Rcpplist(init_list, check_names, "ExponentialOperator::initFromList");
  std::vector<std::string> check_names2 =  {"use.chol"};
  check_Rcpplist(solver_list, check_names2, "ExponentialOperator::initFromList");
  out_list = clone(init_list);
  kappa = Rcpp::as<double>(init_list["kappa"]);
  beta = Rcpp::as<double>(init_list["beta"]);
  Rcpp::Rcout << "beta = " << beta << "\n";
  if(kappa < 0){
    Rcpp::Rcout << "warning kappa negative\n";
  }
  tau = Rcpp::as<double>(init_list["tau"]);
  int nIter = Rcpp::as<double>(init_list["nIter"]);
  kappaVec.resize(nIter);
  tauVec.resize(nIter);
  dkappa = 0;
  dkappa_old = 0;
  dtau = 0;
  dtau_old = 0;
  use_chol = Rcpp::as<int>(solver_list["use.chol"]);
  
  Rcpp::List G_list  = Rcpp::as<Rcpp::List> (init_list["G"]);
  Rcpp::List C_list  = Rcpp::as<Rcpp::List> (init_list["C"]);
  Rcpp::List Ci_list  = Rcpp::as<Rcpp::List> (init_list["Ci"]);
  Rcpp::List h_list  = Rcpp::as<Rcpp::List> (init_list["h"]);
  nop = G_list.size();
  
  d.resize(nop);
  loc.resize(nop);
  h.resize(nop);
  h_average.resize(nop);
  matrix_set.resize(nop);
  tau_trace.resize(nop);
  tau_trace2.resize(nop);
  kappa_trace.resize(nop);
  kappa_trace2.resize(nop);
  
  Q = new Eigen::SparseMatrix<double,0,int>[nop];
  G = new Eigen::SparseMatrix<double,0,int>[nop];
  GCG = new Eigen::SparseMatrix<double,0,int>[nop];
  C = new Eigen::SparseMatrix<double,0,int>[nop];
  d2tauQ = new Eigen::SparseMatrix<double,0,int>[nop];
  dtauQ = new Eigen::SparseMatrix<double,0,int>[nop];
  dkappaQ = new Eigen::SparseMatrix<double,0,int>[nop];
  d2kappaQ = new Eigen::SparseMatrix<double,0,int>[nop];
  
  Qsolver = new solver*[nop];
  Qepssolver = new solver*[nop];
  
  for(int i=0;i<nop;i++){
    //Rcpp::Rcout << "init patient " << i << "\n";
    matrix_set[i] = 0;
    if(use_chol==1){
      Qsolver[i] = new lu_sparse_solver;
      Qepssolver[i] = new lu_sparse_solver;
    } else {
      Qsolver[i] = new iterative_solver;
      Qepssolver[i] = new iterative_solver;
    }
    
    G[i] =  Rcpp::as<Eigen::SparseMatrix<double,0,int>>(G_list[i]);
    C[i] =  Rcpp::as<Eigen::SparseMatrix<double,0,int>>(C_list[i]);
    Eigen::SparseMatrix<double,0,int> Ci =  Rcpp::as<Eigen::SparseMatrix<double,0,int>>(Ci_list[i]);
    GCG[i] =  G[i]*Ci*G[i];
    if(beta==2){
      Q[i] = G[i] + C[i] + GCG[i];  
    } else {
      Q[i] = G[i] + C[i];
    }
    
    d[i] = Q[i].rows();
    //d2tauQ.resize(n[i],n[i]);
    
    h[i]  = Rcpp::as< Eigen::VectorXd >(h_list[i]);
    
    h_average[i] = h[i].sum() / h[i].size();
    (*Qsolver[i]).initFromList(d[i],solver_list);
    (*Qsolver[i]).analyze(Q[i]);
    
    (*Qepssolver[i]).initFromList(d[i],solver_list);
    (*Qepssolver[i]).analyze(Q[i]);
  }
  
  //Rcpp::Rcout << "set matrices\n";
  this->set_matrices();
  //Rcpp::Rcout << "init done\n";
}


void ExponentialOperator::set_matrices()
{
  for(int i=0;i<nop;i++){
    this->set_matrix(i);
  }
}


void ExponentialOperator::set_matrix(int i)
{
  if(matrix_set[i]==0){
    matrix_set[i] = 1;
    
    double c = pow(2,-0.5);
    if(beta==2){
      Q[i] =        c*tau*(pow(kappa,2)*C[i] + 2*kappa*G[i] + GCG[i]);
      dkappaQ[i] =  c*tau*(2*kappa*C[i] + 2*G[i]);
      d2kappaQ[i] = c*tau*(2*C[i]);
      dtauQ[i] =    c*(pow(kappa,2)*C[i] + 2*kappa*G[i] + GCG[i]);
      
      Q[i] =        c*tau*(       pow(kappa,0.5)*C[i] +   2*pow(kappa,-0.5)*G[i] +     pow(kappa,-1.5)*GCG[i]);
      dkappaQ[i] =  c*tau*(  0.5*pow(kappa,-0.5)*C[i] -     pow(kappa,-1.5)*G[i] - 1.5*pow(kappa,-2.5)*GCG[i]);
      d2kappaQ[i] = c*tau*(-0.25*pow(kappa,-1.5)*C[i] + 1.5*pow(kappa,-2.5)*G[i] + 3.75*pow(kappa,-3.5)*GCG[i]);
      dtauQ[i] =    c*(pow(kappa,0.5)*C[i] + 2*pow(kappa,-0.5)*G[i] + pow(kappa,-1.5)*GCG[i]);
      
    } else {
      Q[i] =        c*tau*(     pow(kappa,-0.5)*G[i] +     pow(kappa,0.5)*C[i]);
      dkappaQ[i] =  c*tau*(-0.5*pow(kappa,-1.5)*G[i] + 0.5*pow(kappa,-0.5)*C[i]);
      d2kappaQ[i] = c*tau*(0.75*pow(kappa,-3.5)*G[i] - 0.25*pow(kappa,-1.5)*C[i]);
      dtauQ[i] =    c*(pow(kappa,-0.5)*G[i] + pow(kappa,0.5)*C[i]);  
    }
    
    (*Qsolver[i]).compute(Q[i]);
    tau_trace[i] = Q[i].cols()/tau;
    tau_trace2[i] = -tau_trace[i]/tau;
    
    //To compute the traces for kappa, we use that :
    //trace(dK*K^-1) = sum((diag(dK)./diag(K))) for triangular matrices
    //trace(dK*K^-1*dK*K^-1) = sum((diag(dK)./diag(K)).^2) for triangular matrices
    VectorXd Kd = dkappaQ[i].diagonal();
    Kd = Kd.cwiseQuotient(Q[i].diagonal());
    kappa_trace[i] = Kd.sum();
    
    VectorXd Kd2 = d2kappaQ[i].diagonal();
    Kd2 = Kd2.cwiseQuotient(Q[i].diagonal());
    kappa_trace2[i] = Kd2.sum();
    Kd = Kd.array().square();
    kappa_trace2[i] -= Kd.sum();
  }
}

Rcpp::List ExponentialOperator::output_list()
{
  out_list["tau"] = tauVec[tauVec.size() - 1];
  out_list["kappa"] = kappaVec[tauVec.size() - 1];
  out_list["tauVec"] = tauVec;
  out_list["kappaVec"] = kappaVec;
  out_list["nIter"] = tauVec.size();
  out_list["use.chol"] = use_chol;
  out_list["Cov_theta"]   = Cov_theta;
  return(out_list);
}


void ExponentialOperator::gradient_init(int nsim, int nrep)
{
  dtau =  0;
  ddtau = 0;
  dkappa = 0;
  ddkappa = 0;
  ddtaukappa = 0;
}
Eigen::MatrixXd ExponentialOperator::d2Given( const Eigen::VectorXd & X,
                                              const Eigen::VectorXd & iV,
                                              const Eigen::VectorXd & mean_KX,
                                              int ii,
                                              const double weight)
{
  int i = ii;
  if(nop == 1)
    i = 0;
  this->set_matrix(i);
  Eigen::VectorXd vtmp = Q[i] * X;
  
  
  Eigen::VectorXd KX = Q[i]*X;
  //compute gradients wrt tau
  Eigen::VectorXd dKX = dtauQ[i] * X;
  Eigen::VectorXd d2KX;
  
  Eigen::MatrixXd d2 = Eigen::MatrixXd::Zero(2, 2);
  
  d2(0, 0)  =- weight * tau_trace2[i];
  d2(0, 0) -=- weight * (dKX.dot(iV.asDiagonal() * dKX));
  dKX      = dkappaQ[i] * X;
  d2KX     = d2kappaQ[i] * X;
  d2(0, 1) -=- weight * dKX.dot(iV.asDiagonal() * (2 * KX - mean_KX))/tau;
  d2(1, 0) = d2(0, 1);
  d2(1, 1)  =- weight * kappa_trace2[i];
  d2(1, 1) -=- weight * dKX.dot(iV.asDiagonal()*dKX);
  d2(1, 1) -=- weight * d2KX.dot(iV.asDiagonal()*(KX - mean_KX));
  
  return(d2);
  
}

void ExponentialOperator::gradient_add( const Eigen::VectorXd & X,
                                        const Eigen::VectorXd & iV,
                                        const Eigen::VectorXd & mean_KX,
                                        const int ii,
                                        const double weight)
{
  int i = ii;
  if(nop == 1)
    i = 0;
  
  this->set_matrix(i);
  dtau    += weight * tau_trace[i];
  ddtau   += weight * tau_trace2[i];
  dkappa  += weight * kappa_trace[i];
  ddkappa += weight * kappa_trace2[i];
  
  Eigen::VectorXd KX = Q[i]*X;
  //compute gradients wrt tau
  Eigen::VectorXd dKX = dtauQ[i] * X;
  Eigen::VectorXd d2KX;
  //Eigen::VectorXd d2KX = d2tauQ * X;
  dtau -=  weight * dKX.dot(iV.asDiagonal() * (KX - mean_KX));
  ddtau -= weight * (dKX.dot(iV.asDiagonal() * dKX));// + d2KX.dot(iV.asDiagonal() * (KX - mean_KX)));
  
  //compute gradients wrt kappa
  dKX      = dkappaQ[i] * X;
  d2KX     = d2kappaQ[i] * X;
  dkappa  -= weight * dKX.dot(iV.asDiagonal() * (KX - mean_KX));
  ddkappa -= weight * dKX.dot(iV.asDiagonal()*dKX);
  ddkappa -= weight * d2KX.dot(iV.asDiagonal()*(KX - mean_KX));
  ddtaukappa -= weight * dKX.dot(iV.asDiagonal()*(2*KX - mean_KX))/tau;

}

void ExponentialOperator::gradient( const Eigen::VectorXd & X, const Eigen::VectorXd & iV)
{
  throw(" ExponentialOperator::gradient depricated \n");
}

void ExponentialOperator::print_parameters(){
  Rcpp::Rcout << "tau = " << tau << "\n";
  Rcpp::Rcout << "kappa = " << kappa << "\n";
}


void ExponentialOperator::step_theta(const double stepsize,
                                     const double learning_rate,
                                     const double polyak_rate,
                                     const int burnin)
{
  if(1){ //independent steps
    dtau  /= ddtau;
    dtau_old = learning_rate * dtau_old + dtau;
    double step = stepsize * dtau_old;
    double tau_temp = -1.;
    while(tau_temp < 0)
    {
      step *= 0.5;
      tau_temp = tau - step;
    }
    tau = tau_temp;
    
    dkappa  /= ddkappa;
    dkappa_old = learning_rate * dkappa_old + dkappa;
    step   = stepsize * dkappa_old;
    double kappa_temp = -1.;
    while(kappa_temp < 0)
    {
      step *= 0.5;
      kappa_temp = kappa - step;
    }
    kappa = kappa_temp;  
  } else { //better correlated step
    Eigen::MatrixXd d2 = Eigen::MatrixXd::Zero(2, 2);
    Eigen::VectorXd step_v, dtheta_v, dtheta_old_v, theta_v, theta_v_tmp;
    
    dtheta_v.setZero(2);
    dtheta_old_v.setZero(2);
    theta_v.setZero(2);
    theta_v_tmp.setZero(2);
    
    theta_v(0) = tau;
    theta_v(1) = kappa;
    
    dtheta_v(0) = dtau;
    dtheta_v(1) = dkappa;
    
    dtheta_old_v(0) = dtau_old;
    dtheta_old_v(1) = dkappa_old;
    
    d2(0, 0)  = ddtau;
    d2(1, 0) = ddtaukappa;
    d2(0,1) = ddtaukappa;
    d2(1,1) = ddkappa;
    
    step_v  = d2.ldlt().solve(dtheta_v);
    
    theta_v_tmp(0) = -1;
    theta_v_tmp(1) = -1;
    
    while(theta_v_tmp(0) < 0 || theta_v_tmp(1) < 0)
    {
      step_v *= 0.5;
      theta_v_tmp = theta_v - step_v;
    }
    theta_v = theta_v_tmp;
    tau = theta_v(0);
    kappa = theta_v(1);
  }
  
  
  if(counter == 0 || polyak_rate == -1)
    tauVec[counter] = tau;
  else
    tauVec[counter] = polyak_rate * tau + (1 - polyak_rate) * tauVec[counter-1];
  
  if(counter == 0 || polyak_rate == -1)
    kappaVec[counter] = kappa;
  else
    kappaVec[counter] = polyak_rate * kappa + (1 - polyak_rate) * kappaVec[counter-1];
  counter++;
  clear_gradient();
  ddtau   = 0;
  ddkappa = 0;
  for(int i=0;i<nop;i++){
    matrix_set[i] = 0;
  }
  
  //this->set_matrices();
}

double ExponentialOperator::trace_variance( const Eigen::SparseMatrix<double,0,int> & A, int i)
{
  //return(A.rows() * tau/ h_average);
  return(-1);
}

Eigen::VectorXd  ExponentialOperator::get_gradient()
{
  Eigen::VectorXd g(npars);
  g[0] = dtau;
  g[1] = dkappa;
  return(g);
}
void  ExponentialOperator::clear_gradient()
{
  dkappa = 0;
  dtau   = 0;
};
