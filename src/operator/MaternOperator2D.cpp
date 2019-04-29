#include "operatorMatrix.h"
#include "error_check.h"
#include "eigen_add_on.h"
#include "operator_helper.h"
#include "MatrixAlgebra.h"

MaternOperator2D::~MaternOperator2D()
{
  if(is_initialized == 1){
    delete[] Q;
    delete[] G;
    delete[] C;
    delete[] d2tau1Q;
    delete[] d2tau2Q;
    delete[] dtau1Q;
    delete[] dtau2Q;
    delete[] dkappa1Q;
    delete[] dkappa2Q;
    delete[] d2kappa1Q;
    delete[] d2kappa2Q;
    delete[] d2rhoQ;
    delete[] drhoQ;
    delete[] dthetaQ;
    delete[] d2thetaQ;
    for(int i=0;i<nop;i++){
      delete Qsolver[i];
      delete Qepssolver[i];
    }
    delete[] Qsolver;
    delete[] Qepssolver;
  }
}

void MaternOperator2D::initFromList(Rcpp::List const & init_list, Rcpp::List const & solver_list)
{
  is_initialized = 1;
  npars = 6;
  std::vector<std::string> check_names =  {"C", "G", "kappa1", "kappa2", "tau1", "tau2","rho","theta","h"};
  check_Rcpplist(init_list, check_names, "MaternOperator2D::initFromList");
  std::vector<std::string> check_names2 =  {"use.chol"};
  check_Rcpplist(solver_list, check_names2, "MaternOperator2D::initFromList");
  out_list = clone(init_list);
  kappa1 = Rcpp::as<double>(init_list["kappa1"]);
  kappa2 = Rcpp::as<double>(init_list["kappa2"]);
  if(kappa1 < 0 || kappa2 < 0){
    Rcpp::Rcout << "warning kappa negative\n";
  }
  tau1 = Rcpp::as<double>(init_list["tau1"]);
  tau2 = Rcpp::as<double>(init_list["tau2"]);
  if(tau1 < 0 || tau2 < 0){
    Rcpp::Rcout << "warning tau negative\n";
  }
  rho = Rcpp::as<double>(init_list["rho"]);
  theta = Rcpp::as<double>(init_list["theta"]);
  estimate_theta = Rcpp::as<int>(init_list["estimate_theta"]);
  
  int nIter = Rcpp::as<double>(init_list["nIter"]);
  tau1Vec.resize(nIter);
  tau2Vec.resize(nIter);
  kappa1Vec.resize(nIter);
  kappa2Vec.resize(nIter);
  rhoVec.resize(nIter);
  thetaVec.resize(nIter);
  dkappa1 = 0;
  dkappa1_old = 0;
  dtau1 = 0;
  dtau1_old = 0;
  dkappa2 = 0;
  dkappa2_old = 0;
  dtau2 = 0;
  dtau2_old = 0;
  drho = 0;
  drho_old = 0;
  dtheta = 0;
  dtheta_old = 0;
  
  use_chol = Rcpp::as<int>(solver_list["use.chol"]);
  
  Rcpp::List G_list  = Rcpp::as<Rcpp::List> (init_list["G"]);
  Rcpp::List C_list  = Rcpp::as<Rcpp::List> (init_list["C"]);
  Rcpp::List h_list  = Rcpp::as<Rcpp::List> (init_list["h"]);
  nop = G_list.size();
  
  d.resize(nop);
  loc.resize(nop);
  h.resize(nop);
  h_average.resize(nop);
  matrix_set.resize(nop);
  
  tau1_trace.resize(nop);
  tau1_trace2.resize(nop);
  kappa1_trace.resize(nop);
  kappa1_trace2.resize(nop);
  
  tau2_trace.resize(nop);
  tau2_trace2.resize(nop);
  kappa2_trace.resize(nop);
  kappa2_trace2.resize(nop);
  
  rho_trace.resize(nop);
  rho_trace2.resize(nop);
  
  theta_trace.resize(nop);
  theta_trace2.resize(nop);  
  
  
  
  Q = new Eigen::SparseMatrix<double,0,int>[nop];
  G = new Eigen::SparseMatrix<double,0,int>[nop];
  C = new Eigen::SparseMatrix<double,0,int>[nop];
  
  d2tau1Q = new Eigen::SparseMatrix<double,0,int>[nop];
  dtau1Q = new Eigen::SparseMatrix<double,0,int>[nop];
  dkappa1Q = new Eigen::SparseMatrix<double,0,int>[nop];
  d2kappa1Q = new Eigen::SparseMatrix<double,0,int>[nop];
  
  d2tau2Q = new Eigen::SparseMatrix<double,0,int>[nop];
  dtau2Q = new Eigen::SparseMatrix<double,0,int>[nop];
  dkappa2Q = new Eigen::SparseMatrix<double,0,int>[nop];
  d2kappa2Q = new Eigen::SparseMatrix<double,0,int>[nop];
  
  d2rhoQ = new Eigen::SparseMatrix<double,0,int>[nop];
  drhoQ = new Eigen::SparseMatrix<double,0,int>[nop];
  if(estimate_theta){
    dthetaQ = new Eigen::SparseMatrix<double,0,int>[nop];
    d2thetaQ = new Eigen::SparseMatrix<double,0,int>[nop];
  }
  
  Qsolver = new solver*[nop];
  Qepssolver = new solver*[nop];
  
  for(int i=0;i<nop;i++){
    //Rcpp::Rcout << "init patient " << i << "\n";
    matrix_set[i] = 0;
    if(use_chol==1){
      Qsolver[i] = new cholesky_solver;
      Qepssolver[i] = new cholesky_solver;
    } else {
      Qsolver[i] = new iterative_solver;
      Qepssolver[i] = new iterative_solver;
    }
    
    G[i] =  Rcpp::as<Eigen::SparseMatrix<double,0,int>>(G_list[i]);
    C[i] =  Rcpp::as<Eigen::SparseMatrix<double,0,int>>(C_list[i]);
    
    //initialize Q to [G G;G G]
    d[i] = 2*G[i].rows();
    Q[i].resize(d[i],d[i]);
    setSparseBlock(&Q[i],0,0, G[i]);
    setSparseBlock(&Q[i],d[i]/2,d[i]/2, G[i]);
    setSparseBlock(&Q[i],0,d[i]/2, G[i]);
    setSparseBlock(&Q[i],d[i]/2,0, G[i]);
    
    //initialize derivative matrices
    dtau1Q[i].resize(d[i],d[i]);
    setSparseBlock(&dtau1Q[i],0,0, G[i]);
    setSparseBlock(&dtau1Q[i],d[i]/2,0, G[i]);
    
    dtau2Q[i].resize(d[i],d[i]);
    setSparseBlock(&dtau2Q[i],d[i]/2,d[i]/2, G[i]);
    setSparseBlock(&dtau2Q[i],0,d[i]/2, G[i]);
    
    dkappa1Q[i].resize(d[i],d[i]);
    setSparseBlock(&dkappa1Q[i],0,0, G[i]);
    setSparseBlock(&dkappa1Q[i],d[i]/2,0, G[i]);
    
    dkappa2Q[i].resize(d[i],d[i]);
    setSparseBlock(&dkappa2Q[i],d[i]/2,d[i]/2, G[i]);
    setSparseBlock(&dkappa2Q[i],0,d[i]/2, G[i]);
    
    
    drhoQ[i].resize(d[i],d[i]);
    setSparseBlock(&drhoQ[i],0,0, G[i]);
    setSparseBlock(&drhoQ[i],d[i]/2,d[i]/2, G[i]);
    setSparseBlock(&drhoQ[i],0,d[i]/2, G[i]);
    setSparseBlock(&drhoQ[i],d[i]/2,0, G[i]);
    
    if(estimate_theta == 1){
      dthetaQ[i].resize(d[i],d[i]);
      setSparseBlock(&dthetaQ[i],0,0, G[i]);
      setSparseBlock(&dthetaQ[i],d[i]/2,d[i]/2, G[i]);
      setSparseBlock(&dthetaQ[i],0,d[i]/2, G[i]);
      setSparseBlock(&dthetaQ[i],d[i]/2,0, G[i]);
    }
    
    //initialize second derivatives
    d2tau1Q[i].resize(d[i],d[i]);
    d2tau2Q[i].resize(d[i],d[i]);
    
    d2kappa1Q[i].resize(d[i],d[i]);
    setSparseBlock(&d2kappa1Q[i],0,0, G[i]);
    setSparseBlock(&d2kappa1Q[i],d[i]/2,0, G[i]);
    
    d2kappa2Q[i].resize(d[i],d[i]);                        
    setSparseBlock(&d2kappa2Q[i],d[i]/2,d[i]/2, G[i]);
    setSparseBlock(&d2kappa2Q[i],0,d[i]/2, G[i]);
    
    d2rhoQ[i].resize(d[i],d[i]);
    setSparseBlock(&d2rhoQ[i],d[i]/2,d[i]/2, G[i]);
    setSparseBlock(&d2rhoQ[i],0,d[i]/2, G[i]);
    
    if(estimate_theta == 1){
      d2thetaQ[i].resize(d[i],d[i]);
      setSparseBlock(&d2thetaQ[i],0,0, G[i]);
      setSparseBlock(&d2thetaQ[i],d[i]/2,d[i]/2, G[i]);
      setSparseBlock(&d2thetaQ[i],0,d[i]/2, G[i]);
      setSparseBlock(&d2thetaQ[i],d[i]/2,0, G[i]);
    }
    
    h[i]  = Rcpp::as< Eigen::VectorXd >(h_list[i]);
    
    h_average[i] = h[i].sum() / h[i].size();
    
    //Qsolver will only operate on blocks
    (*Qsolver[i]).initFromList(d[i]/2,solver_list);
    (*Qsolver[i]).analyze(G[i]);
    
    (*Qepssolver[i]).initFromList(d[i]/2,solver_list);
    (*Qepssolver[i]).analyze(G[i]);
  }
  
  //Rcpp::Rcout << "set matrices\n";
  this->set_matrices();
  //Rcpp::Rcout << "init done\n";
}


void MaternOperator2D::set_matrices()
{
  for(int i=0;i<nop;i++){
    this->set_matrix(i);
  }
}


void MaternOperator2D::set_matrix(int i)
{
  if(matrix_set[i]==0){
    matrix_set[i] = 1;
    
    //Set Q[i] from parameters by filling the four blocks:
    MatrixXd D(2,2);
    D(0,0) = cos(theta) + rho*sin(theta);
    D(0,1) = -sin(theta)*pow(1+pow(rho,2),0.5);
    D(1,0) = sin(theta) - rho*cos(theta);
    D(1,1) = cos(theta)*pow(1+pow(rho,2),0.5);
    
    MatrixXd Drho(2,2);
    Drho(0,0) = sin(theta);
    Drho(0,1) = -sin(theta)*rho*pow(1+pow(rho,2),-0.5);
    Drho(1,0) = -cos(theta);
    Drho(1,1) = cos(theta)*rho*pow(1+pow(rho,2),-0.5);
    
    MatrixXd Drho2(2,2);
    Drho2(0,0) = 0;
    Drho2(0,1) = -sin(theta)*pow(1+pow(rho,2),-0.5);
    Drho2(0,1)+= sin(theta)*pow(rho,2)*pow(1+pow(rho,2),-1.5);
    Drho2(1,0) = 0;
    Drho2(1,1) = cos(theta)*pow(1+pow(rho,2),-0.5);
    Drho2(1,1) -= cos(theta)*pow(rho,2)*pow(1+pow(rho,2),-1.5);
    
    MatrixXd Dtheta(2,2);
    MatrixXd Dtheta2(2,2);
    if(estimate_theta == 1){
      Dtheta(0,0) = -sin(theta) + rho*cos(theta);
      Dtheta(0,1) = -cos(theta)*pow(1+pow(rho,2),0.5);
      Dtheta(1,0) = cos(theta) + rho*sin(theta);
      Dtheta(1,1) = -sin(theta)*pow(1+pow(rho,2),0.5);
      
      Dtheta2(0,0) = -cos(theta) - rho*sin(theta);
      Dtheta2(0,1) = sin(theta)*pow(1+pow(rho,2),0.5);
      Dtheta2(1,0) = -sin(theta) + rho*cos(theta);
      Dtheta2(1,1) = -cos(theta)*pow(1+pow(rho,2),0.5);
    }
    
    SparseMatrix<double,0,int> B,K1,K2;
    K1 = (tau1/kappa1)*G[i];
    K1+=tau1*kappa1*C[i];
    K2 = (tau2/kappa2)*G[i];
    K2+=tau2*kappa2*C[i];
    
    B = D(0,0)*K1;
    setSparseBlock_update(&Q[i],0,0, B);
    B = D(1,1)*K2;
    setSparseBlock_update(&Q[i],d[i]/2,d[i]/2, B);
    B = D(0,1)*K2;
    setSparseBlock_update(&Q[i],0,d[i]/2, B);
    B = D(1,0)*K1;
    setSparseBlock_update(&Q[i],d[i]/2,0, B);
    
    // Set kappa1 derivative and trace
    B = (-tau1*pow(kappa1,-2))*G[i];
    B+= tau1*C[i];
    double kappa1_trace_i, kappa2_trace_i;
    //set kappa traces
    (*Qsolver[i]).compute(K1);
    kappa1_trace_i = (*Qsolver[i]).trace(B);
    kappa1_trace[i] = kappa1_trace_i;
    
    B = D(0,0)*B;
    setSparseBlock_update(&dkappa1Q[i],0,0, B);

    B = (-D(1,0)*tau1*pow(kappa1,-2))*G[i];
    B+=D(1,0)*tau1*C[i];
    setSparseBlock_update(&dkappa1Q[i],d[i]/2,0, B);
    
    //set kappa2 derivative and trace
    B = (-tau2*pow(kappa2,-2))*G[i];
    B+= tau2*C[i];
    
    (*Qsolver[i]).compute(K2);
    kappa2_trace_i = (*Qsolver[i]).trace(B);
    kappa2_trace[i] = kappa2_trace_i;
    
    B = D(1,1)*B;
    setSparseBlock_update(&dkappa2Q[i],d[i]/2,d[i]/2, B);
    
    B = (-D(0,1)*tau2*pow(kappa2,-2))*G[i];
    B+=D(0,1)*tau2*C[i];
    setSparseBlock_update(&dkappa2Q[i],0,d[i]/2, B);
    
    //set tau1 derivative
    B = (D(0,0)/kappa1)*G[i];
    B+=D(0,0)*kappa1*C[i];
    setSparseBlock_update(&dtau1Q[i],0,0, B);
    B=(D(1,0)/kappa1)*G[i];
    B+=D(1,0)*kappa1*C[i];
    setSparseBlock_update(&dtau1Q[i],d[i]/2,0, B);
    
    //set tau2 derivative
    B=(D(1,1)/kappa2)*G[i];
    B+=D(1,1)*kappa2*C[i];
    setSparseBlock_update(&dtau2Q[i],d[i]/2,d[i]/2, B);
    B=(D(0,1)/kappa2)*G[i];
    B+=D(0,1)*kappa2*C[i];
    setSparseBlock_update(&dtau2Q[i],0,d[i]/2, B);

    //set rho derivative
    B=Drho(0,0)*K1;
    setSparseBlock_update(&drhoQ[i],0,0, B);
    B=Drho(1,1)*K2;
    setSparseBlock_update(&drhoQ[i],d[i]/2,d[i]/2, B);
    B=Drho(0,1)*K2;
    setSparseBlock_update(&drhoQ[i],0,d[i]/2, B);
    B=Drho(1,0)*K1;
    setSparseBlock_update(&drhoQ[i],d[i]/2,0, B);
    
    //set theta derivative
    if(estimate_theta == 1){
      B=(Dtheta(0,0)*tau1/kappa1)*G[i];
      B+=Dtheta(0,0)*tau1*kappa1*C[i];
      setSparseBlock_update(&dthetaQ[i],0,0, B);
      B=(Dtheta(1,1)*tau2/kappa2)*G[i];
      B+=Dtheta(1,1)*tau2*kappa2*C[i];
      setSparseBlock_update(&dthetaQ[i],d[i]/2,d[i]/2, B);
      B=(Dtheta(0,1)*tau2/kappa2)*G[i];
      B+=Dtheta(0,1)*tau2*kappa2*C[i];
      setSparseBlock_update(&dthetaQ[i],0,d[i]/2, B);
      B=(Dtheta(1,0)*tau1/kappa1)*G[i];
      B+=Dtheta(1,0)*tau1*kappa1*C[i];
      setSparseBlock_update(&dthetaQ[i],d[i]/2,0, B);
    }
    
    //set kappa1 second derivative
    B=(2*D(0,0)*tau1*pow(kappa1,-3))*G[i];
    setSparseBlock_update(&d2kappa1Q[i],0,0, B);
    B=(2*D(1,0)*tau1*pow(kappa1,-3))*G[i];
    setSparseBlock_update(&d2kappa1Q[i],d[i]/2,0, B);
    
    //set kappa2 second derivative
    B=(2*D(1,1)*tau2*pow(kappa2,-3))*G[i];
    setSparseBlock_update(&d2kappa2Q[i],d[i]/2,d[i]/2, B);
    B=(2*D(0,1)*tau2*pow(kappa2,-3))*G[i];
    setSparseBlock_update(&d2kappa2Q[i],0,d[i]/2, B);
    
    //set rho second derivative
    B=Drho2(1,1)*K2;
    setSparseBlock_update(&d2rhoQ[i],d[i]/2,d[i]/2, B);
    B=Drho2(0,1)*K2;
    setSparseBlock_update(&d2rhoQ[i],0,d[i]/2, B);
    
    //set theta second derivative
    if(estimate_theta == 1){
      B=(Dtheta2(0,0)*tau1/kappa1)*G[i];
      B+=Dtheta2(0,0)*tau1*kappa1*C[i];
      setSparseBlock_update(&d2thetaQ[i],0,0, B);
      B=(Dtheta2(1,1)*tau2/kappa2)*G[i];
      B+=Dtheta2(1,1)*tau2*kappa2*C[i];
      setSparseBlock_update(&d2thetaQ[i],d[i]/2,d[i]/2, B);
      B=(Dtheta2(0,1)*tau2/kappa2)*G[i];
      B+=Dtheta2(0,1)*tau2*kappa2*C[i];
      setSparseBlock_update(&d2thetaQ[i],0,d[i]/2, B);
      B=(Dtheta2(1,0)*tau1/kappa1)*G[i];
      B+=Dtheta2(1,0)*tau1*kappa1*C[i];
      setSparseBlock_update(&d2thetaQ[i],d[i]/2,0, B);
    }

    // Set other traces: 
    // trace1 = d log|Q| = trace(dQ*Q^-1)
    // trace2 = d2 log|Q| = trace(d2Q*Q^-1 + dQ*Q^-1*dQ*Q^-1)
    // d log|Q| = d log|DK| = d log|D| +dlog|K|
    // K = diag(tau1*K1,tau2*K2), dK = diag(K1,0)
    MatrixXd tmp = Drho*D.inverse();
    rho_trace[i] = (d[i]/2)*tmp.trace();
    MatrixXd tmp2 = -tmp*tmp + Drho2*D.inverse();
    rho_trace2[i] = (d[i]/2)*tmp2.trace();
    
    if(estimate_theta == 1){
      tmp = Dtheta*D.inverse();
      theta_trace[i] = (d[i]/2)*tmp.trace();
      tmp2 = -tmp*tmp + Dtheta2*D.inverse();
      theta_trace2[i] = (d[i]/2)*tmp2.trace();
    }
    
    tau1_trace[i] = 0.5*d[i]/tau1;
    tau2_trace[i] = 0.5*d[i]/tau2;
    tau1_trace2[i] = -0.5*d[i]/pow(tau1,2);
    tau2_trace2[i] = -0.5*d[i]/pow(tau2,2);
    
    //Numerical approx of trace2 for kappa1 and kappa2
    SparseMatrix<double,0,int> Qeps, dQeps;
    double eps = 0.0001;
    double kappa_eps = kappa1 + eps;
    double trje;
    
    Qeps = (tau1/kappa_eps)*G[i];
    Qeps+=tau1*kappa_eps*C[i];
    (*Qepssolver[i]).compute(Qeps);

    dQeps=(-tau1*pow(kappa_eps,-2))*G[i];
    dQeps+=tau1*C[i];
    
    trje = (*Qepssolver[i]).trace(dQeps);
    kappa1_trace2[i] = (trje - kappa1_trace_i)/eps;
    
    kappa_eps = kappa2 + eps;
    
    Qeps=(tau2/kappa_eps)*G[i];
    Qeps+=tau2*kappa_eps*C[i];
    (*Qepssolver[i]).compute(Qeps);
    
    dQeps=(-tau2*pow(kappa_eps,-2))*G[i];
    dQeps+=tau2*C[i];
    trje = (*Qepssolver[i]).trace(dQeps);
    kappa2_trace2[i] = (trje - kappa2_trace_i)/eps;
    
    
  }
}


Rcpp::List MaternOperator2D::output_list()
{
  out_list["tau1"] = tau1Vec[tau1Vec.size() - 1];
  out_list["kappa1"] = kappa1Vec[tau1Vec.size() - 1];
  out_list["tau1Vec"] = tau1Vec;
  out_list["kappa1Vec"] = kappa1Vec;
  out_list["tau2"] = tau2Vec[tau1Vec.size() - 1];
  out_list["kappa2"] = kappa2Vec[tau1Vec.size() - 1];
  out_list["tau2Vec"] = tau2Vec;
  out_list["kappa2Vec"] = kappa2Vec;
  out_list["rho"] = rhoVec[rhoVec.size() - 1];
  out_list["rhoVec"] = rhoVec;
  out_list["theta"] = thetaVec[thetaVec.size() - 1];
  out_list["thetaVec"] = thetaVec;
  
  out_list["nIter"] = tau1Vec.size();
  out_list["use.chol"] = use_chol;
  out_list["Cov_theta"]   = Cov_theta;
  return(out_list);
}


void MaternOperator2D::gradient_init(int nsim, int nrep)
{
  dtau1 =  0;
  ddtau1 = 0; 
  dkappa1 = 0;
  ddkappa1 = 0;
  dtau2 =  0;
  ddtau2 = 0;
  dkappa2 = 0;
  ddkappa2 = 0;
  drho =  0;
  ddrho = 0;
  dtheta = 0;
  ddtheta = 0;
}


void MaternOperator2D::gradient_add( const Eigen::VectorXd & X,
                                     const Eigen::VectorXd & iV,
                                     const Eigen::VectorXd & mean_KX,
                                     const int ii,
                                     const double weight)
{
  int i = ii;
  if(nop == 1)
    i = 0;
  
  this->set_matrix(i);
  dtau1    += weight * tau1_trace[i];
  ddtau1   += weight * tau1_trace2[i];
  dkappa1  += weight * kappa1_trace[i];
  ddkappa1 += weight * kappa1_trace2[i];
  
  dtau2    += weight * tau2_trace[i];
  ddtau2   += weight * tau2_trace2[i];
  dkappa2  += weight * kappa2_trace[i];
  ddkappa2 += weight * kappa2_trace2[i];
  
  drho    += weight * rho_trace[i];
  ddrho   += weight * rho_trace2[i];
  dtheta  += weight * theta_trace[i];
  ddtheta += weight * theta_trace2[i];
  
  Eigen::VectorXd KX = Q[i]*X;
  //compute gradients wrt tau
  Eigen::VectorXd dKX = dtau1Q[i] * X;
  dtau1 -=  weight * dKX.dot(iV.asDiagonal() * (KX - mean_KX));
  ddtau1 -= weight * (dKX.dot(iV.asDiagonal() * dKX));
  
  dKX = dtau2Q[i] * X;
  dtau2 -=  weight * dKX.dot(iV.asDiagonal() * (KX - mean_KX));
  ddtau2 -= weight * (dKX.dot(iV.asDiagonal() * dKX));
  
  //compute gradients wrt kappa
  Eigen::VectorXd d2KX;
  dKX      = dkappa1Q[i] * X;
  d2KX     = d2kappa1Q[i] * X;
  dkappa1  -= weight * dKX.dot(iV.asDiagonal() * (KX - mean_KX));
  ddkappa1 -= weight * dKX.dot(iV.asDiagonal()*dKX);
  ddkappa1 -= weight * d2KX.dot(iV.asDiagonal()*(KX - mean_KX));
  
  dKX      = dkappa2Q[i] * X;
  d2KX     = d2kappa2Q[i] * X;
  dkappa2  -= weight * dKX.dot(iV.asDiagonal() * (KX - mean_KX));
  ddkappa2 -= weight * dKX.dot(iV.asDiagonal()*dKX);
  ddkappa2 -= weight * d2KX.dot(iV.asDiagonal()*(KX - mean_KX));
  
  //compute gradients wrt rho and theta
  dKX      = drhoQ[i] * X;
  d2KX     = d2rhoQ[i] * X;
  drho  -= weight * dKX.dot(iV.asDiagonal() * (KX - mean_KX));
  ddrho -= weight * dKX.dot(iV.asDiagonal()*dKX);
  ddrho -= weight * d2KX.dot(iV.asDiagonal()*(KX - mean_KX));
  
  if(estimate_theta == 1){
    dKX      = dthetaQ[i] * X;
    d2KX     = d2thetaQ[i] * X;
    dtheta  -= weight * dKX.dot(iV.asDiagonal() * (KX - mean_KX));
    ddtheta -= weight * dKX.dot(iV.asDiagonal()*dKX);
    ddtheta -= weight * d2KX.dot(iV.asDiagonal()*(KX - mean_KX));
  }
}

void MaternOperator2D::gradient( const Eigen::VectorXd & X, const Eigen::VectorXd & iV)
{
  throw(" MaternOperator::gradient depricated \n");
}

void MaternOperator2D::print_parameters(){
  Rcpp::Rcout << "tau1 = " << tau1 << "\n";
  Rcpp::Rcout << "tau2 = " << tau2 << "\n";
  Rcpp::Rcout << "kappa1 = " << kappa1 << "\n";
  Rcpp::Rcout << "kappa2 = " << kappa2 << "\n";
  Rcpp::Rcout << "rho = " << rho << "\n";
  Rcpp::Rcout << "theta = " << theta << "\n";
}


void MaternOperator2D::step_theta(const double stepsize,
                                  const double learning_rate,
                                  const double polyak_rate,
                                  const int burnin)
{
  
    
  //step tau1
  dtau1  /= ddtau1;
  dtau1_old = learning_rate * dtau1_old + dtau1;
  double step = stepsize * dtau1_old;
  double par_temp = -1.;
  while(par_temp < 0)
  {
    step *= 0.5;
    par_temp = tau1 - step;
  }
  tau1 = par_temp;
  
  //step tau2
  dtau2  /= ddtau2;
  dtau2_old = learning_rate * dtau2_old + dtau2;
  step = stepsize * dtau2_old;
  par_temp = -1.;
  while(par_temp < 0)
  {
    step *= 0.5;
    par_temp = tau2 - step;
  }
  tau2 = par_temp;
  
  
  //step kappa1
  dkappa1  /= ddkappa1;
  dkappa1_old = learning_rate * dkappa1_old + dkappa1;
  step   = stepsize * dkappa1_old;
  par_temp = -1.;
  while(par_temp < 0)
  {
    step *= 0.5;
    par_temp = kappa1 - step;
  }
  kappa1 = par_temp;
  
  //step kappa2
  dkappa2  /= ddkappa2;
  dkappa2_old = learning_rate * dkappa2_old + dkappa2;
  step   = stepsize * dkappa2_old;
  par_temp = -1.;
  while(par_temp < 0)
  {
    step *= 0.5;
    par_temp = kappa2 - step;
  }
  kappa2 = par_temp;
  
  //step rho
  drho  /= ddrho;
  drho_old = learning_rate * drho_old + drho;
  step   = stepsize * drho_old;
  rho = rho - step;
  
  if(estimate_theta){
  //step theta
  dtheta  /= ddtheta;
  dtheta_old = learning_rate * dtheta_old + dtheta;
  step   = stepsize * dtheta_old;
  theta = theta - step;
  }
  if(counter == 0 || polyak_rate == -1){
    tau1Vec[counter] = tau1;
    tau2Vec[counter] = tau2;
    kappa1Vec[counter] = kappa1;
    kappa2Vec[counter] = kappa2;
    rhoVec[counter] = rho;
    thetaVec[counter] = theta;
  }else{
    tau1Vec[counter] = polyak_rate * tau1 + (1 - polyak_rate) * tau1Vec[counter-1];
    tau2Vec[counter] = polyak_rate * tau2 + (1 - polyak_rate) * tau2Vec[counter-1];
    kappa1Vec[counter] = polyak_rate * kappa1 + (1 - polyak_rate) * kappa1Vec[counter-1];
    kappa2Vec[counter] = polyak_rate * kappa2 + (1 - polyak_rate) * kappa2Vec[counter-1];
    rhoVec[counter] = polyak_rate * rho + (1 - polyak_rate) * rhoVec[counter-1];
    thetaVec[counter] = polyak_rate * theta + (1 - polyak_rate) * thetaVec[counter-1];
  }

  counter++;
  clear_gradient();
  ddtau1   = 0;
  ddtau2   = 0;
  ddkappa1 = 0;
  ddkappa2 = 0;
  ddrho   = 0;
  ddtheta = 0;
  for(int i=0;i<nop;i++){
    matrix_set[i] = 0;
  }
}

double MaternOperator2D::trace_variance( const Eigen::SparseMatrix<double,0,int> & A, int i)
{
  //return(A.rows() * tau/ h_average);
  return(-1);
}

Eigen::VectorXd  MaternOperator2D::get_gradient()
{
  Eigen::VectorXd g(npars);
    g[0] = dtau1;
    g[1] = dtau2;
    g[2] = dkappa1;
    g[3] = dkappa2;
    g[4] = drho;
    g[5] = dtheta;
  return(g);
}
void  MaternOperator2D::clear_gradient()
{
  dtau1   = 0;
  dtau2   = 0;
  dkappa1 = 0;
  dkappa2 = 0;
  drho = 0;
  dtheta = 0;
};


Eigen::MatrixXd MaternOperator2D::d2Given( const Eigen::VectorXd & X,
                                         const Eigen::VectorXd & iV,
                                         const Eigen::VectorXd & mean_KX,
                                         int ii,
                                         const double weight)
{
  Rcpp::Rcout << "Warning: d2Given not completely implemented for MaternOperator2D\n";
  int i = ii;
  if(nop == 1)
    i = 0;
  this->set_matrix(i);
  Eigen::VectorXd vtmp = Q[i] * X;
  
  
  Eigen::VectorXd KX = Q[i]*X;
  Eigen::VectorXd dKX;
  Eigen::VectorXd d2KX;
  Eigen::MatrixXd d2 = Eigen::MatrixXd::Zero(6, 6);
  
  //compute gradients wrt tau
  //d2tau is zero so no addition from that term
  dKX = dtau1Q[i] * X;
  d2(0, 0)  =- weight * tau1_trace2[i];
  d2(0, 0) -=- weight * (dKX.dot(iV.asDiagonal() * dKX));
  
  dKX = dtau2Q[i] * X;
  d2(1, 1)  =- weight * tau2_trace2[i];
  d2(1, 1) -=- weight * (dKX.dot(iV.asDiagonal() * dKX));
  
  //compute gradients wrt kappa
  dKX      = dkappa1Q[i] * X;
  d2KX     = d2kappa1Q[i] * X;
  d2(2, 2)  =- weight * kappa1_trace2[i];
  d2(2, 2) -=- weight * dKX.dot(iV.asDiagonal()*dKX);
  d2(2, 2) -=- weight * d2KX.dot(iV.asDiagonal()*(KX - mean_KX));
  
  dKX      = dkappa2Q[i] * X;
  d2KX     = d2kappa2Q[i] * X;
  d2(3, 3)  =- weight * kappa2_trace2[i];
  d2(3, 3) -=- weight * dKX.dot(iV.asDiagonal()*dKX);
  d2(3, 3) -=- weight * d2KX.dot(iV.asDiagonal()*(KX - mean_KX));
  
  //compute gradients wrt rho and tau
  dKX      = drhoQ[i] * X;
  d2KX     = d2rhoQ[i] * X;
  d2(4, 4)  =- weight * rho_trace2[i];
  d2(4, 4) -=- weight * dKX.dot(iV.asDiagonal()*dKX);
  d2(4, 4) -=- weight * d2KX.dot(iV.asDiagonal()*(KX - mean_KX));
  
  dKX      = dthetaQ[i] * X;
  d2KX     = d2thetaQ[i] * X;
  d2(5, 5)  =- weight * theta_trace2[i];
  d2(5, 5) -=- weight * dKX.dot(iV.asDiagonal()*dKX);
  d2(5, 5) -=- weight * d2KX.dot(iV.asDiagonal()*(KX - mean_KX));
  
  //compute cross terms
  //All cross term tau1-tau2, kappa1-kappa2, tau1-kappa2 are zero
  
  //kappa-tau
  
  //kappa-rho
  
  //kappa-theta
  
  //tau-rho
  
  //tau-theta
  
  //rho-theta

  
  d2(0, 1) -=- weight * dKX.dot(iV.asDiagonal() * (2 * KX - mean_KX))/tau;
  d2(1, 0) = d2(0, 1);
  
  return(d2);
  
}
