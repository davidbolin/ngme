#include "operatorMatrix.h"
#include "error_check.h"
#include "eigen_add_on.h"
#include "operator_helper.h"

MaternOperator::~MaternOperator()
{
  delete Q;
  delete d2tauQ;
  delete dtauQ;
  delete dkappaQ;
  delete d2kappaQ;
  delete Qsolver;
  delete Qepssolver;
}

void MaternOperator::initFromList(Rcpp::List const & init_list, Rcpp::List const & solver_list)
{
	npars = 2;
  std::vector<std::string> check_names =  {"C", "G", "kappa", "tau","h"};
  check_Rcpplist(init_list, check_names, "MaternOperator::initFromList");
  std::vector<std::string> check_names2 =  {"use.chol"};
  check_Rcpplist(solver_list, check_names2, "MaternOperator::initFromList");

  kappa = Rcpp::as<double>(init_list["kappa"]);
  tau = Rcpp::as<double>(init_list["tau"]);
  int nIter = Rcpp::as<double>(init_list["nIter"]);
  kappaVec.resize(nIter);
  tauVec.resize(nIter);
  dkappa = 0;
  dtau = 0;
  use_chol = Rcpp::as<int>(solver_list["use.chol"]);


  Rcpp::List G_list  = Rcpp::as<Rcpp::List> (init_list["G"]);
  Rcpp::List C_list  = Rcpp::as<Rcpp::List> (init_list["C"]);
  Rcpp::List h_list  = Rcpp::as<Rcpp::List> (init_list["h"]);
  nop = G_list.size();

  d.resize(nop);
  loc.resize(nop);
  h.resize(nop);
  h_average.resize(nop);

  Q = new Eigen::SparseMatrix<double,0,int>[nop];
  G = new Eigen::SparseMatrix<double,0,int>[nop];
  C = new Eigen::SparseMatrix<double,0,int>[nop];
  d2tauQ = new Eigen::SparseMatrix<double,0,int>[nop];
  dtauQ = new Eigen::SparseMatrix<double,0,int>[nop];
  dkappaQ = new Eigen::SparseMatrix<double,0,int>[nop];
  d2kappaQ = new Eigen::SparseMatrix<double,0,int>[nop];

  Qsolver = new solver*[nop];
  Qepssolver = new solver*[nop];

  for(int i=0;i<nop;i++){

    if(use_chol==1){
      Qsolver[i] = new cholesky_solver;
      Qepssolver[i] = new cholesky_solver;
    } else {
      Qsolver[i] = new iterative_solver;
      Qepssolver[i] = new iterative_solver;
    }

    G[i] =  Rcpp::as<Eigen::SparseMatrix<double,0,int>>(G_list[i]);
    C[i] =  Rcpp::as<Eigen::SparseMatrix<double,0,int>>(C_list[i]);

    Q[i] = G[i] + C[i];
    d[i] = Q[i].rows();
    //d2tauQ.resize(n[i],n[i]);

    h[i]  = Rcpp::as< Eigen::VectorXd >(h_list[i]);

    h_average[i] = h[i].sum() / h[i].size();
    (*Qsolver[i]).initFromList(d[i],solver_list);
    (*Qsolver[i]).analyze(Q[i]);

    (*Qepssolver[i]).initFromList(d[i],solver_list);
    (*Qepssolver[i]).analyze(Q[i]);
  }


  this->set_matrices();
}


void MaternOperator::set_matrices()
{
  SparseMatrix<double,0,int> Qeps, dQeps;
  double eps = 0.0001;
  double kappa_eps = kappa + eps;
  double trje, kappa_trace_i;
  double c = 0.5; //sqrt(gamma(alpha))/(sqrt(gamma(nu))*(4*pi)^(d/4))
  tau_trace = 0;
  tau_trace2 = 0;
  kappa_trace = 0;
  kappa_trace2 = 0;

  for(int i=0;i<nop;i++){
    Q[i] = c*tau*(pow(kappa,-1.5)*G[i] + pow(kappa,0.5)*C[i]);
    dkappaQ[i] = c*tau*(-1.5*pow(kappa,-2.5)*G[i] + 0.5*pow(kappa,-0.5)*C[i]);
    d2kappaQ[i] = c*tau*(1.5*2.5*pow(kappa,-3.5)*G[i] - 0.5*0.5*pow(kappa,-1.5)*C[i]);
    dtauQ[i] = c*(pow(kappa,-1.5)*G[i] + pow(kappa,0.5)*C[i]);
    (*Qsolver[i]).compute(Q[i]);
    tau_trace += (*Qsolver[i]).trace(dtauQ[i]);
    tau_trace2 += -tau_trace/tau;
    kappa_trace_i = (*Qsolver[i]).trace(dkappaQ[i]);
    kappa_trace += kappa_trace_i;
    Qeps = c*tau*(pow(kappa_eps,-1.5)*G[i] + pow(kappa_eps,0.5)*C[i]);
    dQeps = c*tau*(-1.5*pow(kappa_eps,-2.5)*G[i] + 0.5*pow(kappa_eps,-0.5)*C[i]);
    (*Qepssolver[i]).compute(Qeps);
    trje = (*Qepssolver[i]).trace(dQeps);
    kappa_trace2 += (trje - kappa_trace_i)/eps;
  }
}


Rcpp::List MaternOperator::output_list()
{
  Rcpp::List  List;

  List["type"] = "Matern";
  List["tau"] = tau;
  List["kappa"] = kappa;
  List["tauVec"] = tauVec;
  List["kappaVec"] = kappaVec;
  //List["G"] = G;
  //List["C"] = C;
  List["nIter"] = tauVec.size();
  List["use.chol"] = use_chol;
  List["Cov_theta"]   = Cov_theta;
  return(List);
}


void MaternOperator::gradient_init(int nsim, int nrep)
{
  dtau =  nsim*tau_trace;
  ddtau = nsim*tau_trace2;
  dkappa = nsim*kappa_trace;
  ddkappa = nsim*kappa_trace2;
}


void MaternOperator::gradient_add( const Eigen::VectorXd & X,
								   const Eigen::VectorXd & iV,
								   const Eigen::VectorXd & mean_KX,
                   const int i)
{
  Eigen::VectorXd KX = Q[i] * X;

  //compute gradients wrt tau
  Eigen::VectorXd dKX = dtauQ[i] * X;
  Eigen::VectorXd d2KX;
  //Eigen::VectorXd d2KX = d2tauQ * X;
  dtau -=  dKX.dot(iV.asDiagonal() * (KX - mean_KX));
  ddtau -= 0.5*(dKX.dot(iV.asDiagonal() * dKX));// + d2KX.dot(iV.asDiagonal() * (KX - mean_KX)));

  //compute gradients wrt kappa
  dKX = dkappaQ[i] * X;
  d2KX = d2kappaQ[i] * X;
  dkappa -= dKX.dot(iV.asDiagonal() * (KX - mean_KX));
  ddkappa -= 0.5*(dKX.dot(iV.asDiagonal() * dKX) + d2KX.dot(iV.asDiagonal() *(KX - mean_KX)));
}

void MaternOperator::gradient( const Eigen::VectorXd & X, const Eigen::VectorXd & iV)
{
 throw(" MaternOperator::gradient depricated \n");
}

void MaternOperator::print_parameters(){
  Rcpp::Rcout << "tau = " << tau << "\n";
  Rcpp::Rcout << "kappa = " << kappa << "\n";
}


void MaternOperator::step_theta(const double stepsize)
{

	dtau  /= ddtau;
  dtau *= stepsize;
	double tau_temp = -1.;
    while(tau_temp < 0)
    {
    	dtau *= 0.5;
        tau_temp = tau - dtau;
    }
	tau = tau_temp;

	dkappa  /= ddkappa;
	dkappa *= stepsize;
	double kappa_temp = -1.;
	while(kappa_temp < 0)
	{
	  dkappa *= 0.5;
	  kappa_temp = kappa - dkappa;
	}
	kappa = kappa_temp;
	tauVec[counter] = tau;
	kappaVec[counter] = kappa;
	counter++;
	clear_gradient();
	ddtau   = 0;
	ddkappa = 0;
  	this->set_matrices();
}

double MaternOperator::trace_variance( const Eigen::SparseMatrix<double,0,int> & A, int i)
{
	//return(A.rows() * tau/ h_average);
  return(-1);
}

Eigen::VectorXd  MaternOperator::get_gradient()
{
	Eigen::VectorXd g(npars);
	g[0] = dtau;
	g[1] = dkappa;
	return(g);
}
void  MaternOperator::clear_gradient()
{
	dkappa = 0;
	dtau   = 0;

};
