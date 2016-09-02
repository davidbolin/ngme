#include "operatorMatrix.h"
#include "error_check.h"
#include "eigen_add_on.h"
#include "operator_helper.h"

void MaternOperator::initFromList(Rcpp::List const & init_list, Rcpp::List const & solver_list)
{
/*
	npars = 2;
  std::vector<std::string> check_names =  {"C", "G", "kappa", "tau","h"};
  check_Rcpplist(init_list, check_names, "MaternOperator::initFromList");
  std::vector<std::string> check_names2 =  {"use.chol"};
  check_Rcpplist(solver_list, check_names2, "MaternOperator::initFromList");
  G  = Rcpp::as<Eigen::SparseMatrix<double,0,int> >(init_list["G"]);
  C  = Rcpp::as<Eigen::SparseMatrix<double,0,int> >(init_list["C"]);
  Q = G + C;
  n = G.rows();
  d = n;
  kappa = Rcpp::as<double>(init_list["kappa"]);
  tau = Rcpp::as<double>(init_list["tau"]);

  int nIter = Rcpp::as<double>(init_list["nIter"]);
  kappaVec.resize(nIter);
  tauVec.resize(nIter);
  d2tauQ.resize(n,n);

  dkappa = 0;
  dtau = 0;
	h = h;
	h_average = h.sum() / h.size();
  use_chol = Rcpp::as<int>(solver_list["use.chol"]);

  if(use_chol==1){
    Qsolver = new cholesky_solver;
    Qepssolver = new cholesky_solver;
  } else {
    Qsolver = new iterative_solver;
    Qepssolver = new iterative_solver;
  }
  (*Qsolver).initFromList(n,solver_list);
  (*Qsolver).analyze(Q);

  (*Qepssolver).initFromList(n,solver_list);
  (*Qepssolver).analyze(Q);

  this->set_matrices();
*/
}


void MaternOperator::set_matrices()
{
/*
  double c = 0.5; //sqrt(gamma(alpha))/(sqrt(gamma(nu))*(4*pi)^(d/4))
  Q = c*tau*(pow(kappa,-1.5)*G + pow(kappa,0.5)*C);
  dkappaQ = c*tau*(-1.5*pow(kappa,-2.5)*G + 0.5*pow(kappa,-0.5)*C);
  d2kappaQ = c*tau*(1.5*2.5*pow(kappa,-3.5)*G - 0.5*0.5*pow(kappa,-1.5)*C);
  dtauQ = c*(pow(kappa,-1.5)*G + pow(kappa,0.5)*C);

  (*Qsolver).compute(Q);
  tau_trace = (*Qsolver).trace(dtauQ);
  tau_trace2 = -tau_trace/tau;
  kappa_trace = (*Qsolver).trace(dkappaQ);

  double eps = 0.0001;
  double kappa_eps = kappa + eps;
  SparseMatrix<double,0,int> Qeps = c*tau*(pow(kappa_eps,-1.5)*G + pow(kappa_eps,0.5)*C);
  SparseMatrix<double,0,int> dQeps = c*tau*(-1.5*pow(kappa_eps,-2.5)*G + 0.5*pow(kappa_eps,-0.5)*C);
  (*Qepssolver).compute(Qeps);
  double trje = (*Qepssolver).trace(dQeps);
  kappa_trace2 = (trje - kappa_trace)/eps;
  */
}


Rcpp::List MaternOperator::output_list()
{
  Rcpp::List  List;
  /*
  List["type"] = "Matern";
  List["tau"] = tau;
  List["kappa"] = kappa;
  List["tauVec"] = tauVec;
  List["kappaVec"] = kappaVec;
  List["G"] = G;
  List["C"] = C;
  List["nIter"] = tauVec.size();
  List["use.chol"] = use_chol;
  List["Cov_theta"]   = Cov_theta;
  */
  return(List);
}


void MaternOperator::gradient_init(int nsim, int nrep)
{
  dtau =  nsim*nrep*tau_trace;
  ddtau = nsim*nrep*tau_trace2;
  dkappa = nsim*nrep*kappa_trace;
  ddkappa = nsim*nrep*kappa_trace2;
}


void MaternOperator::gradient_add( const Eigen::VectorXd & X,
								   const Eigen::VectorXd & iV,
								   const Eigen::VectorXd & mean_KX,
                   const int i)
{
/*
  Eigen::VectorXd KX = Q * X;

  //compute gradients wrt tau
  Eigen::VectorXd dKX = dtauQ * X;
  Eigen::VectorXd d2KX = d2tauQ * X;
  dtau -=  dKX.dot(iV.asDiagonal() * (KX - mean_KX));
  ddtau -= 0.5*(dKX.dot(iV.asDiagonal() * dKX) + d2KX.dot(iV.asDiagonal() * (KX - mean_KX)));

  //compute gradients wrt kappa
  dKX = dkappaQ * X;
  d2KX = d2kappaQ * X;
  dkappa -= dKX.dot(iV.asDiagonal() * (KX - mean_KX));
  ddkappa -= 0.5*(dKX.dot(iV.asDiagonal() * dKX) + d2KX.dot(iV.asDiagonal() *(KX - mean_KX)));
*/
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
