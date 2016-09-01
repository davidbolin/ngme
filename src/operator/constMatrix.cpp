#include "operatorMatrix.h"
#include "error_check.h"

using namespace std;

void constMatrix::initFromList(Rcpp::List const & init_list)
{
	npars  = 1;
 std::vector<std::string> check_names =  {"Q","loc", "h"};
 check_Rcpplist(init_list, check_names, "constMatrix::initFromList");

 tau = 1.;
 if(init_list.containsElementNamed("tau"))
   tau = Rcpp::as<double >( init_list["tau"]);

  Rcpp::List Q_list  = Rcpp::as<Rcpp::List> (init_list["Q"]);
  Rcpp::List loc_list  = Rcpp::as<Rcpp::List> (init_list["loc"]);
  Rcpp::List h_list  = Rcpp::as<Rcpp::List> (init_list["h"]);

  nop = Q_list.size();
  Q = new Eigen::SparseMatrix<double,0,int>[nop];
  d = new int[nop];
  loc = new Eigen::VectorXd[nop];

  for(int i=0;i<nop;i++){
      SEXP tmp = Q_list[i];
      Q[i] =  Rcpp::as<Eigen::SparseMatrix<double,0,int>>(tmp);
      d[i] = Q[i].rows();
      Q[i] *= tau;

      tmp = loc_list[i];
      loc[i]  = Rcpp::as< Eigen::VectorXd >(tmp);
      tmp = h_list[i];
      h[i]  = Rcpp::as< Eigen::VectorXd >(tmp);
      h_average[i] = h[i].sum() / h[i].size();
      m_loc[i] = loc[i].minCoeff();
  }

  int nIter = Rcpp::as<double>(init_list["nIter"]);
  tauVec.resize(nIter+1);
  npars = 0;
  v.setZero(1);
  m.resize(1,1);

  dtau  = 0.;
  ddtau = 0.;
  counter = 0;
  tauVec[counter] = tau;
  counter++;

}

void constMatrix::initFromList(Rcpp::List const & init_List, Rcpp::List const & solver_list)
{
  this->initFromList(init_List);
}

void constMatrix::gradient_init(int nsim, int nrep)
{
  dtau   = 0;
  ddtau  = 0;
}

void constMatrix::gradient_add( const Eigen::VectorXd & X,
								   const Eigen::VectorXd & iV,
								   const Eigen::VectorXd & mean_KX,
								   const int ii)
{
  Eigen::VectorXd vtmp = Q[ii] * X;

  double xtQx =  vtmp.dot( iV.asDiagonal() * vtmp);
  double xtQmean = - vtmp.dot( iV.asDiagonal() * mean_KX);
  dtau +=  (d[ii] - xtQx - xtQmean)/ tau;
  ddtau -=  (d[ii] + xtQx) / pow(tau, 2);
}

void constMatrix::gradient( const Eigen::VectorXd & X, const Eigen::VectorXd & iV)
{
  throw(" constMatrix::gradient depricated \n");
}

void constMatrix::print_parameters(){
  Rcpp::Rcout << "tau = " << tau << "\n";
}

void constMatrix::step_theta(const double stepsize)
{

	dtau  /= ddtau;
  dtau *= stepsize;
	double tau_temp = -1.;
    while(tau_temp < 0)
    {
    	dtau *= 0.5;
        tau_temp = tau - dtau;
    }
  for(int i=0;i<nop;i++){
    Q[i] *= tau_temp/tau;
  }

	tau = tau_temp;
	tauVec[counter] = tau;

	counter++;
	clear_gradient();
	ddtau  = 0;
}

Rcpp::List constMatrix::output_list()
{
  Rcpp::List  List;
  List["tau"] = tau;
  List["tauVec"] = tauVec;
  List["Q"] = Q;
  List["loc"] = loc;
  List["nIter"] = tauVec.size();
  List["h"] = h;
  List["Cov_theta"]   = Cov_theta;
  return(List);
}
