#ifndef __CONSTMATRIX__H__
#define __CONSTMATRIX__H__


#include "operatorMatrix.h"
#include "error_check.h"

using namespace std;



constMatrix::~constMatrix()
{
  Rcpp::Rcout << "in constMatrix::~constMatrix()\n";
	//for(int i = 0; i < nop; i++)
	 // Q[i].~sparseMatrix<double,0,int>();
	delete Q;
}

void constMatrix::initFromList(Rcpp::List const & init_list)
{
  npars  = 1;
 std::vector<std::string> check_names =  {"Q","loc", "h"};
 check_Rcpplist(init_list, check_names, "constMatrix::initFromList");
	dtau_old = 0;
 tau = 1.;
 if(init_list.containsElementNamed("tau"))
   tau = Rcpp::as<double >( init_list["tau"]);
  Rcpp::List Q_list  = Rcpp::as<Rcpp::List> (init_list["Q"]);
  Rcpp::List loc_list  = Rcpp::as<Rcpp::List> (init_list["loc"]);
  Rcpp::List h_list  = Rcpp::as<Rcpp::List> (init_list["h"]);
  nop = Q_list.size();
  Q = new Eigen::SparseMatrix<double,0,int>[nop];
  d.resize(nop);
  loc.resize(nop);
  h.resize(nop);
  h_average.resize(nop);
  m_loc.resize(nop);
  for(int i=0;i<nop;i++){
      //SEXP tmp = Q_list[i];
      Q[i] =  Rcpp::as<Eigen::SparseMatrix<double,0,int>>(Q_list[i]);
      d[i] = Q[i].rows();
      Q[i] *= tau;


      loc[i]  = Rcpp::as< Eigen::VectorXd >( loc_list[i]);
      h[i]  = Rcpp::as< Eigen::VectorXd >(h_list[i]);
      h_average[i] = h[i].sum() / h[i].size();
      m_loc[i] = loc[i].minCoeff();
  }

  int nIter = Rcpp::as<double>(init_list["nIter"]);
  tauVec.resize(nIter+1);
  v.setZero(1);
  m.resize(1,1);

  dtau  = 0.;
  ddtau = 0.;
  counter = 0;
  tauVec[counter] = tau;
  counter++;
  term1 = 0;
  term2 = 0;
  term3 = 0;

}

void constMatrix::initFromList(Rcpp::List const & init_List, Rcpp::List const & solver_list)
{
  this->initFromList(init_List);
}

void constMatrix::gradient_init(int nsim, int nrep)
{
  dtau   = 0;
  ddtau  = 0;
  term1 = 0;
  term2 = 0;
  term3 = 0;
}
Eigen::MatrixXd constMatrix::d2Given( const Eigen::VectorXd & X,
                   const Eigen::VectorXd & iV,
                   const Eigen::VectorXd & mean_KX,
                  int ii,
                  const double weight)
{
  if(nop == 1)
    ii = 0;
  Eigen::VectorXd vtmp = Q[ii] * X;

  double xtQx =  vtmp.dot( iV.asDiagonal() * vtmp);
  double xtQmean = - vtmp.dot( iV.asDiagonal() * mean_KX);
  Eigen::MatrixXd d2 = Eigen::MatrixXd::Zero(1,1);
  d2(0, 0) = weight * (d[ii] + xtQx ) / pow(tau, 2);

  return(d2);

}
void constMatrix::gradient_add( const Eigen::VectorXd & X,
								   const Eigen::VectorXd & iV,
								   const Eigen::VectorXd & mean_KX,
								  int ii,
								  const double weight)
{
	if(nop == 1)
		ii = 0;
  Eigen::VectorXd vtmp = Q[ii] * X;

  double xtQx =  vtmp.dot( iV.asDiagonal() * vtmp);
  double xtQmean = - vtmp.dot( iV.asDiagonal() * mean_KX);
  dtau  -=  weight *  (d[ii] - xtQx - xtQmean)/ tau;
  ddtau += weight * (d[ii] + xtQx) / pow(tau, 2);
  term1 += weight * xtQx/pow(tau,2);
  term2 += weight * xtQmean/tau;
  term3 -= weight * d[ii];
}

void constMatrix::gradient( const Eigen::VectorXd & X, const Eigen::VectorXd & iV)
{
  throw(" constMatrix::gradient depricated \n");
}

void constMatrix::print_parameters(){
  Rcpp::Rcout << "tau = " << tau ;
}

void constMatrix::step_theta(const double stepsize,
							 const double learning_rate,
							 const double polyak_rate,
							 const int burnin)
{
	dtau  /= ddtau;
  	dtau_old = learning_rate * dtau_old + dtau;
  	double step = stepsize * dtau_old;
	double tau_temp = -1.;
    while(tau_temp < 0)
    {
    	step *= 0.5;
        tau_temp = tau - step;
    }

    if(burnin == 1){
      tau_temp = -term2 + pow(term2*term2 - 4*term1*term3,0.5);
      tau_temp /= 2*term1;
      if(tau < 0)
        tau = 1;

      term1 = 0;
      term2 = 0;
      term3 = 0;
    }
  for(int i=0;i<nop;i++){
    Q[i] *= tau_temp/tau;
  }

	tau = tau_temp;

	if(counter == 0 || polyak_rate == -1)
		tauVec[counter] = tau;
	else
		tauVec[counter] = polyak_rate * tau + (1 - polyak_rate) * tauVec[counter-1];

	counter++;
	clear_gradient();
	ddtau  = 0;
}

Rcpp::List constMatrix::output_list()
{
  Rcpp::List  List;
  List["tau"] = tauVec(tauVec.size() -1 );
  List["tauVec"] = tauVec;
  //List["Q"] = Q;
  List["loc"] = loc;
  List["nIter"] = tauVec.size();
  List["h"] = h;
  List["Cov_theta"]   = Cov_theta;
  return(List);
}

#endif
