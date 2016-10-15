#include "MixedEffect.h"
#include "error_check.h"
#include  <cmath>
#include <unsupported/Eigen/KroneckerProduct>

#include <Eigen/LU>


Eigen::VectorXd sampleNormalCan(const Eigen::VectorXd & b,const Eigen::MatrixXd & Q)
{
  return b;
}



NormalMixedEffect::NormalMixedEffect(){
  counter = 0;
  noise = "Normal";
  npars = 0;
  store_param  = 0;
  //dlog_sigma2  = 0;
  //ddlog_sigma2 = 0;
}


void NormalMixedEffect::printIter()
{
	if(Bf.size() > 0)
		Rcpp::Rcout << "beta_f = " << beta_fixed.transpose() << "\n";


	if(Br.size() > 0){
		Rcpp::Rcout << "beta_r = " << beta_random.transpose() << "\n";
	}
}

void NormalMixedEffect::setupStoreTracj(const int Niter)
{

	if(Bf.size() > 0)
		betaf_vec.resize(Niter, Bf[0].cols());

	if(Br.size() > 0){
		betar_vec.resize(Niter, Br[0].cols());
		Sigma_vec.resize(Niter, pow(Br[0].cols(), 2));
	}
	vec_counter = 0;
	store_param = 1;
}


Rcpp::List NormalMixedEffect::toList()
{
  Rcpp::List out;
  out["B_fixed"]          = Bf;
  out["B_random"]          = Br;
  out["beta_random"] = beta_random;
  out["beta_fixed"]  = beta_fixed;
  out["Sigma"]       = Sigma;
  out["U"]           = U;
  out["noise"]       = noise;
  out["Sigma_epsilon"]       = Sigma_epsilon;
  out["Cov_theta"]   = Cov_theta;
  if(store_param){
  if(Bf.size() > 0)
  	if(Bf.size() > 0){
		 out["betaf_vec"] = betaf_vec;
		 out["beta_fixed"] = betaf_vec.row(betaf_vec.rows() - 1);
	}

   if(Br.size() > 0){
		out["betar_vec"]   = betar_vec;
		out["beta_random"] = betar_vec.row(betar_vec.rows() - 1);
		out["Sigma_vec"]   = Sigma_vec;
		Eigen::VectorXd temp = Sigma_vec.row(betar_vec.rows() - 1);
		out["Sigma"]       = veci(temp, Sigma.rows(), Sigma.cols());
	}
  }
  return(out);
}

void NormalMixedEffect::initFromList(Rcpp::List const &init_list)
{
  int count =0;
  if(init_list.containsElementNamed("B_fixed"))
  {
    Rcpp::List Bf_list = init_list["B_fixed"];
    Bf.resize(Bf_list.length());
    for( Rcpp::List::iterator it = Bf_list.begin(); it != Bf_list.end(); ++it ) {
      Bf[count++] = Rcpp::as < Eigen::MatrixXd >( it[0]);
    }
    grad_beta_f.setZero(Bf[0].cols());

    if(init_list.containsElementNamed("beta_fixed"))
      beta_fixed = Rcpp::as < Eigen::VectorXd >( init_list["beta_fixed"]);
    else
      beta_fixed.setZero(Bf[0].cols());
     npars += Bf[0].cols();
     H_beta_fixed.setZero(Bf[0].cols(), Bf[0].cols());
	
	dbeta_f_old.setZero(Bf[0].cols());

  }else{ Bf.resize(0);}
  count = 0;
  if(init_list.containsElementNamed("B_random"))
  {
    Rcpp::List Br_list = init_list["B_random"];
    Br.resize(Br_list.length());
    for( Rcpp::List::iterator it = Br_list.begin(); it != Br_list.end(); ++it ) {
      Br[count++] = Rcpp::as < Eigen::MatrixXd >( it[0]);
    }
    grad_beta_r.setZero(Br[0].cols());
    grad_beta_r2.setZero(Br[0].cols());
    if(init_list.containsElementNamed("beta_random"))
      beta_random = Rcpp::as < Eigen::VectorXd >( init_list["beta_random"]);
    else
      beta_random.setZero(Br[0].cols());

	
	dbeta_r_old.setZero(Br[0].cols());
    npars += Br[0].cols();
  }else{ Br.resize(0);}

  if(Br.size() > 0)
  	H_beta_random.setZero(Br[0].cols(), Br[0].cols());

  if(Br.size() > 0){
    D = duplicatematrix(Br[0].cols());
    Dd = D.cast <double> ();
  }

  if(Br.size() > 0){
    if(init_list.containsElementNamed("Sigma"))
      Sigma     =  Rcpp::as< Eigen::MatrixXd > (init_list["Sigma"]) ;
    else
      Sigma.setIdentity(Br[0].cols(), Br[0].cols());

    Sigma_vech = vech(Sigma);
    npars += Sigma_vech.size();
    UUt.setZero(Sigma.cols() * Sigma.rows());
    dSigma_vech.setZero(Sigma_vech.size());
    dSigma_vech_old.setZero(Sigma_vech.size());
    invSigma  = Sigma.inverse();
    if( init_list.containsElementNamed("U" ))
      U = Rcpp::as< Eigen::MatrixXd > (init_list["U"]);
    else
      U.setZero(Br[0].cols(), Br.size());
  }

  Sigma_epsilon = 0;
  if(init_list.containsElementNamed("Sigma_epsilon"))
  	Sigma_epsilon  =1;

}



void NormalMixedEffect::sampleU2(const int i,
                                const Eigen::VectorXd& res,
                                const Eigen::VectorXd& iV,
                                const double log_sigma2_noise //= 0
                                )
{
    if(Br.size() == 0)
      return;


    Eigen::VectorXd b   = exp( - log_sigma2_noise) * (Br[i].transpose() * iV.cwiseProduct(res));
    Eigen::MatrixXd Q   = exp( - log_sigma2_noise)  * Br[i].transpose() * iV.asDiagonal() * Br[i];
    Q        +=  invSigma;
    U.col(i) =   sample_Nc(b, Q);
    if(Sigma_epsilon){
      Eigen::VectorXd Z = Rcpp::as< Eigen::VectorXd >(Rcpp::rnorm( U.col(i).size()) );
      Z.array() *= beta_random.array().abs() * 1e-4 + 1e-14;
      U.col(i) += Z;
    }
}


void NormalMixedEffect::sampleU2_par(const int i,
                                 const Eigen::VectorXd& res,
                                 const Eigen::VectorXd& iV,
                                 std::mt19937 & random_engine,
                                 const double log_sigma2_noise)
{
  if(Br.size() == 0)
    return;


  Eigen::VectorXd b   = exp( - log_sigma2_noise) * (Br[i].transpose() * iV.cwiseProduct(res));
  Eigen::MatrixXd Q   = exp( - log_sigma2_noise)  * Br[i].transpose() * iV.asDiagonal() * Br[i];
  Q        +=  invSigma;
  U.col(i) =   sample_Nc_par(b, Q,random_engine);
  if(Sigma_epsilon){
    std::normal_distribution<double> normal;
    Eigen::VectorXd Z;
    Z.setZero(U.col(i).size());
    for(int j =0; j < U.col(i).size(); j++)
      Z[j] =  normal(random_engine);

    Z.array() *= beta_random.array().abs() * 1e-4 + 1e-14;
    U.col(i) += Z;
  }
}

void NormalMixedEffect::simulate( )
{
	if(Br.size() == 0)
      return;

	Eigen::VectorXd b; b.setZero(U.rows());
	for(int i = 0; i < U.cols(); i++)
		U.col(i) =   sample_Nc(b, invSigma);
}

void NormalMixedEffect::simulate(std::vector< Eigen::VectorXd > & Y )
{
	if(Bf.size() >0)
	{
		for(int i = 0; i < Bf.size(); i++)
			Y[i] += Bf[i]* beta_fixed;
	}
	if(Br.size() == 0)
      return;

	Eigen::VectorXd b; b.setZero(U.rows());
	for(int i = 0; i < U.cols(); i++) {
		U.col(i) =   sample_Nc(b, invSigma);
		Y[i] += Br[i]* (U.col(i) + beta_random);
	}
}

void NormalMixedEffect::simulate(Eigen::VectorXd  & Y, const int i )
{
	if(Bf.size() >0)
	{
		Y += Bf[i]* beta_fixed;
	}
	if(Br.size() == 0)
      return;

	Eigen::VectorXd b; b.setZero(U.rows());
	U.col(i) =   sample_Nc(b, invSigma);
	Y += Br[i] * (U.col(i) + beta_random);
}


void NormalMixedEffect::sampleU(const int i,
                                const Eigen::VectorXd& res,
                                const double log_sigma2_noise )
{
    if(Br.size() == 0)
      return;

    Eigen::VectorXd b   = exp( - log_sigma2_noise) * (Br[i].transpose() * res);
    Eigen::MatrixXd Q   = exp( - log_sigma2_noise)  * Br[i].transpose() * Br[i];
    Q        +=  invSigma;
    U.col(i) =   sample_Nc(b, Q);
    if(Sigma_epsilon){
      Eigen::VectorXd Z = Rcpp::as< Eigen::VectorXd >(Rcpp::rnorm( U.col(i).size()) );
      Z.array() *= beta_random.array().abs() * 1e-4 + 1e-14;
      U.col(i) += Z;
    }
}

void NormalMixedEffect::sampleU_par(const int i,
                                const Eigen::VectorXd& res,
                                const double log_sigma2_noise,
                                std::mt19937 & random_engine)
{
  if(Br.size() == 0)
    return;

  Eigen::VectorXd b   = exp( - log_sigma2_noise) * (Br[i].transpose() * res);
  Eigen::MatrixXd Q   = exp( - log_sigma2_noise)  * Br[i].transpose() * Br[i];
  Q        +=  invSigma;
  U.col(i) =   sample_Nc_par(b, Q,random_engine);
  if(Sigma_epsilon){
    std::normal_distribution<double> normal;
    Eigen::VectorXd Z;
    Z.setZero(U.col(i).size());
    for(int j =0; j < U.col(i).size(); j++)
      Z[j] =  normal(random_engine);

    Z.array() *= beta_random.array().abs() * 1e-4 + 1e-14;
    U.col(i) += Z;
  }
}

void NormalMixedEffect::gradient2(const int i,
                                 const Eigen::VectorXd& res,
                                 const Eigen::VectorXd& iV,
                                 const double log_sigma2_noise,  // = 0
                                 const double EiV // = 0
                                 )
{
    counter++;
    Eigen::VectorXd res_  = res;
    if(Br.size() > 0){
      res_ -= Br[i] * U.col(i);
      Eigen::MatrixXd UUT = U.col(i) * U.col(i).transpose();
      UUt += vec( UUT);
      grad_beta_r  += exp( - log_sigma2_noise) * (Br[i].transpose() * iV.cwiseProduct(res_));
      grad_beta_r2 += (invSigma * U.col(i));
      //H_beta_random +=   exp( - log_sigma2_noise) * (Br[i].transpose() * iV.asDiagonal()* Br[i]);
      H_beta_random +=  EiV *exp( - log_sigma2_noise) * (Br[i].transpose()* Br[i]);
    }
    if(Bf.size() > 0){
      grad_beta_f   +=  exp( - log_sigma2_noise) * (Bf[i].transpose() *  iV.cwiseProduct(res_));
      //H_beta_fixed  +=  exp( - log_sigma2_noise) * (Bf[i].transpose() *iV.asDiagonal()* Bf[i]);
      H_beta_fixed  +=  EiV * exp( - log_sigma2_noise) * (Bf[i].transpose() * Bf[i]);
    }
}
void NormalMixedEffect::gradient(const int i,
                                 const Eigen::VectorXd& res,
                                 const double log_sigma2_noise)
{
    counter++;
    Eigen::VectorXd res_  = res;
    if(Br.size() > 0){
      res_ -= Br[i] * U.col(i);
      Eigen::MatrixXd UUT = U.col(i) * U.col(i).transpose();
      UUt += vec( UUT);
      grad_beta_r  += exp( - log_sigma2_noise) * (Br[i].transpose() * res_);
      grad_beta_r2 +=  (invSigma * U.col(i));
      H_beta_random +=  exp( - log_sigma2_noise) * (Br[i].transpose() * Br[i]);

    }
    if(Bf.size() > 0){
      grad_beta_f   +=  exp( - log_sigma2_noise) * (Bf[i].transpose() * res_);
      H_beta_fixed  +=  exp( - log_sigma2_noise) * (Bf[i].transpose() * Bf[i]);
    }
}
void NormalMixedEffect::step_theta(const double stepsize,
								   const double learning_rate,
								   const double polyak_rate)
{
  if(Br.size() > 0){
    step_beta_random(stepsize, learning_rate);
    step_Sigma(stepsize, learning_rate);
  }
  if(Bf.size() > 0)
    step_beta_fixed(stepsize, learning_rate);

if(store_param){
    if(Bf.size() > 0){
      if(vec_counter == 0 || polyak_rate == -1)
      	betaf_vec.row(vec_counter)  = beta_fixed;
      else{
      	betaf_vec.row(vec_counter).array()  = polyak_rate * beta_fixed.array();
      	betaf_vec.row(vec_counter).array()  += (1. - polyak_rate) * betaf_vec.row(vec_counter - 1).array();
      }
    }
    if(Br.size() > 0)
    {
    	Eigen::Map<Eigen::VectorXd> temp(Sigma.data(),Sigma.size());
    	if(vec_counter == 0 || polyak_rate == -1){	
      		betar_vec.row(vec_counter)  = beta_random;
      		Sigma_vec.row(vec_counter)  = temp;
      	
      	}else{
      		betar_vec.row(vec_counter).array()   = polyak_rate * beta_random.array();
      		betar_vec.row(vec_counter).array()  += (1 - polyak_rate) * betar_vec.row(vec_counter - 1).array();
      	
      		Sigma_vec.row(vec_counter).array()   = polyak_rate * temp.array();
      		Sigma_vec.row(vec_counter).array()  +=  + (1 - polyak_rate) * Sigma_vec.row(vec_counter - 1).array();
      	}	
    }	
    vec_counter++;
  }

  clear_gradient();
  counter = 0;
}
void NormalMixedEffect::step_beta_fixed(const double stepsize,const double learning_rate)
{
	dbeta_f_old.array() *= learning_rate;
	dbeta_f_old += H_beta_fixed.ldlt().solve(grad_beta_f);
    beta_fixed += stepsize *  dbeta_f_old;
    H_beta_fixed.setZero(Bf[0].cols(), Bf[0].cols());

}
void NormalMixedEffect::step_beta_random(const double stepsize,const double learning_rate)
{
	
	dbeta_r_old.array() *= learning_rate;
	dbeta_r_old += 0.5 *  H_beta_random.ldlt().solve(grad_beta_r);
	dbeta_r_old += 0.5 * (Sigma * grad_beta_r2)/ counter;
    beta_random += stepsize * dbeta_r_old;
    grad_beta_r2.setZero(Br[0].cols());
    H_beta_random.setZero(Br[0].cols(), Br[0].cols());
}

void NormalMixedEffect::step_Sigma(const double stepsize,const double learning_rate)
{
    double pos_def = 0;
  iSkroniS = kroneckerProduct(invSigma, invSigma);
  UUt -= counter*vec(Sigma);
  dSigma_vech = 0.5 * Dd.transpose() * iSkroniS * UUt;
  ddSigma = 0.5 * counter * Dd.transpose() * iSkroniS * Dd;
  dSigma_vech = ddSigma.ldlt().solve(dSigma_vech);
  dSigma_vech_old *= learning_rate;
  dSigma_vech_old += dSigma_vech;
  
  double stepsize_temp  = stepsize;
  while(pos_def <= 0){
    Eigen::VectorXd Sigma_vech_temp = Sigma_vech;
    Sigma_vech_temp += stepsize_temp * dSigma_vech_old;
    Eigen::VectorXd temp = Dd*Sigma_vech_temp;
    Sigma = veci(temp, Sigma.rows(), Sigma.cols());
    stepsize_temp *= 0.5;
    SelfAdjointEigenSolver<MatrixXd> eig(Sigma,EigenvaluesOnly);
    pos_def = eig.eigenvalues().minCoeff();
    if(stepsize_temp <= 1e-16){
        Rcpp::Rcout << "Sigma = \n" << Sigma << "\n";
        Rcpp::Rcout << "pos_def = " << pos_def <<"\n";
        throw("in midexeffect not pos def \n");
    }


  }

    UUt.setZero(Sigma.cols() * Sigma.rows());
    invSigma  = Sigma.inverse();
    Sigma_vech = vech(Sigma);
}


void NormalMixedEffect::clear_gradient()
{

	if(Bf.size() > 0)
		grad_beta_f.setZero(Bf[0].cols());

	if(Br.size() > 0){
		dSigma_vech.setZero(Sigma_vech.size());
		grad_beta_r.setZero(Br[0].cols());
	}
}

Eigen::VectorXd NormalMixedEffect::get_gradient()
{
	Eigen::VectorXd g(npars);
	int start = 0;
	if(Bf.size() > 0 ){
		g.segment(0, Bf[0].cols()) = grad_beta_f;
		start += Bf[0].cols();
	}
	if(Br.size() > 0)
	{
		g.segment(start, Br[0].cols()) = grad_beta_r;
		start += Br[0].cols();

  		iSkroniS = kroneckerProduct(invSigma, invSigma);
  		Eigen::VectorXd UUt_temp = UUt;
  		UUt_temp -= counter*vec(Sigma);
  		dSigma_vech = 0.5 * Dd.transpose() * iSkroniS * UUt_temp;
		g.segment(start, dSigma_vech.size()) = dSigma_vech;
	}
	return(g);
}
