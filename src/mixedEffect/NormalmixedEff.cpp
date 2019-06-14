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
  weight_total = 0;
  calc_grad = 1;
  //dlog_sigma2  = 0;
  //ddlog_sigma2 = 0;
}

NormalMixedEffect::~NormalMixedEffect(){

  for(int i=0;i<Bf.size();i++) {
    Bf[i].resize(0,0);
  }
  Bf.resize(0);
}

void NormalMixedEffect::get_param(std::vector<double> & param_in){

  MixedEffect::get_param(param_in);
  if(Br.size() > 0 )
  {
    int n_s = Br[0].cols() * (Br[0].cols() +1) /2;
    for (int i = 0; i < n_s; ++i)
      param_in.push_back(Sigma_vech[i]);
  }


}

void NormalMixedEffect::get_param_names(Rcpp::StringVector & names){

  MixedEffect::get_param_names(names);
  if(Br.size() > 0 )
  {

    int n_s = Br[0].cols() * (Br[0].cols() +1) /2;
    for (int i = 0; i < n_s; ++i)
      names.push_back("Sigma_random_" + std::to_string(i+1));

  }
}
void NormalMixedEffect::printIter()
{
	if(Bf.size() > 0)
		Rcpp::Rcout << "beta_f = " << beta_fixed.transpose() << "\n";


	if(Br.size() > 0){
		Rcpp::Rcout << "beta_r = " << beta_random.transpose() << "\n";
    Rcpp::Rcout << "D(sigma) = " << Sigma.diagonal().transpose() << "\n";
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
  Eigen::VectorXd temp_r = - beta_random_constrainted;
  temp_r.array() += 1;
  Eigen::VectorXd temp_f = - beta_fixed_constrainted;
  temp_f.array() += 1;
  
  Rcpp::List out;
  if(Bf.size() > 0){
    out["B_fixed"]          = Bf;
    out["beta_fixed"]  = beta_fixed;
    out["beta_fixed_constrained"]  = temp_f;
  }
  if(Br.size() > 0){
    out["B_random"]          = Br;
    out["beta_random"] = beta_random;
    out["beta_random_constrained"] = temp_r;
    out["Sigma"]       = Sigma;
    out["U"]           = U;
  }
  
  
  out["noise"]       = noise;
  out["Sigma_epsilon"]       = Sigma_epsilon;
  out["Cov_theta"]   = Cov_theta;
  if(store_param){
    if(Bf.size() > 0){
		  out["betaf_vec"] = betaf_vec;
      if(betaf_vec.rows() > 1)
		    out["beta_fixed"] = betaf_vec.row(betaf_vec.rows() - 1);
	}

   if(Br.size() > 0){
		  out["betar_vec"]   = betar_vec;
      if(betar_vec.rows() > 1)
		    out["beta_random"] = betar_vec.row(betar_vec.rows() - 1);
		  out["Sigma_vec"]   = Sigma_vec;
		  Eigen::VectorXd temp = Sigma_vec.row(betar_vec.rows() - 1);
      if(betar_vec.rows() > 1)
        out["Sigma"]       = veci(temp, Sigma.rows(), Sigma.cols());
	 }
  }
  return(out);
}

void NormalMixedEffect::initFromList(Rcpp::List const &init_list)
{

  int count =0;
  n_f = 0;
  if(init_list.containsElementNamed("B_fixed"))
  {
    Rcpp::List Bf_list = init_list["B_fixed"];
    Bf.resize(Bf_list.length());
    for( Rcpp::List::iterator it = Bf_list.begin(); it != Bf_list.end(); ++it ) {
      Bf[count++] = Rcpp::as < Eigen::MatrixXd >( it[0]);
    }
    grad_beta_f.setZero(Bf[0].cols());
    n_f = Bf[0].cols();
    if(init_list.containsElementNamed("beta_fixed"))
      beta_fixed = Rcpp::as < Eigen::VectorXd >( init_list["beta_fixed"]);
    else
      beta_fixed.setZero(Bf[0].cols());
    if(beta_fixed.size() != Bf[0].cols())
    {
      Rcpp::Rcout << "\nERROR: \n";
      Rcpp::Rcout << "beta_fixed.size = " << beta_fixed.size() << "\n";
      Rcpp::Rcout << "B_fixed.cols    = " << Bf[0].cols() << "\n";
      throw("input error\n");
    }

    if(init_list.containsElementNamed("beta_fixed_constrained")){
        beta_fixed_constrainted = Rcpp::as < Eigen::VectorXd >( init_list["beta_fixed_constrained"]);
        beta_fixed_constrainted =  1 - beta_fixed_constrainted.array();
    }else{
        beta_fixed_constrainted.setOnes(Bf[0].cols());;
    }

    npars += Bf[0].cols();
    H_beta_fixed.setZero(Bf[0].cols(), Bf[0].cols());
	  dbeta_f_old.setZero(Bf[0].cols());
  }else{ Bf.resize(0);}
  count = 0;
  n_r = 0;

  if(init_list.containsElementNamed("B_random"))
  {
    Rcpp::List Br_list = init_list["B_random"];
    Br.resize(Br_list.length());
    for( Rcpp::List::iterator it = Br_list.begin(); it != Br_list.end(); ++it ) {
      Br[count++] = Rcpp::as < Eigen::MatrixXd >( it[0]);
    }

    grad_beta_r.setZero(Br[0].cols());
    grad_beta_r2.setZero(Br[0].cols());

    n_r = Br[0].cols();
    V.setOnes(Br.size());
    if(init_list.containsElementNamed("beta_random"))
      beta_random = Rcpp::as < Eigen::VectorXd >( init_list["beta_random"]);
    else
      beta_random.setZero(Br[0].cols());
    if(beta_random.size() != Br[0].cols())
    {
      Rcpp::Rcout << "\nERROR:\n ";
      Rcpp::Rcout << "beta_random.size = " << beta_random.size() << "\n";
      Rcpp::Rcout << "B_random.cols    = " << Br[0].cols() << "\n";
      throw("input error\n");
    }
    if(init_list.containsElementNamed("beta_random_constrained")){
        beta_random_constrainted = Rcpp::as < Eigen::VectorXd >( init_list["beta_random_constrained"]);
        beta_random_constrainted =  1 - beta_random_constrainted.array();
    }else{
        beta_random_constrainted.setOnes(Br[0].cols());;
    }

	   dbeta_r_old.setZero(Br[0].cols());
    npars += Br[0].cols();
  }else{ Br.resize(0);}

  if(Br.size() > 0)
  	H_beta_random.setZero(Br[0].cols(), Br[0].cols());


  if(Br.size() > 0){
    if(init_list.containsElementNamed("Sigma"))
      Sigma     =  Rcpp::as< Eigen::MatrixXd > (init_list["Sigma"]) ;
    else
      Sigma.setIdentity(Br[0].cols(), Br[0].cols());

    D = duplicatematrix(Br[0].cols());
    Dd = D.cast <double> ();
    invSigma  = Sigma.inverse();
    iSkroniS = kroneckerProduct(invSigma, invSigma);
    ddSigma = 0.5 * weight_total * Dd.transpose() * iSkroniS * Dd;


    Sigma_vech = vech(Sigma);
    npars += Sigma_vech.size();
    UUt.setZero(Sigma.cols() * Sigma.rows());
    grad_Sigma.setZero(Dd.cols());
    dSigma_vech.setZero(Sigma_vech.size());
    dSigma_vech_old.setZero(Sigma_vech.size());
    invSigma  = Sigma.inverse();

    if( init_list.containsElementNamed("U" ))
      U = Rcpp::as< Eigen::MatrixXd > (init_list["U"]);
    else
      U.setZero(Br[0].cols(), Br.size());
    mu0.resize(Br.size());
    for(int i = 0; i < Br.size(); i++){
      mu0[i].setZero(Br[0].cols());
    }
    EU.setZero(Br[0].cols(), Br.size());
  }else{
    mu0.resize(0);
  }

  Sigma_epsilon = 0;
  if(init_list.containsElementNamed("Sigma_epsilon"))
  	Sigma_epsilon  =1;

  grad_beta.setZero(n_f + n_r);
  H_beta.setZero(n_f + n_r, n_r + n_f);
  nfr = n_f + n_r;
  
  Hessian.setZero(nfr, nfr);
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
    EU.col(i) = Q.ldlt().solve(b);
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
  EU.col(i) = Q.ldlt().solve(b);
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
    EU.col(i) = Q.ldlt().solve(b);
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
  EU.col(i) = Q.ldlt().solve(b);
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
/*
  Caculating the gradient for measurement error which is Normal
  two special casses are of importance
  single measurement error sigma or vector type
*/
void NormalMixedEffect::gradient(const int i,
                                 const Eigen::VectorXd& res,
                                 const double log_sigma2_noise,
                                 const double weight,
                                 const int use_EU// =1
                                 )
{

    counter++;
    Eigen::VectorXd res_  = res;
    Eigen::VectorXd isigmas2;
    gradientVec.setZero(npars);
    if(Br.size() > 0){
      if(use_EU)
        res_ -= Br[i] * EU.col(i);
      else
        res_ -= Br[i] * U.col(i);
    
      Eigen::MatrixXd UUT = U.col(i) * U.col(i).transpose();
      Eigen::VectorXd  dSigma_ = 0.5 * Dd.transpose() * iSkroniS * (vec( UUT) - vec(Sigma));
      gradientVec.tail(dSigma_.size()) = dSigma_; 
      UUt += weight * vec( UUT);//gradient.tail(UUt.length());
      if(calc_grad){
        grad_beta_r   += weight * exp( - log_sigma2_noise) * (Br[i].transpose() * res_);
        if(use_EU)
          grad_beta_r2  += weight * (invSigma * EU.col(i));
        else
          grad_beta_r2  += weight * (invSigma * U.col(i));
      }
      
      H_beta_random += weight * exp( - log_sigma2_noise) * (Br[i].transpose() * Br[i]);
      

    }

  if(Bf.size() > 0){
      if(calc_grad)
        grad_beta_f   += weight * exp( - log_sigma2_noise) * (Bf[i].transpose() * res_);
      H_beta_fixed  += weight * exp( - log_sigma2_noise) * (Bf[i].transpose() * Bf[i]);
        
  }
  if(1){
      if(Br.size() > 0){
        
        if(calc_grad){
          gradientVec.segment(n_f, n_r) = exp( - log_sigma2_noise) * (Br[i].transpose() * res_); 
          grad_beta.head(n_r)   += weight * gradientVec.segment(n_f, n_r) ;
          grad_beta_r2         += weight * (invSigma * U.col(i));
        }
      //
      H_beta.topLeftCorner(n_r, n_r)  += weight * exp( - log_sigma2_noise) * (Br[i].transpose() * Br[i]);
        
      }
      if(Bf.size() > 0){
        
        if(calc_grad){
          gradientVec.head(n_f) = exp( - log_sigma2_noise) * (Bf[i].transpose() * res_);
          grad_beta.tail(n_f) += weight * gradientVec.head(n_f);
        }
        H_beta.bottomRightCorner(n_f, n_f)  += weight * exp( - log_sigma2_noise) * (Bf[i].transpose() * Bf[i]);
        
        if(Br.size() > 0){
            H_beta.topRightCorner(n_r, n_f) +=  weight * exp( - log_sigma2_noise) * (Br[i].transpose() * Bf[i]);
            H_beta.bottomLeftCorner(n_f, n_r) +=  weight * exp( - log_sigma2_noise) * (Bf[i].transpose() * Br[i]);
           
          }
      }
    }
     weight_total += weight;

}

void NormalMixedEffect::gradient2(const int i,
                                 const Eigen::VectorXd& res,
                                 const Eigen::VectorXd& iV,
                                 const Eigen::VectorXd& sigmas,  // =0
                                 const double log_sigma2_noise,  // = 0
                                 const double EiV, // = 0
                                 const double weight, // = 1
                                 const int use_EU, // =1,
                                 const int nssigma            //  = 0
                                 )
{
    counter++;
    Eigen::VectorXd res_  = res;
    
    weight_total += weight;

    /*
        Hessian part
    */

    Eigen::VectorXd isigmas2;
    if(nssigma)
      isigmas2 = sigmas.array().inverse().pow(2);
    Eigen::MatrixXd Hbr, Hbf;
    if(Br.size() > 0){
        if(nssigma){
            Hbr = weight * EiV *exp( - log_sigma2_noise) * (Br[i].transpose() * isigmas2.asDiagonal() * Br[i]);
        }else{
            Hbr = weight * EiV *exp( - log_sigma2_noise) * (Br[i].transpose() * Br[i]);
        }
        
        H_beta_random                   += Hbr;
        H_beta.topLeftCorner(n_r, n_r)  += Hbr;
    }

    if(Bf.size() > 0){
      if(nssigma){
        Hbf = weight * EiV * exp( - log_sigma2_noise) * (Bf[i].transpose() * isigmas2.asDiagonal() * Bf[i]);
      }else{
        Hbf = weight * EiV * exp( - log_sigma2_noise) * (Bf[i].transpose() * Bf[i]);
      }
        H_beta_fixed                        += Hbf;
        H_beta.bottomRightCorner(n_f, n_f)  += Hbf;
        if(Br.size() > 0){
          if(nssigma){
            H_beta.topRightCorner(n_r, n_f)   +=  weight * EiV * exp( - log_sigma2_noise) * (Br[i].transpose() * isigmas2.asDiagonal() *  Bf[i]);
          }else{
            H_beta.topRightCorner(n_r, n_f)   +=  weight * EiV * exp( - log_sigma2_noise) * (Br[i].transpose() * Bf[i]);
          }
          H_beta.bottomLeftCorner(n_f, n_r) +=  H_beta.topRightCorner(n_r, n_f).transpose();
        }
    }
    /*
        Gradient part

    */
          Eigen::VectorXd gf, gr;
      if(Br.size() > 0){
        if(use_EU)
          res_ -= Br[i] * EU.col(i);
        else
          res_ -= Br[i] * U.col(i);

        res_ = iV.cwiseProduct(res_);
        Eigen::MatrixXd UUT = U.col(i) * U.col(i).transpose();
        UUt += weight * vec( UUT);
        if(calc_grad){
          gr            = weight * exp( - log_sigma2_noise) * (Br[i].transpose() * res_);
          grad_beta_r         += gr;
          grad_beta.head(n_r) += gr;
          if(use_EU)
            grad_beta_r2 += weight * (invSigma * EU.col(i));
          else
            grad_beta_r2 += weight * (invSigma * U.col(i));
        }
      }else{
        res_ = iV.cwiseProduct(res_);
      }

      if(calc_grad){
        if(Bf.size() > 0){
          gf                   =  weight * exp( - log_sigma2_noise) * (Bf[i].transpose() *  res_);
          grad_beta_f         +=  gf;
          grad_beta.tail(n_f) +=  gf;
        }
      }

}

Eigen::MatrixXd NormalMixedEffect::d2Given2(const int i,
                                 const Eigen::VectorXd& res,
                                 const Eigen::VectorXd& iV,
                                 const double log_sigma2_noise,  // = 0
                                 const double EiV, // = 0
                                 const double weight) // = 1)
{

  Eigen::VectorXd res_  = res;

  int n_s  =0;

  if(Br.size()>0)
     n_s = n_r * (n_r +1) /2;
  Eigen::MatrixXd d2            = Eigen::MatrixXd::Zero(n_s + n_f + n_r + 1,n_s + n_f+n_r + 1);
  /*
  if(Br.size()>0){
    res_ -= Br[i] * U.col(i);
    d2.block(n_r +n_f , n_r +n_f, n_s, n_s)  -=  0.5* weight * Dd.transpose() * iSkroniS * Dd;
    Eigen::MatrixXd UUT = U.col(i) * U.col(i).transpose();
    d2.block(n_r +n_f , n_r +n_f, n_s, n_s)  += weight *  Dd.transpose() * kroneckerProduct(invSigma, invSigma * UUT * invSigma ) * Dd;
  }
  d2(n_s +  n_r +n_f, n_s +  n_r +n_f ) =  3  * weight * exp( - 2   * log_sigma2_noise)  * (res_.array().square() *iV.array()).sum();
  d2(n_s +  n_r +n_f, n_s +  n_r +n_f ) +=  -1 * weight * res_.size()  * exp( - log_sigma2_noise);
  res_ = iV.cwiseProduct(res_);
  if(calc_grad){
    d2.block(  n_f      , n_f       , n_r, n_r)  =  weight * exp( - log_sigma2_noise) * (Br[i].transpose() * iV.asDiagonal() * Br[i]);
    if(Br.size() * Bf.size()>0){
      d2.block(  0      , n_f     , n_f, n_r)  =  weight * exp( - log_sigma2_noise) * (Bf[i].transpose() * iV.asDiagonal() * Br[i]);
      d2.block(n_f      , 0       , n_r, n_f)  =  weight * exp( - log_sigma2_noise) * (Br[i].transpose() * iV.asDiagonal() * Bf[i]);
    }
    if(Bf.size() > 0 )
      d2.block(0      , 0     , n_f, n_f)  =  weight * exp( - log_sigma2_noise) * (Bf[i].transpose() * iV.asDiagonal() * Bf[i]);
  }
  if(Br.size()>0){

    d2.block(n_f              , n_s + n_r +n_f , n_r , 1 )   =  2 * weight * exp( - 1.5 * log_sigma2_noise)  * (Br[i].transpose() * res_);
    d2.block(n_s +  n_r + n_f, n_f              , 1   , n_r ) = d2.block(n_f              , n_s + n_r +n_f , n_r , 1 ).transpose();
  }
  
 if(Bf.size() > 0){
    d2.block(0            , n_s + n_r + n_f , n_f , 1 )   =  2 * weight * exp( - 1.5 * log_sigma2_noise)  * (Bf[i].transpose() * res_);
    d2.block(n_s + n_r + n_f, 0             , 1   , n_f ) = d2.block(0            , n_s + n_r + n_f , n_f , 1 ).transpose();
  }
  */
   d2(n_s +  n_r +n_f, n_s +  n_r +n_f ) +=  1;
  return(d2);
}


Eigen::MatrixXd NormalMixedEffect::d2Given(const int i,
                                           const Eigen::VectorXd& res,
                                           const double log_sigma2_noise,
                                           const double weight)
{

  Eigen::VectorXd res_  = res;

  int n_s  =0;
  if(Br.size()>0)
     n_s = n_r * (n_r +1) /2;
  Eigen::MatrixXd d2            = Eigen::MatrixXd::Zero(n_s + n_f + n_r + 1,n_s + n_f+n_r + 1);
  /*
  if(Br.size()>0){
    res_ -= Br[i] * U.col(i);
    d2.block(n_r +n_f , n_r +n_f, n_s, n_s)  +=  0.5* weight * Dd.transpose() * iSkroniS * Dd;
    Eigen::MatrixXd UUT = U.col(i) * U.col(i).transpose();
    d2.block(n_r +n_f , n_r +n_f, n_s, n_s)  -= weight *  Dd.transpose() * kroneckerProduct(invSigma, invSigma * UUT * invSigma ) * Dd;
  }
  
  if(calc_grad){
    if(Br.size() * Bf.size()>0){
      d2.block(  0      , n_f     , n_f, n_r)  =  weight * exp( - log_sigma2_noise) * (Bf[i].transpose() * Br[i]);
      d2.block(n_f      , 0       , n_r, n_f)  =  weight * exp( - log_sigma2_noise) * (Br[i].transpose() * Bf[i]);
    }
    if(Br.size()>0){
      d2.block(  n_f      , n_f       , n_r, n_r)  =  weight * exp( - log_sigma2_noise) * (Br[i].transpose() * Br[i]);
    }
    if(Bf.size() > 0 )
      d2.block(0      , 0     , n_f, n_f)  =  weight * exp( - log_sigma2_noise) * (Bf[i].transpose() * Bf[i]);
  }
  if(Br.size()>0){

    d2.block(n_f              , n_s + n_r +n_f , n_r , 1 )   =  2 * weight * exp( - 1.5 * log_sigma2_noise)  * (Br[i].transpose() * res_);
    d2.block(n_s +  n_r + n_f, n_f              , 1   , n_r ) = d2.block(n_f              , n_s + n_r +n_f , n_r , 1 ).transpose();
  }
 if(Bf.size() > 0){
    d2.block(0            , n_s + n_r + n_f , n_f , 1 )   =  2 * weight * exp( - 1.5 * log_sigma2_noise)  * (Bf[i].transpose() * res_);
    d2.block(n_s + n_r + n_f, 0             , 1   , n_f ) = d2.block(0            , n_s + n_r + n_f , n_f , 1 ).transpose();
  }
  
  d2(n_s +  n_r +n_f, n_s +  n_r +n_f ) =  3  * weight * exp( - 2   * log_sigma2_noise)  * res_.array().square().sum();
  d2(n_s +  n_r +n_f, n_s +  n_r +n_f ) +=  -1 * weight * res_.size()  * exp( - log_sigma2_noise);
  */
   d2(n_s +  n_r +n_f, n_s +  n_r +n_f ) +=  1;
  return(d2);
}



void NormalMixedEffect::step_theta(const double stepsize,
								   const double learning_rate,
								   const double polyak_rate,
								   const int burnin)
{
  if(0){
    if(Br.size() > 0){
      step_beta_random(stepsize, learning_rate,0);
      step_Sigma(stepsize, learning_rate,0);
    }
    if(Bf.size() > 0)
      step_beta_fixed(stepsize, learning_rate,0);
}else{
    step_beta(stepsize, learning_rate,0);
    if(0){
      Eigen::VectorXd grad0(nfr);
      grad0 << grad_beta_f, grad_beta_r;
      Eigen::VectorXd step0  = -Hessian.ldlt().solve(grad0);
      Eigen::VectorXd betamu(nfr);
      betamu << beta_fixed, beta_random ;
      Eigen::VectorXd step1=  betamu + stepsize * step0;
      beta_fixed   = step1.head(beta_fixed.size());
      beta_random  = step1.tail(beta_random.size());
      Hessian *= 0.4;
      grad_beta_f *= 0.;
      grad_beta_r *= 0.;
      grad_beta   *= 0.;
    
  }
if(Br.size() > 0)
      step_Sigma(stepsize, learning_rate,0);
}
  weight_total = 0;
  clear_gradient();

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

  counter = 0;
}

 void NormalMixedEffect::add_gradient(Eigen::VectorXd & grad){
      int nF = 0;
      if(Bf.size() > 0){
        grad_beta_f += grad.head(Bf[0].cols());
        grad_beta.tail(Bf[0].cols()) += grad.head(Bf[0].cols());
        nF = Bf[0].cols();
      }
      
      if(Br.size() > 0 ){
        int nR = Br[0].cols();
        int nSigma = nR * (nR +1) /2;
        grad_beta_r += grad.segment(nF, nR);
        grad_beta.head(Br[0].cols()) += grad.segment(nF, nR);
        grad_Sigma  += grad.segment(nF + nR, nSigma);
      }

}

void NormalMixedEffect::step_beta(const double stepsize,const double learning_rate,const int burnin)
{
  grad_beta.tail( n_f) = grad_beta.tail( n_f).cwiseProduct( beta_fixed_constrainted);
  grad_beta.head( n_r) = grad_beta.head( n_r).cwiseProduct(beta_random_constrainted);
  

  Eigen::VectorXd step1;
  step1.setZero(n_f + n_r);
  int n_unconstrained = beta_fixed_constrainted.sum() + beta_random_constrainted.sum() ;
  if( n_unconstrained < n_f + n_r)
  {
      Eigen::VectorXd constrianed;
      constrianed.setZero(n_f + n_r);
      for(int i = 0; i < n_r ; i++)
        constrianed(i) = beta_random_constrainted(i);
      for(int i = 0; i < n_f ; i++)
        constrianed(i + n_r) = beta_fixed_constrainted(i);

      solve_const_x_Ab(step1,
                       constrianed,
                       grad_beta,
                       H_beta);


  }else{
    step1  = H_beta.ldlt().solve(grad_beta);
  }
  if(Bf.size() > 0){
    dbeta_f_old.array() *= learning_rate;
    dbeta_f_old += 0.5 * step1.tail(n_f);
    solve_const_x_Ab(dbeta_f_old,
                     beta_fixed_constrainted,
                     0.5 * grad_beta.tail( n_f),
                     H_beta.bottomRightCorner(n_f, n_f));
    dbeta_f_old =  dbeta_f_old.cwiseProduct(beta_fixed_constrainted);
    beta_fixed.array() += stepsize * dbeta_f_old.array();

  }

  if(Br.size() > 0){
    dbeta_r_old.array() *= learning_rate;
    dbeta_r_old += 0.5 *  step1.head(n_r);
    dbeta_r_old += 0.5 * (Sigma * beta_random_constrainted.cwiseProduct(grad_beta_r2))/ weight_total;
    dbeta_r_old = beta_random_constrainted.cwiseProduct(dbeta_r_old);
    beta_random += stepsize * dbeta_r_old;
    grad_beta_r2.setZero(Br[0].cols());
  }

  H_beta.setZero(n_f + n_r, n_f + n_r);

}
void NormalMixedEffect::step_beta_fixed(const double stepsize,const double learning_rate,const int burnin)
{
	dbeta_f_old.array() *= learning_rate;

  solve_const_x_Ab(dbeta_f_old,
                   beta_fixed_constrainted,
                   grad_beta_f,
                   H_beta_fixed);

  beta_fixed +=  stepsize * (beta_fixed_constrainted.cwiseProduct(dbeta_f_old));
  H_beta_fixed.setZero(Bf[0].cols(), Bf[0].cols());

}
void NormalMixedEffect::step_beta_random(const double stepsize,const double learning_rate,const int burnin)
{
	dbeta_r_old.array() *= learning_rate;

  solve_const_x_Ab(dbeta_r_old,
                   beta_random_constrainted,
                   0.5 * grad_beta_r,
                   H_beta_random);
	dbeta_r_old = 0.5 * (Sigma *  beta_random_constrainted.cwiseProduct(grad_beta_r2))/ weight_total;

  beta_random += stepsize *  beta_random_constrainted.cwiseProduct(dbeta_r_old);
  grad_beta_r2.setZero(Br[0].cols());
  H_beta_random.setZero(Br[0].cols(), Br[0].cols());
}

void NormalMixedEffect::step_Sigma(const double stepsize,const double learning_rate,const int burnin)
{
  double pos_def = 0;
  ddSigma = 0.5 * weight_total * Dd.transpose() * iSkroniS * Dd;
  UUt -= weight_total*vec(Sigma);
  dSigma_vech =  dSigma_vech = grad_Sigma;// 0.5 * Dd.transpose() * iSkroniS * UUt;
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
    iSkroniS = kroneckerProduct(invSigma, invSigma);

}


void NormalMixedEffect::clear_gradient()
{

	if(Bf.size() > 0)
		grad_beta_f.setZero(Bf[0].cols());

	if(Br.size() > 0){
		grad_beta_r.setZero(Br[0].cols());
	}
  grad_beta.setZero(n_f+n_r);
  UUt.setZero(Sigma.cols() * Sigma.rows());
  weight_total = 0.;
  grad_Sigma.array() *= 0.;

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

  		Eigen::VectorXd UUt_temp = UUt;
  		UUt_temp -= weight_total * vec(Sigma);
  		dSigma_vech = 0.5 * Dd.transpose() * iSkroniS * UUt_temp;
		g.segment(start, dSigma_vech.size()) = dSigma_vech;
	}
	return(g);
}
