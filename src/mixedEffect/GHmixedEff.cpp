#include "MixedEffect.h"
#include "error_check.h"
#include <chrono>

void GHMixedEffect::get_param_names(Rcpp::StringVector & names){

  MixedEffect::get_param_names(names);
  if(Br.size() > 0 )
  {
    for (int i = 0; i < Br[0].cols(); ++i)
      names.push_back("mu_random_" + std::to_string(i+1));

    int n_s = Br[0].cols() * (Br[0].cols() +1) /2;
    for (int i = 0; i < n_s; ++i)
      names.push_back("Sigma_random_" + std::to_string(i+1));
  }


}

void GHMixedEffect::printIter()
{
  
  Rcpp::Rcout << "sample_MALA = " << sample_MALA << "\n";
  if(Bf.size() > 0)
    Rcpp::Rcout << "beta_f = " << beta_fixed.transpose() << "\n";


  if(Br.size() > 0){
    Rcpp::Rcout << "mu = " << mu.transpose() << "\n";
    Rcpp::Rcout << "beta_r = " << beta_random.transpose() << "\n";
    Rcpp::Rcout << "D(sigma) = " << Sigma.diagonal().transpose() << "\n";
  }
}

void GHMixedEffect::setupStoreTracj(const int Niter)
{

  if(Bf.size() > 0)
    betaf_vec.resize(Niter, Bf[0].cols());

  if(Br.size() > 0){
    betar_vec.resize(Niter, Br[0].cols());
    mu_vec.resize(Niter, Br[0].cols());
    Sigma_vec.resize(Niter, pow(Br[0].cols(), 2));
  }


  vec_counter = 0;
  store_param = 1;
}

void GHMixedEffect::get_param(std::vector<double> & param_in ){

  MixedEffect::get_param(param_in);
  if(Br.size() > 0 )
  {
    for (int i = 0; i < Br[0].cols(); ++i)
      param_in.push_back(mu[i]);

    int n_s = Br[0].cols() * (Br[0].cols() +1) /2;
    for (int i = 0; i < n_s; ++i)
      param_in.push_back(Sigma_vech[i]);
  }
}



Rcpp::List GHMixedEffect::toList()
{
  Rcpp::List out;
  out["B_fixed"]          = Bf;
  out["B_random"]          = Br;
  out["beta_random"] = beta_random;
  out["beta_fixed"]  = beta_fixed;
  Eigen::VectorXd temp_r = - beta_random_constrainted;
  temp_r.array() += 1;
  Eigen::VectorXd temp_f = - beta_fixed_constrainted;
  temp_f.array() += 1;
  out["beta_random_constrained"] = temp_r;
  out["beta_fixed_constrained"]  = temp_f;
  out["Sigma"]  = Sigma;
  out["U"]      = U;
  out["V"]      = V;
  out["mu"]     = mu;
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
      out["beta_random"] = betar_vec.row(mu_vec.rows() - 1);
    out["mu_vec"]      = mu_vec;
    if(betar_vec.rows() > 1)
      out["mu"]          = mu_vec.row(mu_vec.rows() - 1);
    out["Sigma_vec"]   = Sigma_vec;
    Eigen::VectorXd temp = Sigma_vec.row(mu_vec.rows() - 1);
    if(betar_vec.rows() > 1)
      out["Sigma"]       = veci(temp, Sigma.rows(), Sigma.cols());
  }

  }
  return(out);
}

void GHMixedEffect::initFromList(Rcpp::List const &init_list)
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


    if(beta_fixed.size() != Bf[0].cols())
    {
      Rcpp::Rcout << "\nERROR:\n ";
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

    dbeta_f_old.setZero(Bf[0].cols());
    H_beta_fixed.setZero(Bf[0].cols(), Bf[0].cols());
    npars += Bf[0].cols();

  }else{ Bf.resize(0);}
  count = 0;

  if(init_list.containsElementNamed("B_random"))
  {
    Rcpp::List Br_list = init_list["B_random"];
    Br.resize(Br_list.length());
    for( Rcpp::List::iterator it = Br_list.begin(); it != Br_list.end(); ++it ) {
      Br[count++] = Rcpp::as < Eigen::MatrixXd >( it[0]);
    }
    if(init_list.containsElementNamed("beta_random"))
        beta_random = Rcpp::as < Eigen::VectorXd >( init_list["beta_random"]);
    else
      beta_random.setZero(Br[0].cols());

    if(beta_random.size() != Br[0].cols())
    {
      Rcpp::Rcout << "ERROR:\n ";
      Rcpp::Rcout << "beta_random.size = " << beta_random.size() << "\n";
      Rcpp::Rcout << "B_random.cols    = " << Br[0].cols() << "\n";
      throw("input error\n");
    }
    npars += Br[0].cols();
    grad_beta_r.setZero(Br[0].cols());
    grad_beta_r2.setZero(Br[0].cols());
    term2_mu.setZero(Br[0].cols());
    term1_mu = 0;
    if(init_list.containsElementNamed("beta_random_constrained")){
        beta_random_constrainted = Rcpp::as < Eigen::VectorXd >( init_list["beta_random_constrained"]);
        beta_random_constrainted =  1 - beta_random_constrainted.array();
    }else{
        beta_random_constrainted.setOnes(Br[0].cols());;
    }

    dbeta_r_old.setZero(Br[0].cols());
    if(Br.size() > 0){
      H_beta_random.setZero(Br[0].cols(), Br[0].cols());
      D = duplicatematrix(Br[0].cols());
      Dd = D.cast <double> ();
    }

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

    iSkroniS = kroneckerProduct(invSigma, invSigma);
    if( init_list.containsElementNamed("U" ))
      U = Rcpp::as< Eigen::MatrixXd > (init_list["U"]);
    else
      U.setZero(Br[0].cols(), Br.size());

    if( init_list.containsElementNamed("mu" ))
      mu = Rcpp::as< Eigen::VectorXd > (init_list["mu"]) ;
    else
      mu.setZero(Br[0].cols(), 1);

    term1 = 0.;
    term2 = 0.;
    gradMu.setZero(Br[0].cols(), 1);
    gradMu_old.setZero(Br[0].cols(), 1);
    npars += Br[0].cols();
    gradMu_2.setZero(Br[0].cols(), 1);

    
    if(init_list.containsElementNamed("seed"))
      rgig.seed( Rcpp::as< unsigned long > (init_list["seed"]));
    else
      rgig.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());

    Sigma_epsilon = 0;
    if( init_list.containsElementNamed("V" ))
      V = Rcpp::as< Eigen::VectorXd > (init_list["V"]) ;
    else{
       V.setZero(Br.size());
    }

    SelfAdjointEigenSolver<MatrixXd> eig(Sigma,EigenvaluesOnly);
    double pos_def = eig.eigenvalues().minCoeff();
    if(pos_def/mu.array().abs().maxCoeff() < 0.0005)
      sample_MALA = 1;
    else
      sample_MALA = 0;

  }else{ Br.resize(0);}



}

void GHMixedEffect::sampleV(const int i)
{
  if(Br.size() == 0)
      return;
      //  GIG (p, a, b)
      double p  = get_p_GIG() - 0.5 * Br[0].cols();
      Eigen::VectorXd U_ = U.col(i) + mu;
      double b  =  U_.transpose() * invSigma * U_;
      b += get_b_GIG();
      V(i) = rgig.sample(p, get_a_GIG(), b);
      
}

void GHMixedEffect::simulate( )
{
  if(Br.size() == 0)
      return;

  Eigen::VectorXd b; b.setZero( U.rows());
  for(int i = 0; i < U.cols(); i++) {
    V(i)     =   rgig.sample(get_p_GIG(), 
                             get_a_GIG(),
                             get_b_GIG());

    U.col(i) =   sample_Nc(b, invSigma/V(i));
    U.col(i) += -mu + mu * V(i);
  }
}


void GHMixedEffect::simulate(std::vector< Eigen::VectorXd > & Y )
{
  if(Bf.size() >0)
  {
    for(int i = 0; i < Bf.size(); i++)
      Y[i] += Bf[i]* beta_fixed;
  }
  if(Br.size() == 0)
      return;

  Eigen::VectorXd b; b.setZero( U.rows());
  for(int i = 0; i < U.cols(); i++) {
    V(i)     =   rgig.sample(get_p_GIG(),
                             get_a_GIG(),
                             get_b_GIG());
    U.col(i) =   sample_Nc(b, invSigma/V(i));
    U.col(i) += -mu + mu * V(i);
    Y[i] += Br[i]* (U.col(i) + beta_random);
  }
}

void GHMixedEffect::simulate(Eigen::VectorXd  & Y, const int i )
{

  if(Bf.size() >0)
  {
    Y += Bf[i]* beta_fixed;
  }
  if(Br.size() == 0)
      return;

  Eigen::VectorXd b; b.setZero(U.rows());
  
  V(i)     =   rgig.sample(get_p_GIG(),
                           get_a_GIG(),
                           get_b_GIG());
  U.col(i) =   sample_Nc(b, invSigma/V(i));
  U.col(i) += -mu + mu * V(i);
  Y += Br[i] * (U.col(i) + beta_random);
}



void GHMixedEffect::sampleU_MALA_(const int i,
                                  const Eigen::VectorXd & res,
                                  const Eigen::VectorXd & b_prior,
                                  const Eigen::MatrixXd   & Q_noise)
{
    Eigen::MatrixXd Q_prior   =   Br[i].transpose() * (Q_noise * Br[i]);
    double scale = 2.48/pow(U.col(i).size(),0.33);

    Eigen::VectorXd dU;
    Eigen::MatrixXd ddU;
    dU_ddU_GH( dU,
               ddU,
               U.col(i),
               invSigma,
               -mu,
                mu,
                get_p_GIG(),
                get_a_GIG(),
                get_b_GIG(),
                res,
                Q_noise,
                Br[i]);
    ddU.array() /= scale;

    Eigen::LLT<Eigen::MatrixXd> chol_ddU(ddU);
    Eigen::VectorXd U_mean = U.col(i) - 0.5*chol_ddU.solve(dU);

    Eigen::VectorXd Z = Rcpp::as< Eigen::VectorXd >(Rcpp::rnorm( U.col(i).size()) );
    Eigen::VectorXd Ustar = U_mean + chol_ddU.matrixU().solve(Z);

    Eigen::VectorXd dUstar;
    Eigen::MatrixXd ddUstar;
    dU_ddU_GH( dUstar,
               ddUstar,
               Ustar,
               invSigma,
               -mu,
                mu,
                get_p_GIG(),
                get_a_GIG(),
                get_b_GIG(),
                res,
                Q_noise,
                Br[i]);
    ddUstar.array() /= scale;

    Eigen::LLT<Eigen::MatrixXd> chol_ddUstar(ddUstar);
    Eigen::VectorXd Ustar_mean = Ustar - 0.5*chol_ddUstar.solve(dUstar);

    double qstar = -0.5 * (Ustar - U_mean).transpose()* ( ddU * (Ustar - U_mean));
    Eigen::MatrixXd L_Q = chol_ddUstar.matrixL();
    qstar += L_Q.diagonal().array().log().sum();
    double q0 = -0.5 * (U.col(i) - Ustar_mean).transpose()* ( ddU * (U.col(i) - Ustar_mean));
    L_Q = chol_ddU.matrixL();
    q0 += L_Q.diagonal().array().log().sum();

    double loglik = logdensity( U.col(i));
    loglik -= 0.5 * U.col(i).transpose()*Q_prior*U.col(i) - U.col(i).dot(b_prior);


    double loglikstar = logdensity( Ustar);
    loglikstar -= 0.5 * Ustar.transpose()*Q_prior*Ustar - Ustar.dot(b_prior);

    double alpha = loglikstar - qstar + q0 - loglik;
    if(log( Rcpp::as<double>( Rcpp::runif( 1))) < alpha)
    {
      U.col(i) = Ustar;
      accept_MALA++;
    }
      count_MALA++;
}

