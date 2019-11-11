#include "measError.h"
#include "error_check.h"


void NormalVarianceMixtureBaseError::printIter()
{
	Rcpp::Rcout << "sigma = " << sigma;
  if(assymetric)
    Rcpp::Rcout << "\nmu = " << mu;

}
void NormalVarianceMixtureBaseError::setupStoreTracj(const int Niter) // setups to store the tracjetory
{
	sigma_vec.resize(Niter);
  mu_vec.resize(Niter);
	vec_counter = 0;
	store_param = 1;
}




NormalVarianceMixtureBaseError::NormalVarianceMixtureBaseError() : MeasurementError(){

  store_param = 0;
  counter     = 0;
  sigma       = 1;
  dsigma      = 0;
  ddsigma     = 0;
  dsigma_old  = 0;
  dmu_old  = 0;
  dmu  = 0;
  common_V    = 0;
  npars       = 0;
  nsSigma     = 0;
  assymetric  = 0;
  VV          = 1;
  ddmu        = 0;
  rgig.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());
}
Rcpp::List NormalVarianceMixtureBaseError::toList()
{
  Rcpp::List out;
  out["noise"]       = noise;
  out["sigma"]       = sigma;
  out["Vs"]          = Vs;
  out["common_V"]    = common_V;
  out["Cov_theta"]   = Cov_theta;
  out["assymetric"]  = assymetric;

  if(assymetric)
    out["mu"] = mu;
  if(store_param){
  	out["sigma_vec"] = sigma_vec;
	  out["sigma"]       = sigma_vec(sigma_vec.size() - 1);
    out["mu_vec"]    = mu_vec;
   }



  return(out);
}
void NormalVarianceMixtureBaseError::initFromList(Rcpp::List const &init_list)
{
  if(init_list.containsElementNamed("sigma"))
    sigma = Rcpp::as < double >( init_list["sigma"]);
  else
  	sigma  = 1.;

   npars += 1;
  if(init_list.containsElementNamed("common_V"))
    common_V = Rcpp::as < int >( init_list["common_V"]);
  else
    common_V  = 0;

  dsigma_old = 0;
 int i = 0;

 if( init_list.containsElementNamed("Vs" )){
 	Rcpp::List Vs_list = init_list["Vs"];
 	Vs.resize(Vs_list.length());
    for( Rcpp::List::iterator it = Vs_list.begin(); it != Vs_list.end(); ++it ) {
      Vs[i++] = Rcpp::as < Eigen::VectorXd >( it[0]);
    }
 }else{
    Rcpp::Rcout << "in NormalVarianceMixtureBaseError:: missing Vs\n";
 	  throw("in NormalVarianceMixtureBaseError::initFromList Vs must be set! \n");
  }

  if( init_list.containsElementNamed("assymetric" ))
    assymetric = Rcpp::as < int>(init_list["assymetric"]);
  
  npars += assymetric;

  mu = 0;
  dmu_old = 0;
  if( init_list.containsElementNamed("mu" ))
      mu = Rcpp::as < double>(init_list["mu"]);
}

double NormalVarianceMixtureBaseError::simulate_V()
{
	return -1;
}
Eigen::VectorXd  NormalVarianceMixtureBaseError::simulate(const Eigen::VectorXd & Y)
{
	Eigen::VectorXd residual =  sigma * (Rcpp::as< Eigen::VectorXd >(Rcpp::rnorm( Y.size()) ));
  
	if(common_V == 0){
		for(int ii = 0; ii < residual.size(); ii++)
	    {
	      double V = simulate_V();
	      residual[ii] *=  sqrt(V);
        if(assymetric)
          residual[ii] += (-1 + V) * mu;
        
	    }

	}else{
		  double V = simulate_V();
	      residual.array() *=  sqrt(V);
        if(assymetric)
          residual.array() += (-1 + V) * mu;
	}

  // 
	return(residual);
}


Eigen::VectorXd  NormalVarianceMixtureBaseError::simulate_par(const Eigen::VectorXd & Y,std::mt19937 & random_engine)
{
  
  std::normal_distribution<double> normal;
  Eigen::VectorXd residual;
  residual.setZero(Y.size());
  for(int j =0; j < Y.size(); j++)
    residual[j] =  sigma*normal(random_engine);

  //Eigen::VectorXd residual =  sigma * (Rcpp::as< Eigen::VectorXd >(Rcpp::rnorm( Y.size()) ));
  if(common_V == 0){
    for(int ii = 0; ii < residual.size(); ii++)
    {
      double V = simulate_V();
      residual[ii] *=  sqrt(V);
      if(assymetric)
        residual[ii] += (-1 + V) * mu;
    }

  }else{
    double V = simulate_V();
    residual.array() *=  sqrt(V);
    if(assymetric)
          residual.array() += (-1 + V) * mu;
  }

  return(residual);
}

Eigen::VectorXd  NormalVarianceMixtureBaseError::simulate_par(const int i,std::mt19937 & random_engine, int nsim){
  std::normal_distribution<double> normal;
  Eigen::VectorXd residual;
  residual.setZero(nsim);
  for(int j =0; j < nsim; j++)
    residual[j] =  sigma*normal(random_engine);
  
  //Eigen::VectorXd residual =  sigma * (Rcpp::as< Eigen::VectorXd >(Rcpp::rnorm( Y.size()) ));
  if(common_V == 0){
    for(int ii = 0; ii < residual.size(); ii++)
    {
      double V = 1;//simulate_V();
      residual[ii] *=  sqrt(V);
    }
    
  }else{
    double V = 1;//simulate_V();
    residual.array() *=  sqrt(V);
  }
  
  return(residual);
}

Eigen::VectorXd  NormalVarianceMixtureBaseError::getSigma(const int i,  const int n )
{ 
  Eigen::VectorXd Result(n);
  Result.setOnes(n);
  Result *= pow(sigma,2);

  if(common_V == 0){
    Result.array() *= Vs[i].array();
  }else{
    Result *= Vs[i];
  }

  return(Result);
}; 


std::vector< Eigen::VectorXd > NormalVarianceMixtureBaseError::simulate(std::vector< Eigen::VectorXd > Y)
{
  
	std::vector< Eigen::VectorXd > residual( Y.size());
	for(int i = 0; i < Y.size(); i++){
		residual[i] =   sigma * (Rcpp::as< Eigen::VectorXd >(Rcpp::rnorm( Y[i].size()) ));
	  if(common_V == 0){
	    for(int ii = 0; ii < residual[i].size(); ii++)
	    {
	      double V = simulate_V();
	      Vs[i][ii] = V;
	      residual[i][ii] *=  sqrt(V);
        if(assymetric)
          residual[i][ii] += (-1 + V) * mu;
	    }
	  } else {
	    double V = simulate_V();
	    for(int ii = 0; ii < residual[i].size(); ii++)
	    {
	      Vs[i][ii] = V;
	      residual[i][ii] *=  sqrt(V);

        if(assymetric)
          residual[i][ii] += (-1 + V) * mu;
	    }
	  }


	}
	return(residual);
}

void NormalVarianceMixtureBaseError::sampleV(const int i, const Eigen::VectorXd& res, int n_s )
{
  // add assymetric
	if(n_s == -1)
		n_s = Vs[i].size();
	if(common_V == 0){
    	for(int j = 0; j < n_s; j++){
        double res2 = res[j];
    	  double mu2  = 0.;
        if(assymetric){
          res2 -=  -mu;
          mu2 = pow(mu/sigma, 2);
        }
    		res2 = pow(res2/sigma, 2);

      	Vs[i](j) = sample_V(res2, -1, mu2);
		}

	} else {
    double tmp; 
    double mu2 = 0;
    if(assymetric){
      tmp = (res.array() + mu).square().sum()/pow(sigma, 2);
      mu2  = pow(mu/sigma, 2)*res.size();
    }else
      tmp = res.array().square().sum()/pow(sigma, 2); 

	  double cv = sample_V(tmp, n_s, mu2);
	  for(int j = 0; j < n_s; j++)
	    Vs[i](j) = cv;
  }
}



void NormalVarianceMixtureBaseError::remove_asym(const int i, Eigen::VectorXd & Y){
  if(assymetric)
    Y.array() -= (-1. + Vs[i].array()) * mu;

}
void NormalVarianceMixtureBaseError::add_asym(const int i, Eigen::VectorXd & Y) {
  if(assymetric)
    Y.array() += (-1. + Vs[i].array()) * mu;
}


void NormalVarianceMixtureBaseError::gradient(const int i,
                                 const Eigen::VectorXd& res,
                                 const double weight)
{
    counter++;
    Eigen::VectorXd res_ = res;
    Eigen::VectorXd iV = Vs[i].cwiseInverse();
    //res_.array() *= iV.array();
    dsigma += weight * (- res.size()/sigma + (res_.array().square()*iV.array()).sum() / pow(sigma, 3) );
    // Expected fisher infromation
    // res.size()/pow(sigma, 2) - 3 * E[res.array().square().sum()] /pow(sigma, 4);
    ddsigma += weight * (- 2 * res.size()/pow(sigma, 2));

    // add assymetric
    if(assymetric){
      dmu    +=  weight * (res_.array()  *  (-iV.array() + 1.)).sum()/pow(sigma, 2);
      ddmu   +=  weight * VV * res_.size()/ pow(sigma,2);  
    }
    weight_total += weight;
}

void NormalVarianceMixtureBaseError::step_theta(const double stepsize,
												const double learning_rate,
												const double polyak_rate,
												const int burnin)
{
  step_sigma(stepsize, learning_rate, burnin, polyak_rate);
  if(assymetric)
    step_mu(stepsize, learning_rate, burnin, polyak_rate);

  NormalVarianceMixtureBaseError::clear_gradient();

	if(store_param)
      vec_counter++;
  	

}

void NormalVarianceMixtureBaseError::step_sigma(const double stepsize, 
                                                const double learning_rate,
                                                const int burnin,
                                                const double polyak_rate)
{

  double sigma_temp = -1;
  dsigma /= ddsigma;
  dsigma_old = learning_rate * dsigma_old + dsigma;
  double stepsize_t = stepsize;
  while(sigma_temp < 0)
  {
    sigma_temp = sigma - stepsize_t * dsigma_old;
    stepsize_t *= 0.5;
    if(stepsize_t <= 1e-16)
        throw("in NormalVarianceMixtureBaseError:: can't make sigma it positive \n");
  }
  sigma = sigma_temp;
  ddsigma = 0;
  if(store_param){
  if(vec_counter == 0 || polyak_rate == -1)
    sigma_vec(vec_counter) = sigma;
  else
    sigma_vec(vec_counter) = polyak_rate * sigma + (1 - polyak_rate) * sigma_vec(vec_counter-1);
    
  }


}

void NormalVarianceMixtureBaseError::step_mu(const double stepsize, 
                                             const double learning_rate,
                                             const int burnin,
                                             const double polyak_rate)
{

  dmu     /= ddmu;
  dmu_old = learning_rate * dmu_old + dmu;
  //Rcpp::Rcout << "mu  = " << mu << "\n";
  //Rcpp::Rcout << "dmu  = " << dmu << "\n";
  mu      = mu + stepsize * dmu_old;
  //Rcpp::Rcout << "mu  = " << mu << "\n";
  if(store_param){
    if(vec_counter == 0 || polyak_rate == -1)
      mu_vec(vec_counter) = mu;
    else
      mu_vec(vec_counter) = polyak_rate * mu + (1 - polyak_rate) * mu_vec(vec_counter-1); 
  }

}


void NormalVarianceMixtureBaseError::clear_gradient()
{
	dsigma = 0;
  dmu    = 0;
  ddmu    = 0;
  weight_total = 0;
}

Eigen::VectorXd NormalVarianceMixtureBaseError::get_gradient()
{
	Eigen::VectorXd g(npars);
	g(0) = dsigma;
  if(assymetric)
    g(1) = dmu;
	return(g);
}
Eigen::MatrixXd NormalVarianceMixtureBaseError::d2Given(const int i ,
                    const Eigen::VectorXd& res,
                    const double weight){
  Eigen::MatrixXd d2 = Eigen::MatrixXd::Zero(npars - 1, npars - 1);

  
  if(assymetric){
    Eigen::VectorXd B_mu =  (-1 + Vs[i].array());
    Eigen::VectorXd iV = Vs[i].cwiseInverse();
    d2(0, 0) =   B_mu.transpose()  * (iV.asDiagonal() * B_mu);  
    d2(0, 0) *= weight/pow(sigma,2);
  }

  return d2;
}


void NormalVarianceMixtureBaseError::get_param(std::vector<double> & param_in)
{
 MeasurementError::get_param(param_in); 
 if(assymetric)
    param_in.push_back(mu);
};
void NormalVarianceMixtureBaseError::get_param_names(Rcpp::StringVector & names){
  MeasurementError::get_param_names(names); 
  if(assymetric)
    names.push_back("mu_error");
};
