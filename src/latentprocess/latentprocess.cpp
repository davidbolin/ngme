#include "latentprocess.h"
#include "error_check.h"
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>
#include <Eigen/OrderingMethods>



double Digamma(double x)
{
  return(R::digamma(x));
}

double Trigamma(double x)
{
  return(R::trigamma(x));
}

void GaussianProcess::initFromList(const Rcpp::List & init_list,const std::vector<Eigen::VectorXd >& h_in)
{
  npars = 0;

  //iV = h.cwiseInverse();
  std::vector<std::string> check_names =  {"X"};
  check_Rcpplist(init_list, check_names, "GaussianProcess::initFromList");
  Rcpp::List X_list = Rcpp::as<Rcpp::List>  (init_list["X"]);
  nindv = X_list.length();

  h.resize(nindv);
  for(int i =0; i < nindv; i++){
  	if(h_in.size() > 1)
    	h[i] = h_in[i];
    else
    	h[i] = h_in[0];

  }

  Xs.resize(nindv);
  Vs.resize(nindv);
  for(int i = 0; i < nindv; i++ ){
    	Xs[i] = Rcpp::as<Eigen::VectorXd>( X_list[i]);
    	Vs[i] = h[i];
  	}



}

void GaussianProcess::simulate(const int i,
  		                         Eigen::VectorXd & Z,
			                         const Eigen::SparseMatrix<double,0,int> & A,
                               const Eigen::SparseMatrix<double,0,int> & K,
			                         Eigen::VectorXd& Y,
        cholesky_solver       & solver)
{
	Eigen::SparseMatrix<double,0,int> Q = Eigen::SparseMatrix<double,0,int>(K.transpose());
    Q =  Q * iV.asDiagonal();
    Q =  Q * K;

    Eigen::VectorXd b;
    b.setZero(K.rows());
  	Xs[i] = solver.rMVN(b, Z);
	Y += A * Xs[i];

}

void GHProcess::initFromList(const Rcpp::List & init_list,const  std::vector<Eigen::VectorXd >& h_in)
{
	std::vector<std::string> check_names =  {"X","V"};
  check_Rcpplist(init_list, check_names, "GHProcess::initFromList");
  Rcpp::List V_list           = Rcpp::as<Rcpp::List>  (init_list["V"]);
  Rcpp::List X_list = Rcpp::as<Rcpp::List>  (init_list["X"]);
  nindv = X_list.size();
  Xs.resize(nindv);
  Vs.resize(nindv);
  term1 = 0;
  term2 = 0;


  h.resize(nindv);
  for(int i =0; i < nindv; i++){
  	if(h_in.size() > 1)
    	h[i] = h_in[i];
    else
    	h[i] = h_in[0];

  }



  h2.resize(nindv);
  h_sum.resize(nindv);
  h_min.resize(nindv);

  H_mu.resize(nindv);
  h3_mean.resize(nindv);
  Vv_mean.resize(nindv);
  EiV.resize(nindv);
  EV.resize(nindv);
  h_digamma.resize(nindv);
  h_trigamma.resize(nindv);

  for(int i =0; i < nindv; i++){
    EV[i]      = h[i];
    h2[i]      = h[i].cwiseProduct(h[i]);
    h_sum[i]   = h[i].sum();
    h_min[i]   = h[i].minCoeff();
    h3_mean[i] = h[i].array().pow(3).sum()/h[i].size();
  }
  h_MIN = *std::min_element(h_min.begin(), h_min.end());




  for(int i = 0; i < nindv; i++ ){

    	Xs[i] = Rcpp::as<Eigen::VectorXd>( X_list[i]);
    	//Vs[i] = h;
      	Vs[i] = Rcpp::as<Eigen::VectorXd>( V_list[i]);
  	}


  	type_process = Rcpp::as<std::string> (init_list["noise"]);
if(type_process == "CH"){
	  Vv_mean.resize(h2.size());
		for(int i =0; i < nindv; i++){
			EV[i] = h2[i];
			EV[i].array() *= 0.25 / (0.5+1);
			Vv_mean[i] = h2[i].sum()/h2[i].size();
			Vv_mean[i] *=  0.25 / (0.5+1);
			}
	}
  	nu = 0;
  	mu = 0;
  	dnu_prev = 0;
  	dmu_prev = 0;
  	npars = 1;
  	if(type_process != "CH"){
		npars = 2;
  		if(init_list.containsElementNamed("nu"))
    		nu = Rcpp::as < double >( init_list["nu"]);
  		else
    		nu = 1.;
	}
  	if(init_list.containsElementNamed("mu"))
    	mu = Rcpp::as < double >( init_list["mu"]);
  	else
    	mu = 0.;


  	dmu    = 0;
	ddmu_1 = 0;
  	if(type_process != "CH"){
  		update_nu();
  		clear_gradient();
  		}
  	counter = 0;
  	store_param = 0;

}

void GaussianProcess::sample_X(const int i,
              Eigen::VectorXd & Z,
              const Eigen::VectorXd & Y,
              const Eigen::SparseMatrix<double,0,int> & Q,
              const Eigen::SparseMatrix<double,0,int> & K,
              const Eigen::SparseMatrix<double,0,int> & A,
              const double sigma,
              cholesky_solver       & solver)
{
  double sigma2  =pow(sigma, 2);
  Eigen::SparseMatrix<double,0,int> Qi = Q + (A.transpose()*A)/ sigma2;
  solver.compute(Qi);
  Eigen::VectorXd b = A.transpose()*Y/ sigma2;
  Xs[i] = solver.rMVN(b, Z);
}



void GHProcess::sample_X(const int i,
              Eigen::VectorXd & Z,
              const Eigen::VectorXd & Y,
              const Eigen::SparseMatrix<double,0,int> & Q,
              const Eigen::SparseMatrix<double,0,int> & K,
              const Eigen::SparseMatrix<double,0,int> & A,
              const double sigma,
              cholesky_solver       & solver)
{
  iV = Vs[i].cwiseInverse();
  double sigma2  =pow(sigma, 2);
  Eigen::SparseMatrix<double,0,int> Qi = Q + (A.transpose()*A)/ sigma2;
  Eigen::VectorXd DQ_12 = Qi.diagonal().cwiseInverse().cwiseSqrt();
  Eigen::SparseMatrix<double,0,int> DQD = DQ_12.asDiagonal() * Qi * DQ_12.asDiagonal();
  for (int i = 0; i < Qi.rows(); i++)
    DQD.coeffRef(i,i) += 1e-16;
  solver.compute(DQD);
  Eigen::VectorXd b = A.transpose()*Y / sigma2;
  Eigen::VectorXd temp  =  - h[i];
  temp = temp.cwiseProduct(iV);
  temp.array() += 1.;
  temp.array() *= mu;
  b +=  K.transpose() * temp;
  b = b.cwiseProduct(DQ_12);
  Xs[i] = DQ_12.cwiseProduct(solver.rMVN(b, Z));
}

	//simulate from prior distribution
void GHProcess::simulate(const int i,
			  Eigen::VectorXd & Z,
			  const Eigen::SparseMatrix<double,0,int> & A,
              const Eigen::SparseMatrix<double,0,int> & K,
			  Eigen::VectorXd& Y,
			  cholesky_solver  &  solver)
{

  Z *= Vs[i].cwiseSqrt();
  for(int ii = 0; ii < Z.size(); ii++)
  	Z[ii] += - mu * h[i][ii] + Vs[i][ii] * mu;
      Eigen::SparseLU< Eigen::SparseMatrix<double,0,int> > LU(K);  // performs a Cholesky factorization of A

 Xs[i] = LU.solve(Z);         // use the factorization to solve for the given right hand side


  Y += A * Xs[i];

}


void GaussianProcess::sample_Xv2( const int i,
              Eigen::VectorXd & Z,
              const Eigen::VectorXd & Y,
              const Eigen::SparseMatrix<double,0,int> & Q,
              const Eigen::SparseMatrix<double,0,int> & K,
              const Eigen::SparseMatrix<double,0,int> & A,
              const double sigma,
              cholesky_solver       & solver,
              const Eigen::VectorXd & iV_noise)
{
  double sigma2  =pow(sigma, 2);
  Eigen::SparseMatrix<double,0,int> Qi = Q + (A.transpose()*iV_noise.asDiagonal()*A)/ sigma2;
  solver.compute(Qi);
  Eigen::VectorXd b = A.transpose()* (iV_noise.cwiseProduct(Y) )/ sigma2;
  Xs[i] = solver.rMVN(b, Z);
}


void GHProcess::sample_Xv2(  const int i,
              Eigen::VectorXd & Z,
              const Eigen::VectorXd & Y,
              const Eigen::SparseMatrix<double,0,int> & Q,
              const Eigen::SparseMatrix<double,0,int> & K,
              const Eigen::SparseMatrix<double,0,int> & A,
              const double sigma,
              cholesky_solver       & solver,
              const Eigen::VectorXd & iV_noise )
{
	iV = Vs[i].cwiseInverse();
  double sigma2  =pow(sigma, 2);
  Eigen::SparseMatrix<double,0,int> Qi = Q + (A.transpose()*iV_noise.asDiagonal()*A)/ sigma2;


  solver.compute(Qi);
  Eigen::VectorXd b = A.transpose()* (iV_noise.cwiseProduct(Y) )/ sigma2;
  Eigen::VectorXd temp  =  - h[i];
  temp = temp.cwiseProduct(iV);
  temp.array() += 1.;
  temp *= mu;
  b +=  K.transpose() * temp;
   Xs[i] =solver.rMVN(b, Z);
}

void GHProcess::sample_V(const int i ,
    					              gig & rgig,
                            const Eigen::SparseMatrix<double,0,int> & K)
{
	double nu_in = nu;
	if( type_process == "NIG")
		nu_in = sqrt(nu_in);
 	Vs[i] = sampleV_post(rgig,
                 h[i],
                 K * Xs[i],
                 1.,
                 mu,
                 nu_in,
                 type_process);

}
void GHProcess::simulate_V(const int i ,
    					   gig & rgig)
{
	double nu_in = nu;
	if( type_process == "NIG")
		nu_in = sqrt(nu_in);
 	Vs[i] = sampleV_pre(rgig,
                 h[i],
                 nu_in,
                 type_process);

}

void GHProcess::gradient( const int i ,
			   			  const Eigen::SparseMatrix<double,0,int> & K,
			   			  const Eigen::SparseMatrix<double,0,int> & A,
			   			  const Eigen::VectorXd& res,
			   			  const double sigma,
			   			  const double trace_var,
			   			  const double weight)
{

  	counter++;
  	//if( type_process == "CH")
  	//	return;

	if( type_process == "NIG"){
		gradient_mu_centered(i, K, weight);
	}else if(type_process=="GAL"){
  			iV = Vs[i].cwiseInverse();
		    Eigen::VectorXd temp_1  =  Vs[i];
		    temp_1 -= h[i];


		    Eigen::SparseLU< Eigen::SparseMatrix<double,0,int> > LU(K);  // performs a LU factorization of K
      		temp_1 = LU.solve(temp_1);         // use the factorization to solve for the given right hand side


      		Eigen::VectorXd temp_2 = A * temp_1;

      		Eigen::VectorXd temp_3 = - A * Xs[i];

      		temp_3 += res;
      		dmu    += weight * temp_2.dot(temp_3) / pow(sigma,2);
      		ddmu_1 -= weight  *Vv_mean[i] * (trace_var / pow(sigma, 2));

	}else if( type_process=="CH"){

	}

  if( type_process == "CH")
  		return;


  grad_nu(i, weight);
}

void GHProcess:: gradient_v2( const int i ,
			   			  const Eigen::SparseMatrix<double,0,int> & K,
			   			  const Eigen::SparseMatrix<double,0,int> & A,
			   			  const Eigen::VectorXd& res,
			   			  const double sigma,
			   			  const Eigen::VectorXd& iV_noise,
			   			  const double EiV_noise,
			   			  const double trace_var,
			   			  const double weight)
{

  	counter++;


	if( type_process == "NIG"){
		gradient_mu_centered(i, K, weight);
	}else if(type_process=="GAL" ){
			iV = Vs[i].cwiseInverse();
		    Eigen::VectorXd temp_1  =  Vs[i];
		    temp_1 -= h[i];


		    Eigen::SparseLU< Eigen::SparseMatrix<double,0,int> > LU(K);  // performs a LU factorization of K
      		temp_1 = LU.solve(temp_1);         // use the factorization to solve for the given right hand side
      		Eigen::VectorXd temp_2 = A * temp_1;
      		temp_2.cwiseProduct(iV_noise );
      		 Eigen::VectorXd temp_3 = - A * Xs[i];
      		temp_3 += res;
      		dmu    += weight * temp_2.dot(temp_3) / pow(sigma,2);
      		ddmu_1 -= weight * EiV_noise * Vv_mean[i] * (trace_var / pow(sigma, 2));

	}else if( type_process=="CH"){

	}

  	if( type_process == "CH")
  		return;
	grad_nu(i, weight);
}



void GHProcess::gradient_mu_centered(const int i,
									 const Eigen::SparseMatrix<double,0,int> & K,
									 const double weight)
{
		iV = Vs[i].cwiseInverse();
		Eigen::VectorXd temp_1  =  - h[i];
  		temp_1 = temp_1.cwiseProduct(iV);

  		temp_1.array() += 1.;
		Eigen::VectorXd temp_2;
		Eigen::VectorXd temp_3 =  Vs[i] ;
  		temp_3 -= h[i];
		temp_3.array() *= mu;
		temp_2 = K * Xs[i];
		temp_2 -= temp_3;
		dmu    += weight * temp_1.dot(temp_2);
		ddmu_1 += weight * H_mu[i];
}


void GHProcess::grad_nu(const int i, const double weight)
{
 	iV = Vs[i].cwiseInverse();
	// dnu
	if(type_process == "NIG"){
		dnu  +=  0.5 *( h[i].size() / nu -   h2[i].dot(iV) - Vs[i].sum() + 2 * h_sum[i]);
    ddnu -=   0.5 * h[i].size()/ pow(nu,2);
    term2 += h[i].size();
    term1 += 	h2[i].dot(iV) + Vs[i].sum() - 2 * h_sum[i];
	}else if(type_process == "GAL"){
    	Eigen::VectorXd temp(Vs[i].size());
    	temp.array() = Vs[i].array().log();
		dnu  += weight * ( h_sum[i] * (1. + log(nu)) + h[i].dot(temp) - Vs[i].sum() - h_digamma[i]);
    	ddnu += weight * ( h_sum[i]/ nu - h_trigamma[i] );
	}


}


void GHProcess::step_theta(const double stepsize,
						   const double learning_rate,
						   const double polyak_rate,
						   const int burnin)
{

	if( type_process == "CH")
  		return;
	step_mu(stepsize, learning_rate,burnin);
	counter = 0;
	if(store_param){
		if(vec_counter == 0 || polyak_rate == -1)
			mu_vec[vec_counter] = mu;
		else
			mu_vec[vec_counter] = polyak_rate * mu + (1- polyak_rate) * mu_vec[vec_counter-1];

	}


	step_nu(stepsize, learning_rate,burnin);

	counter = 0;

	clear_gradient();

	if(store_param)
	{
		if(vec_counter == 0 || polyak_rate == -1)
			nu_vec[vec_counter] = nu;
		else
			nu_vec[vec_counter] = polyak_rate * nu + (1- polyak_rate) * nu_vec[vec_counter-1];

		vec_counter++;
	}
}

void GHProcess::step_mu(const double stepsize, const double learning_rate,const int burnin)
{
	dmu /= -  ddmu_1;

  if(std::abs(dmu_prev) > 10*std::abs(dmu))
  	dmu_prev *= 0.1;

  dmu_prev = learning_rate * dmu_prev +  dmu;
  double step  = stepsize *  dmu_prev;

	mu +=  step ;//(sqrt(cache_mu) + 1);
	ddmu_1 = 0;
	ddmu_2 = 0;

}

void GHProcess::step_nu(const double stepsize, const double learning_rate, const int burnin)
{
  double nu_temp = -1;
  dnu /=  ddnu;

  dnu_prev = learning_rate * dnu_prev +    dnu;
  double step  = dnu_prev;

  double stepsize_temp  = stepsize;
  while(nu_temp < 0)
  {
    nu_temp = nu - stepsize_temp * step;
    stepsize_temp *= 0.5;
    if(stepsize_temp <= 1e-16)
        throw("in GHProcess:: can't make nu it positive \n");
  }
  // checking min value
  if(type_process=="GAL")
  {
  	if(nu_temp * h_MIN < 0.1)
  		nu_temp = 0.1/h_MIN;
  		dnu_prev = 0;

  }else if(type_process == "NIG"){
    if(burnin == 1){
      nu_temp = term1/term2;
      if(nu_temp * pow(h_MIN,2) < 5e-06){
        nu_temp =5e-06/pow(h_MIN,2);
        dnu_prev = 0;
      }
    }
  }
  nu = nu_temp;

  update_nu();
}


Eigen::VectorXd GHProcess::get_gradient()
{

	Eigen::VectorXd  g(npars);
	g[0] = dmu;

  	if( type_process == "CH")
  		return(g);


	g[1] = dnu;
	return(g);
}

void GHProcess::clear_gradient()
{
	dnu = 0;
	dmu = 0;
	term1 = 0;
	term2 = 0;
}

void GHProcess::setupStoreTracj(const int Niter)
{

	vec_counter = 0;
	store_param = 1;


	mu_vec.resize(Niter);

  	if( type_process == "CH")
  		return;
	nu_vec.resize(Niter);
}


void GHProcess::printIter()
{
	if( type_process != "CH")
		Rcpp::Rcout << "(nu, mu) = " << nu << ", " << mu;
	else
		Rcpp::Rcout << "mu = " <<  mu << "\n";
}

Rcpp::List GHProcess::toList()
{
  Rcpp::List out;
  out["X"]      = Xs;
  out["V"]      = Vs;
  out["noise"]  = type_process;
  out["Cov_theta"]   = Cov_theta;

  if( type_process == "CH")
	return(out);


  out["nu"]     = nu;
  out["mu"]     = mu;
  if(store_param)
  {
  	out["mu_vec"]     = mu_vec;
  	out["nu_vec"]     = nu_vec;
  	out["nu"]         = nu_vec[nu_vec.size() - 1];
    out["mu"]         = mu_vec[mu_vec.size() - 1];
  }

  return(out);
}

Rcpp::List GaussianProcess::toList()
{
  Rcpp::List out;
  out["X"] = Xs;
  out["V"] = Vs;
  out["noise"]  = "Normal";
  out["Cov_theta"]   = Cov_theta;
  return(out);
}

void GHProcess::update_nu()
{

  	ddnu = 0;

  if(type_process == "NIG")
  	{
      for(int i =0; i < h2.size(); i++){
  		  EiV[i] =  h2[i].cwiseInverse();
      	  EiV[i] *=  1./nu;
  		  EiV[i] += h[i].cwiseInverse();

  		  Vv_mean[i]  = h_sum[i] / ( nu * h[i].size());
      }
  	}else if(type_process == "GAL"){


  		for(int i =0; i < h.size() ; i++){
        h_digamma[i]  = 0;
    	  h_trigamma[i] = 0;
        for(int ii = 0; ii < h[i].size() ; ii++){
  			  double h_nu =  h[i][ii] * nu;
  			  h_digamma[i]  += h[i][ii]  * Digamma( h_nu);
  			  h_trigamma[i] += h2[i][ii] * Trigamma(h_nu);
        };
        Vv_mean[i]  = h_sum[i] / (h[i].size() * nu);
        if(h_min[i] * nu > 1.2){
    		  EiV[i] = h[i];
          EiV[i].array() -= nu;
          EiV[i] = EiV[i].cwiseInverse();
  		  }else{
  			  EiV[i].setOnes(h[i].size());
  			  EiV[i].array() *= std::numeric_limits<double>::infinity();
  		  }
  		}
  	}else{ throw("no processes");}
  	 for(int i =0; i < h2.size(); i++)
  	 	H_mu[i] =  (EV[i].sum() - EiV[i].dot(h2[i]));

}

Eigen::VectorXd  GHProcess::mean_X(const int i){

 Eigen::VectorXd mean = -h[i] + Vs[i];
 mean.array() *= mu;
 return mean;
 }
