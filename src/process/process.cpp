#include "process.h"
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

void GaussianProcess::initFromList(const Rcpp::List & init_list,const Eigen::VectorXd & h_in)
{
  npars = 0;
  h = h_in;
  iV = h.cwiseInverse();
  std::vector<std::string> check_names =  {"X"};
  check_Rcpplist(init_list, check_names, "GaussianProcess::initFromList");
  Rcpp::List X_list = Rcpp::as<Rcpp::List>  (init_list["X"]);
  nindv = X_list.length();
  Xs.resize(nindv);
  Vs.resize(nindv);
  for(int i = 0; i < nindv; i++ ){
    	Xs[i] = Rcpp::as<Eigen::VectorXd>( X_list[i]);
    	Vs[i] = h;
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

void GHProcess::initFromList(const Rcpp::List & init_list,const  Eigen::VectorXd & h_in)
{
  
  h = h_in;
  h2 = h.cwiseProduct(h);
  h_sum = h.sum();
  h_min = h.minCoeff();
  h3_mean = h.array().pow(3).sum()/h.size();
  std::vector<std::string> check_names =  {"X","V"};
  check_Rcpplist(init_list, check_names, "GHProcess::initFromList");
  Rcpp::List V_list           = Rcpp::as<Rcpp::List>  (init_list["V"]);
  Rcpp::List X_list = Rcpp::as<Rcpp::List>  (init_list["X"]);
  nindv = X_list.length();
  Xs.resize(nindv);
  Vs.resize(nindv);
  for(int i = 0; i < nindv; i++ ){

    	Xs[i] = Rcpp::as<Eigen::VectorXd>( X_list[i]);
    	//Vs[i] = h;
      	Vs[i] = Rcpp::as<Eigen::VectorXd>( V_list[i]);
  	}
  	
  	
  	type_process = Rcpp::as<std::string> (init_list["noise"]);
  	
  	nu = 0;
  	mu = 0;
  	npars = 0;
  	if(type_process != "CH"){
		npars = 2;
  		if(init_list.containsElementNamed("nu"))
    		nu = Rcpp::as < double >( init_list["nu"]);
  		else
    		nu = 1.;
    
  		if(init_list.containsElementNamed("mu"))
    		mu = Rcpp::as < double >( init_list["mu"]);
  		else
    		mu = 1.;
    }

  	
  	dmu  = 0;

  	EV = h;
  	if(type_process != "CH")
  		update_nu();
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
  solver.compute(Qi);
  Eigen::VectorXd b = A.transpose()*Y/ sigma2;
  Eigen::VectorXd temp  =  - h;
  temp = temp.cwiseProduct(iV);
  temp.array() += 1.;
  temp *= mu;
  b +=  K.transpose() * temp;
  Xs[i] = solver.rMVN(b, Z);
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
  	Z[ii] += - mu * h[ii] + Vs[i][ii] * mu;
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
  Eigen::VectorXd temp  =  - h;
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
                 h,
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
                 h,
                 nu_in,
                 type_process);

}

void GHProcess::gradient( const int i ,
			   			  const Eigen::SparseMatrix<double,0,int> & K,
			   			  const Eigen::SparseMatrix<double,0,int> & A,
			   			  const Eigen::VectorXd& res,
			   			  const double sigma,
			   			  const double trace_var)
{

	
  	counter++;
  	
  	if( type_process == "CH")
  		return;
  	
	if( type_process == "NIG"){
		gradient_mu_centered(i, K);
	}else if(type_process=="GAL"){
  			iV = Vs[i].cwiseInverse();
		    Eigen::VectorXd temp_1  =  Vs[i];
		    temp_1 -= h;


		    Eigen::SparseLU< Eigen::SparseMatrix<double,0,int> > LU(K);  // performs a LU factorization of K
      		temp_1 = LU.solve(temp_1);         // use the factorization to solve for the given right hand side
      		
      		
      		Eigen::VectorXd temp_2 = A * temp_1;
			
      		 Eigen::VectorXd temp_3 = - A * Xs[i];
	
      		temp_3 += res;
      		dmu    += temp_2.dot(temp_3) / pow(sigma,2);
      		ddmu_1 -= Vv_mean * (trace_var / pow(sigma, 2));
	}
  grad_nu(i);
}

void GHProcess:: gradient_v2( const int i ,
			   			  const Eigen::SparseMatrix<double,0,int> & K,
			   			  const Eigen::SparseMatrix<double,0,int> & A,
			   			  const Eigen::VectorXd& res,
			   			  const double sigma,
			   			  const Eigen::VectorXd& iV_noise,
			   			  const double EiV_noise,
			   			  const double trace_var)
{

  	counter++;

  	if( type_process == "CH")
  		return;
  		
	if( type_process == "NIG"){
		gradient_mu_centered(i, K);
	}else if(type_process=="GAL"){
			iV = Vs[i].cwiseInverse();
		    Eigen::VectorXd temp_1  =  Vs[i];
		    temp_1 -= h;


		    Eigen::SparseLU< Eigen::SparseMatrix<double,0,int> > LU(K);  // performs a LU factorization of K
      		temp_1 = LU.solve(temp_1);         // use the factorization to solve for the given right hand side
      		Eigen::VectorXd temp_2 = A * temp_1;
      		temp_2 *= iV_noise;
      		 Eigen::VectorXd temp_3 = - A * Xs[i];
      		temp_3 += res;
      		dmu    += temp_2.dot(temp_3) / pow(sigma,2);
      		ddmu_1 -= EiV_noise * Vv_mean * (trace_var / pow(sigma, 2));
	}
	grad_nu(i);
}



void GHProcess::gradient_mu_centered(const int i, const Eigen::SparseMatrix<double,0,int> & K)
{
		iV = Vs[i].cwiseInverse();
		Eigen::VectorXd temp_1  =  - h;
  		temp_1 = temp_1.cwiseProduct(iV);

  		temp_1.array() += 1.;

		Eigen::VectorXd temp_2;
		Eigen::VectorXd temp_3 =  Vs[i] ;
  		temp_3 -= h;
		temp_3.array() *= mu;
		temp_2 = K * Xs[i];
		temp_2 -= temp_3;
		dmu    += temp_1.dot(temp_2);
		ddmu_1 += H_mu;

}

void GHProcess::grad_nu(const int i)
{
	iV = Vs[i].cwiseInverse();
	// dnu
	if(type_process == "NIG"){
		dnu  +=  0.5 *( h.size() / nu -   h2.dot(iV) - Vs[i].sum() + 2 * h_sum);
	}else if(type_process == "GAL"){
    	Eigen::VectorXd temp(Vs[i].size());
    	temp.array() = Vs[i].array().log();
		dnu  +=  h_sum * (1. + log(nu)) + h.dot(temp) - Vs[i].sum() - h_digamma;
	}

}


void GHProcess::step_theta(const double stepsize)
{


	
  	if( type_process == "CH"){
  		counter = 0;
  		return;
  	}
  		
	step_mu(stepsize);
	step_nu(stepsize);
	counter = 0;

	if(store_param)
	{
		mu_vec[vec_counter] = mu;
		nu_vec[vec_counter] = nu;
		vec_counter++;
	}
	clear_gradient();
}

void GHProcess::step_mu(const double stepsize)
{
	mu -= (stepsize / ddmu_1 ) * dmu;
	ddmu_1 = 0;
	ddmu_2 = 0;
}

void GHProcess::step_nu(const double stepsize)
{
  double nu_temp = -1;
  if(type_process == "NIG"){
  	ddnu = -  h.size()/ (0.5 * pow(nu,2));
  	ddnu *= counter;
  }else if(type_process == "GAL"){
  	ddnu = h_sum/ nu - h_trigamma;
  	ddnu *= counter;
  }
  dnu /= ddnu;
  double stepsize_temp  = stepsize;
  while(nu_temp < 0)
  {
    nu_temp = nu - stepsize_temp * dnu;
    stepsize_temp *= 0.5;
    if(stepsize_temp <= 1e-16)
        throw("in GHProcess:: can't make nu it positive \n");
  }
  nu = nu_temp;

  update_nu();
}


Eigen::VectorXd GHProcess::get_gradient()
{

	Eigen::VectorXd  g(npars);
	
	
  	if( type_process == "CH")
  		return(g);
	
	g[0] = dmu;
	g[1] = dnu;	
	return(g);
}

void GHProcess::clear_gradient()
{
	dnu = 0;
	dmu = 0;
}

void GHProcess::setupStoreTracj(const int Niter)
{
	
	vec_counter = 0;
	store_param = 1;
	
  	if( type_process == "CH")
  		return;

	mu_vec.resize(Niter);
	nu_vec.resize(Niter);
}


void GHProcess::printIter()
{
	if( type_process != "CH")
		Rcpp::Rcout << "(nu, mu) = " << nu << ", " << mu;
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

  		EiV=  h2.cwiseInverse();
      	EiV *=  1./nu;
  		EiV += h.cwiseInverse();

  		Vv_mean  = h_sum / ( nu * h.size());
  	}else if(type_process == "GAL"){
  		h_digamma  = 0;
  		h_trigamma = 0;

  		for(int i =0; i < h.size() ; i++){
  			double h_nu =  h[i] * nu;
  			h_digamma  += h[i]  * Digamma( h_nu);
  			h_trigamma += h2[i] * Trigamma(h_nu);
  		}
  		Vv_mean  = h_sum / (h.size() * nu);
  		if(h_min * nu > 1.2){
  			EiV = h;
        EiV.array() -= nu;
        EiV = EiV.cwiseInverse();
  		}else{
  			EiV.setOnes(h.size());
  			EiV.array() *= std::numeric_limits<double>::infinity();
  		}


  	}
  	H_mu =  (EV.sum() - EiV.dot(h2));
}

 Eigen::VectorXd  GHProcess::mean_X(const int i){
 
 Eigen::VectorXd mean = -h + Vs[i]; 
 mean.array() *= mu;
 return mean;
 }