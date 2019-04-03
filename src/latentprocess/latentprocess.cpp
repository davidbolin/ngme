#include "latentprocess.h"
#include "error_check.h"
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>
#include <Eigen/OrderingMethods>

GaussianProcess::~GaussianProcess(){
  for(int i =0; i < nindv; i++){
    h[i].resize(0);
    Xs[i].resize(0);
    Ws[i].resize(0);
    Vs[i].resize(0);
  }
  h.resize(0);
  Xs.resize(0);
  Ws.resize(0);
  Vs.resize(0);
}

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
  type_process = "Normal";
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
  Ws.resize(nindv);
  Vs.resize(nindv);
  mu0.resize(nindv);
  for(int i = 0; i < nindv; i++ ){
    	Xs[i]  = Rcpp::as<Eigen::VectorXd>( X_list[i]);
      Ws[i]  = Xs[i];
    	Vs[i]  = h[i];
      mu0[i].setZero(h[i].size());
  	}



}

void GaussianProcess::simulate(const int i,
  		                         Eigen::VectorXd & Z,
			                         const Eigen::SparseMatrix<double,0,int> & A,
                               const Eigen::SparseMatrix<double,0,int> & K,
			                         Eigen::VectorXd& Y,
                               solver& solver)
{
  Z = Z.cwiseProduct(h[i].cwiseSqrt());
  Eigen::VectorXd b;
  b.setZero(K.rows());
  Ws[i] = Z; //Save noise
  Eigen::SparseLU< Eigen::SparseMatrix<double,0,int> > LU(K);  // performs a Cholesky factorization of A
  Xs[i] = LU.solve(Z);         // use the factorization to solve for the given right hand side
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
  Ws.resize(nindv);
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
    	Ws[i] = Xs[i];
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
    useEV = 0;
    if(type_process == "GAL")
      useEV = 0;
    if(init_list.containsElementNamed("useEV"))
      useEV = Rcpp::as < double >( init_list["int"]);

    if(useEV)
    {
      ElogV_post.resize(nindv);
      EV_post.resize(nindv);
    }


    mu0.resize(nindv);
    for(int i =0; i < nindv; i++)
      mu0[i] = mu * (-h[i] + Vs[i]);

    sampleVbool.resize(nindv);
    for(int i =0; i < nindv; i++)
      sampleVbool[i].setOnes(h[i].size());

}




void MGHProcess::initFromList(const Rcpp::List & init_list,const  std::vector<Eigen::VectorXd >& h_in)
{
  std::vector<std::string> check_names =  {"X","V","Bmu","Bnu"};
  check_Rcpplist(init_list, check_names, "GHProcess::initFromList");
  Rcpp::List V_list           = Rcpp::as<Rcpp::List>  (init_list["V"]);
  Rcpp::List X_list = Rcpp::as<Rcpp::List>  (init_list["X"]);
  nindv = X_list.size();
  Xs.resize(nindv);
  Ws.resize(nindv);
  Vs.resize(nindv);
  toSampleV.resize(nindv);
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
    EiV[i].resize(h[i].size());
  }
  h_MIN = *std::min_element(h_min.begin(), h_min.end());



  Rcpp::List Bmu_in = Rcpp::as<Rcpp::List>  (init_list["Bmu"]);
  Rcpp::List Bnu_in = Rcpp::as<Rcpp::List>  (init_list["Bnu"]);
  Bmu.resize(nindv);
  Bnu.resize(nindv);
  for(int i = 0; i < nindv; i++ ){

      Xs[i]  = Rcpp::as<Eigen::VectorXd>( X_list[i]);
      Ws[i]  = Xs[i];
      Vs[i]  = Rcpp::as<Eigen::VectorXd>( V_list[i]);
      Bnu[i] = Rcpp::as<Eigen::MatrixXd>( Bnu_in[i]);
      Bmu[i] = Rcpp::as<Eigen::MatrixXd>( Bmu_in[i]);
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
    npars = Bnu[0].cols() + Bmu[0].cols();
      if(init_list.containsElementNamed("nu"))
        Nu = Rcpp::as < Eigen::VectorXd >( init_list["nu"]);
      else
        Nu.resize(Bnu[0].cols());
  }
    if(init_list.containsElementNamed("mu"))
      Mu = Rcpp::as < Eigen::VectorXd >( init_list["mu"]);
    else
      Mu.resize(Bmu[0].cols());

    Mu *= 0;
    Nu *= 0;
    dMu.resize(Mu.size());
    dMu *= 0;
    ddMu.resize(Mu.size(), Mu.size());
    ddMu *=0;
    dMu_prev.resize(Mu.size());
    dMu_prev *=0;
    dNu.resize(Nu.size());
    dNu *=0;
    ddNu.resize(Nu.size(), Nu.size());
    dNu_prev.resize(Nu.size());
    dNu_prev *=0;
    ddNu *= 0;
    dmu    = 0;
  ddmu_1 = 0;
    if(type_process != "CH"){
      update_nu();
      clear_gradient();
      }
    counter = 0;
    store_param = 0;
    useEV = 0;
    if(type_process == "GAL")
      useEV = 0;
    if(init_list.containsElementNamed("useEV"))
      useEV = Rcpp::as < double >( init_list["int"]);

    if(useEV)
    {
      ElogV_post.resize(nindv);
      EV_post.resize(nindv);
    }


    mu0.resize(nindv);
    for(int i =0; i < nindv; i++)
      mu0[i] = mu * (-h[i] + Vs[i]);




    sampleVbool.resize(nindv);
    if(init_list.containsElementNamed("toSampleV")){
      Rcpp::List toSampleV_list = Rcpp::as<Rcpp::List>  (init_list["toSampleV"]);
      for(int i = 0; i < nindv; i++)
        sampleVbool[i] = Rcpp::as<Eigen::VectorXd>( toSampleV_list[i]);

    }else{
      for(int i =0; i < nindv; i++)
        sampleVbool[i].setOnes(h[i].size());
    }

}



//Sample X|Y when Y = AX + sigma*E
//X|Y = N(Qpost^-1 AY/sigma^2, Qpost^-1)
// X = K^-1W -> W = K*X
void GaussianProcess::sample_X(const int i,
              Eigen::VectorXd & Z,
              const Eigen::VectorXd & Y,
              const Eigen::SparseMatrix<double,0,int> & Q,
              const Eigen::SparseMatrix<double,0,int> & K,
              const Eigen::SparseMatrix<double,0,int> & A,
              const double sigma,
              solver       & solver)
{
  double sigma2  =pow(sigma, 2);
  Eigen::SparseMatrix<double,0,int> Qi = Q + (A.transpose()*A)/ sigma2;
  solver.compute(Qi);
  Eigen::VectorXd b = A.transpose()*Y/ sigma2;
  Xs[i] = solver.rMVN(b, Z);
  Ws[i] = K*Xs[i];

}

Eigen::VectorXd GHProcessBase::get_mean_prior(const int i, const Eigen::SparseMatrix<double,0,int> & K){

    Eigen::VectorXd temp  =  - h[i];
    temp += Vs[i];
    temp *= mu;
    Eigen::SparseLU< Eigen::SparseMatrix<double,0,int> > LU(K);  // performs a Cholesky factorization of A
    Eigen::VectorXd mu_prior = LU.solve(temp);         // use the factorization to solve for the given right hand side
    return(mu_prior);
};





	//simulate from prior distribution
void GHProcessBase::simulate(const int i,
			  Eigen::VectorXd & Z,
			  const Eigen::SparseMatrix<double,0,int> & A,
              const Eigen::SparseMatrix<double,0,int> & K,
			  Eigen::VectorXd& Y,
			  solver  &  solver) //Solver not used
{
  Z = Z.cwiseProduct(Vs[i].cwiseSqrt());
  Eigen::VectorXd b = mean_X(i);
  for(int ii = 0; ii < Z.size(); ii++)
  Z[ii] += b[ii];
  
  Ws[i] = Z; //Save noise
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
              solver       & solver,
              const Eigen::VectorXd & iV_noise)
{
  double sigma2  =pow(sigma, 2);
  Eigen::SparseMatrix<double,0,int> Qi = Q + (A.transpose()*iV_noise.asDiagonal()*A)/ sigma2;
  solver.compute(Qi);
  Eigen::VectorXd b = A.transpose()* (iV_noise.cwiseProduct(Y) )/ sigma2;
  Xs[i] = solver.rMVN(b, Z);
  Ws[i] = K*Xs[i];
}




void GHProcessBase::sample_V(const int i ,
    					              gig & rgig,
                            const Eigen::SparseMatrix<double,0,int> & K)
{

  Eigen::VectorXd KX = K * Xs[i];
	double nu_in = nu;

	if( type_process == "NIG")
		nu_in = sqrt(nu_in);
 	Vs[i] = sampleV_post(rgig,
                 h[i],
                 KX,
                 1.,
                 mu,
                 nu_in,
                 type_process,
                 sampleVbool[i]);


  mu0[i] = mu * (-h[i] + Vs[i]);
  if(useEV)
  {
    if(type_process != "GAL"){
      Rcpp::Rcout << "useEV only implimented for GAL\n";
      exit(-1);
    }
    Eigen::VectorXd  p, b;
    double sigma = 1.;
    b = (KX +  mu * h[i]) / sigma;
    b = b.array().square();
    double b_adj = 1e-14;
    double a  =  pow(mu / sigma, 2);
    p = h[i];
    p.array() *= nu;
    p.array() -= 0.5;
    a += 2 * nu;
    b.array() += b_adj;
    EV_post[i].resize(Vs[i].size());
    ElogV_post[i].resize(Vs[i].size());
    for(int j = 0; j < KX.size(); j++ )
    {
      EV_post[i][j]    = ElogV_GIG(p[j], a, b[j]);
      ElogV_post[i][j] = ElogV_GIG(p[j], a, b[j]);
    }

  }

}

void MGHProcess::sample_V(const int i ,
                            gig & rgig,
                            const Eigen::SparseMatrix<double,0,int> & K)
{
  Eigen::VectorXd KX = K * Xs[i];
  Eigen::VectorXd  Nu_in = Bnu[i]*Nu;
  Nu_in.array() = Nu_in.array().exp();
  if( type_process == "NIG"){
    for(int j = 0; j < Nu_in.size(); j++)
      Nu_in(j) = sqrt(Nu_in(j));
    }
    Eigen::VectorXd  Bmu_in = Bmu[i]*Mu;
  Vs[i] = sampleV_post(rgig,
                 h[i],
                 KX,
                 1.,
                 Bmu_in,
                 Nu_in,
                 type_process,
                 sampleVbool[i]);

  //mu0[i] = (Bmu[i] * Mu) * (-h[i] + Vs[i]);
  if(useEV)
  {
    if(type_process != "GAL"){
      Rcpp::Rcout << "useEV only implimented for GAL\n";
      exit(-1);
    }
    Eigen::VectorXd  p, b;
    double sigma = 1.;
    Eigen::VectorXd Bmumu = Bmu[i]*Mu;
    b = (KX +  (Bmumu * h[i]) )/ sigma;
    b = b.array().square();
    double b_adj = 1e-14;
    Eigen::VectorXd a(h[i].size());
    for(int j = 0; j < Nu_in.size(); j++)
      a(i)  =  pow(Bmumu(j) / sigma, 2);
    p = h[i];
    p *= Nu_in;
    p.array() -= 0.5;
    a += 2 * Nu_in;
    b.array() += b_adj;
    EV_post[i].resize(Vs[i].size());
    ElogV_post[i].resize(Vs[i].size());
    for(int j = 0; j < KX.size(); j++ )
    {
      EV_post[i](j)    = ElogV_GIG(p(j), a(j), b(j));
      ElogV_post[i](j) = ElogV_GIG(p(j), a(j), b(j));
    }

  }

}

void GHProcessBase::simulate_V(const int i ,
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


Eigen::MatrixXd  GHProcessBase::d2Given( const int i ,
                                     const Eigen::SparseMatrix<double,0,int> & K,
                                     const Eigen::SparseMatrix<double,0,int> & A,
                                     const Eigen::VectorXd& res,
                                     const double sigma,
                                     const double trace_var,
                                     const double weight)
{

  Eigen::MatrixXd d2 = Eigen::MatrixXd::Zero(npars, npars);
  iV = Vs[i].cwiseInverse();
  if(type_process == "NIG")
  {
    Eigen::VectorXd B_mu  =  Vs[i] - h[i];
    d2(0, 0) = weight * B_mu.dot(B_mu.cwiseProduct(iV));
    d2(1, 1) = weight  *  0.5 * h[i].size()/ pow(nu,2);
  }else if(type_process == "GAL")
  {
    Eigen::VectorXd B_mu  =  Vs[i] - h[i];
    Eigen::SparseLU< Eigen::SparseMatrix<double,0,int> > LU(K);  // performs a LU factorization of K
    Eigen::VectorXd KB_mu = LU.solve(B_mu);         // use the factorization to solve for the given right hand side
    d2(0, 0) = weight * KB_mu.dot(KB_mu) / pow(sigma, 2);
    d2(1, 1) = - weight * ( h_sum[i]/ nu - h_trigamma[i] );
  }
  return(d2);
}

void GHProcessBase::gradient( const int i ,
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
      		ddmu_1 -= weight  * temp_2.dot(temp_2) / pow(sigma,2);

	}

  if( type_process == "CH")
  		return;

  grad_nu(i, weight);
}

Eigen::VectorXd GHProcessBase::d2Given_cross(  const int i ,
                              const Eigen::SparseMatrix<double,0,int> & K,
                              const Eigen::SparseMatrix<double,0,int> & A,
                              const Eigen::VectorXd& res,
                              const double sigma,
                              const Eigen::MatrixXd Bf,
                              const Eigen::MatrixXd Br,
                              const double weight)
    {
      Eigen::VectorXd d2 = Eigen::VectorXd::Zero(Bf.cols() + Br.cols() + 1);
      if(type_process == "GAL"){

        Eigen::VectorXd B_mu  =  Vs[i] - h[i];
        Eigen::SparseLU< Eigen::SparseMatrix<double,0,int> > LU(K);  // performs a LU factorization of K
        Eigen::VectorXd KB_mu = LU.solve(B_mu);
        if(Bf.cols() > 0)
          d2.head(Bf.cols()) = weight *  (Bf.transpose() * KB_mu)/ pow(sigma,2);

        if(Br.cols() > 0)
          d2.segment(Bf.cols(), Br.cols()) = weight *  (Br.transpose() *  KB_mu)/ pow(sigma,2);

        d2(Bf.cols() +  Br.cols()) = 2 * weight * B_mu.dot( res) /  pow(sigma,3);
      }
      return(d2);
    }
    Eigen::VectorXd GHProcessBase::d2Given_v2_cross(  const int i ,
                              const Eigen::SparseMatrix<double,0,int> & K,
                              const Eigen::SparseMatrix<double,0,int> & A,
                              const Eigen::VectorXd& res,
                              const double sigma,
                              const Eigen::MatrixXd Bf,
                              const Eigen::MatrixXd Br,
                              const Eigen::VectorXd& iV_noise,
                              const double weight)
    {
      Eigen::VectorXd d2 = Eigen::VectorXd::Zero(Bf.cols() + Br.cols() + 1);
      if(type_process == "GAL"){

        Eigen::VectorXd B_mu  =  Vs[i] - h[i];
        Eigen::SparseLU< Eigen::SparseMatrix<double,0,int> > LU(K);  // performs a LU factorization of K
        Eigen::VectorXd KB_mu = LU.solve(B_mu);
        if(Bf.cols() > 0)
          d2.head(Bf.cols()) = weight *  (Bf.transpose() * KB_mu.cwiseProduct(iV_noise))/ pow(sigma,2);

        if(Br.cols() > 0)
          d2.segment(Bf.cols(), Br.cols()) = weight *  (Br.transpose() *  KB_mu.cwiseProduct(iV_noise))/ pow(sigma,2);

        d2(Bf.cols() +  Br.cols()) = 2 * weight * B_mu.dot( res.cwiseProduct(iV_noise)) /  pow(sigma,3);
      }
      return(d2);

    }

Eigen::MatrixXd  GHProcessBase::d2Given_v2( const int i ,
                                        const Eigen::SparseMatrix<double,0,int> & K,
                                        const Eigen::SparseMatrix<double,0,int> & A,
                                        const Eigen::VectorXd& res,
                                        const double sigma,
                                        const Eigen::VectorXd& iV_noise,
                                        const double EiV_noise,
                                        const double trace_var,
                                        const double weight)
{

  Eigen::MatrixXd d2 = Eigen::MatrixXd::Zero(npars, npars);
  iV = Vs[i].cwiseInverse();
  if(type_process == "NIG")
  {
    Eigen::VectorXd B_mu  =  Vs[i] - h[i];
    d2(0, 0) = weight * B_mu.dot(B_mu.cwiseProduct(iV));
    d2(1, 1) = weight  *  0.5 * h[i].size()/ pow(nu,2);
  }else if(type_process == "GAL")
  {
    Eigen::VectorXd B_mu  =  Vs[i] - h[i];
    Eigen::SparseLU< Eigen::SparseMatrix<double,0,int> > LU(K);  // performs a LU factorization of K
    Eigen::VectorXd KB_mu = LU.solve(B_mu);         // use the factorization to solve for the given right hand side
    d2(0, 0) = weight * KB_mu.dot(KB_mu.cwiseProduct(iV_noise)) / pow(sigma, 2);
    d2(1, 1) = -weight * ( h_sum[i]/ nu - h_trigamma[i] );
  }
  return(d2);
}


void GHProcessBase::gradient_v2( const int i ,
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
          Eigen::VectorXd temp_2iV = temp_2;
      		temp_2iV.cwiseProduct(iV_noise );
      		 Eigen::VectorXd temp_3 = - A * Xs[i];
      		temp_3 += res;
      		dmu    += weight * temp_2iV.dot(temp_3) / pow(sigma,2);
      		ddmu_1 -= weight * temp_2iV.dot(temp_2) / pow(sigma, 2);

	}

  	if( type_process == "CH")
  		return;
	grad_nu(i, weight);
}


void GHProcessBase::gradient_mu_centered(const int i,
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

//  (( (KX) - (-h + V) * (Bmu * mu) )^T diag(iV) * (-h + V)) * Bmu 
void MGHProcess::gradient_mu_centered(const int i,
                   const Eigen::SparseMatrix<double,0,int> & K,
                   const double weight)
{
    iV = Vs[i].cwiseInverse();
    Eigen::VectorXd mX = mean_X(i);
    Eigen::VectorXd temp_1  =   (K * Xs[i] ) - mX;
    temp_1 = temp_1.cwiseProduct(iV);

    temp_1.array() *= (-h[i] + Vs[i]).array();
    Eigen::VectorXd h_  =  iV.cwiseProduct(-h[i] + Vs[i]);
    h_.array() *= (-h[i] + Vs[i]).array();
    dMu    += weight * Bmu[i].transpose() * temp_1 ;
    ddMu   -= weight * Bmu[i].transpose() * h_.asDiagonal() * Bmu[i];
    }


void GHProcessBase::grad_nu(const int i, const double weight)
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
      if(useEV){
          temp.array() = ElogV_post[i].array();
          dnu  += weight * ( h_sum[i] * (1. + log(nu)) + h[i].dot(temp) - EV_post[i].sum() - h_digamma[i]);

      }else{

        temp.array() = Vs[i].array().log();
        dnu  += weight * ( h_sum[i] * (1. + log(nu)) + h[i].dot(temp) - Vs[i].sum() - h_digamma[i]);
        //Rcpp::Rcout << "dnu , Vs.sum() , h.dit(temp) = " << dnu << "," << Vs[i].sum() << "," << h[i].dot(temp) << "\n";
        //Rcpp::Rcout << "nu  = " << nu << "\n";


      }

    	ddnu += weight * ( h_sum[i]/ nu - h_trigamma[i] );
              if(Vs[i].sum() > 1e4){
          Rcpp::Rcout << "ddnu = " << ddnu << "\n";
          Rcpp::Rcout << "Vs[i] = \n" << Vs[i] << "\n";
          Rcpp::Rcout << "X[i] = \n" << Xs[i].transpose() << "\n";
          Rcpp::Rcout << "i   = " << i << "\n";
          }
	}


}

void MGHProcess::grad_nu(const int i, const double weight)
{
  iV = Vs[i].cwiseInverse();
  Eigen::VectorXd Nu_vec = Bnu[i] * Nu;
  Nu_vec.array() = Nu_vec.array().exp(); 
  // dnu
  if(type_process == "NIG"){
    Eigen::VectorXd temp(h[i].size());
    Eigen::VectorXd Etemp(h[i].size());
    for(int j = 0; j < h[i].size(); j++){
      //1/V * (V - h)^2
      temp(j)  =  0.5 - 0.5 * ( h2[i](j) * iV(j)  + Vs[i](j) - 2 * h[i](j)) * Nu_vec(j);
      Etemp(j) =  0.5 * (h2[i](j)  * EiV[i](j) - h[i](j)) * Nu_vec(j);
    }
    //ddNu -=  weight * Bnu[i].transpose() *temp.asDiagonal() * Bnu[i];
    ddNu -=  weight * Bnu[i].transpose() *Etemp.asDiagonal()*  Bnu[i];
    dNu  +=  weight * Bnu[i].transpose() * temp;
    
  }else if(type_process == "GAL"){
      
  }


}


void GHProcessBase::step_theta(const double stepsize,
						   const double learning_rate,
						   const double polyak_rate,
						   const int burnin)
{

	if( type_process == "CH")
  		return;
  
	step_mu(stepsize,
          learning_rate,
          burnin,
          polyak_rate);
	counter = 0;

	step_nu(stepsize,
          learning_rate,
           burnin,
           polyak_rate);

	counter = 0;

	clear_gradient();
  vec_counter++;
}

void GHProcessBase::step_mu(const double stepsize, const double learning_rate,const int burnin, const double polyak_rate)
{
  dmu /= -  ddmu_1;

  if(std::abs(dmu_prev) > 10*std::abs(dmu))
    dmu_prev *= 0.1;

  dmu_prev = learning_rate * dmu_prev +  dmu;
  double step  = stepsize *  dmu_prev;

  mu +=  step ;//(sqrt(cache_mu) + 1);
  ddmu_1 = 0;
  ddmu_2 = 0;
if(store_param){
    if(vec_counter == 0 || polyak_rate == -1)
      mu_vec(vec_counter,0) = mu;
    else
      mu_vec(vec_counter,0) = polyak_rate * mu + (1- polyak_rate) * mu_vec(vec_counter-1,0);
  }
}

void MGHProcess::step_mu(const double stepsize, const double learning_rate,const int burnin, const double polyak_rate)
{
  dMu = -ddMu.ldlt().solve(dMu);
  dMu_prev *= learning_rate;
  dMu_prev +=  dMu;
	Mu +=  stepsize *  dMu_prev ;//(sqrt(cache_mu) + 1);
  dMu  *= 0;
  ddMu *=0;
  if(store_param){
      if(vec_counter == 0 || polyak_rate == -1)
        mu_vec.row(vec_counter) = Mu;
      else
        mu_vec.row(vec_counter) = polyak_rate * Mu + (1- polyak_rate) * mu_vec.row(vec_counter-1);
    }
}




void MGHProcess::step_nu(const double stepsize, const double learning_rate, const int burnin, const double polyak_rate)
{
  dNu = -ddNu.ldlt().solve(dNu);
  dNu_prev = learning_rate * dNu_prev +    dNu;
  Nu += stepsize * dNu_prev;
  update_nu();
  if(store_param)
  {
    if(vec_counter == 0 || polyak_rate == -1)
      nu_vec.row(vec_counter) = Nu;
    else
      nu_vec.row(vec_counter) = polyak_rate * Nu + (1- polyak_rate) * nu_vec.row(vec_counter-1);
    
  }
}


void GHProcessBase::step_nu(const double stepsize, const double learning_rate, const int burnin, const double polyak_rate)
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
  	if(nu_temp * h_MIN < 0.01)
  		nu_temp = 0.1/h_MIN;
  		dnu_prev = 0;

  }else if(type_process == "NIG"){
    if (nu_temp <0.00001)
      nu_temp = 0.00001;
    if(burnin == 1){
      nu_temp = term1/term2;
      if(nu_temp * pow(h_MIN,2) < 5e-06){
        nu_temp =5e-06/pow(h_MIN,2);
        dnu_prev = 0;
      }
    }
  }
  //nu_temp = term1/term2;
  nu = nu_temp;
  update_nu();
  if(store_param)
  {
    if(vec_counter == 0 || polyak_rate == -1)
      nu_vec(vec_counter,0) = nu;
    else
      nu_vec(vec_counter,0) = polyak_rate * nu + (1- polyak_rate) * nu_vec(vec_counter-1,0);
    
  }
}


Eigen::VectorXd GHProcessBase::get_gradient()
{

	Eigen::VectorXd  g(npars);
  if( type_process == "CH")
      return(g);
	g[0] = dmu;
	g[1] = dnu;
	return(g);
}

Eigen::VectorXd MGHProcess::get_gradient()
{

  Eigen::VectorXd  g(npars);
  if( type_process == "CH")
      return(g);
  g.head(Bmu[0].cols()) = dMu;
  g.tail(Bnu[0].cols()) = dNu;
  return(g);
}

void GHProcessBase::clear_gradient()
{
	dnu = 0;
	dmu = 0;
	term1 = 0;
	term2 = 0;
}
void MGHProcess::clear_gradient()
{
  dNu *= 0;
  dMu *= 0;
  ddMu *= 0;
  ddNu *= 0;
}

void GHProcessBase::setupStoreTracj(const int Niter)
{

	vec_counter = 0;
	store_param = 1;


	mu_vec.resize(Niter,1);

  	if( type_process == "CH")
  		return;
	nu_vec.resize(Niter,1);
}

void MGHProcess::setupStoreTracj(const int Niter)
{

  vec_counter = 0;
  store_param = 1;


  mu_vec.resize(Niter,Mu.size());

    if( type_process == "CH")
      return;
  nu_vec.resize(Niter,Nu.size());
}


void GHProcessBase::printIter()
{
	if( type_process != "CH")
		Rcpp::Rcout << "(nu, mu) = " << nu << ", " << mu;
	else
		Rcpp::Rcout << "mu = " <<  mu << "\n";
}

Rcpp::List GHProcessBase::toList()
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
  	out["nu"]         = nu_vec.row(nu_vec.size() - 1);
    out["mu"]         = mu_vec.row(mu_vec.size() - 1);
  }

  return(out);
}
Rcpp::List MGHProcess::toList(){
Rcpp::List out;
  out["X"]      = Xs;
  out["V"]      = Vs;
  out["noise"]  = type_process;
  out["Cov_theta"]   = Cov_theta;

  if( type_process == "CH")
  return(out);


  out["nu"]     = Nu;
  out["mu"]     = Mu;
  out["Bmu"]     = Bmu;
  out["Bnu"]     = Bnu;
  if(store_param)
  {
    out["mu_vec"]     = mu_vec;
    out["nu_vec"]     = nu_vec;
    out["nu"]         = nu_vec.row(nu_vec.rows() - 1);
    out["mu"]         = mu_vec.row(mu_vec.rows() - 1);
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

void GHProcessBase::update_nu()
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
        Vv_mean[i]  = h_sum[i] * nu / (h[i].size() * nu);
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

void MGHProcess::update_nu()
{
    ddNu *= 0;
    if(type_process == "NIG")
    {

      for(int i =0; i < h2.size(); i++){
        Eigen::VectorXd Nu_vec = Bnu[i] * Nu;
        Nu_vec.array() = Nu_vec.array().exp(); 
        for(int ii = 0; ii < h2[i].size(); ii++){
          EiV[i](ii) =  1./h2[i](ii);
          EiV[i](ii) *=  1./Nu_vec(ii);
          EiV[i](ii) += 1/h[i](ii);
        }
      }
    }
}



Eigen::VectorXd  MGHProcess::mean_X(const int i){

 Eigen::VectorXd mean = -h[i] + Vs[i];
 mean.array() *= (Bmu[i] * Mu).array();
 return mean;
 }


Eigen::VectorXd  GHProcessBase::mean_X(const int i){

 Eigen::VectorXd mean = -h[i] + Vs[i];
 mean.array() *= mu;
 return mean;
 }

void GHProcessBase::sample_X(const int i,
              Eigen::VectorXd & Z,
              const Eigen::VectorXd & Y,
              const Eigen::SparseMatrix<double,0,int> & Q,
              const Eigen::SparseMatrix<double,0,int> & K,
              const Eigen::SparseMatrix<double,0,int> & A,
              const double sigma,
              solver       & solver)
{
  Eigen::VectorXd X_prev = Xs[i];
  iV = Vs[i].cwiseInverse();
  double sigma2  =pow(sigma, 2);
  Eigen::SparseMatrix<double,0,int> Qi = Q + (A.transpose()*A)/ sigma2;
  Eigen::VectorXd DQ_12 = Qi.diagonal().cwiseInverse().cwiseSqrt();
  Eigen::SparseMatrix<double,0,int> DQD = DQ_12.asDiagonal() * Qi * DQ_12.asDiagonal();
  for (int j = 0; j < Qi.rows(); j++)
    DQD.coeffRef(j, j) += 1e-12;

  solver.compute(DQD);
  Eigen::VectorXd b = A.transpose()*Y / sigma2;
  // change
  Eigen::VectorXd priorMean = mean_X(i);
  Eigen::VectorXd temp  =  - h[i];
  priorMean = priorMean.cwiseProduct(iV);
  b +=  K.transpose() * priorMean;
  b = b.cwiseProduct(DQ_12);
  Xs[i] = DQ_12.cwiseProduct(solver.rMVN(b, Z));
  Ws[i] = K*Xs[i];
}

void GHProcessBase::sample_Xv2(  const int i,
              Eigen::VectorXd & Z,
              const Eigen::VectorXd & Y,
              const Eigen::SparseMatrix<double,0,int> & Q,
              const Eigen::SparseMatrix<double,0,int> & K,
              const Eigen::SparseMatrix<double,0,int> & A,
              const double sigma,
              solver       & solver,
              const Eigen::VectorXd & iV_noise )
{
  iV = Vs[i].cwiseInverse();
  double sigma2  =pow(sigma, 2);
  Eigen::SparseMatrix<double,0,int> Qi = Q + (A.transpose()*iV_noise.asDiagonal()*A)/ sigma2;


  solver.compute(Qi);
  Eigen::VectorXd b = A.transpose()* (iV_noise.cwiseProduct(Y) )/ sigma2;
  Eigen::VectorXd priorMean = mean_X(i);
  Eigen::VectorXd temp  =  - h[i];
  priorMean = priorMean.cwiseProduct(iV);
  b +=  K.transpose() * priorMean;
   Xs[i] =solver.rMVN(b, Z);
   Ws[i] = K*Xs[i];
}