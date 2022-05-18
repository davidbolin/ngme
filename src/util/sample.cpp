#include "sample.h"
#include <algorithm>    // std::max

/* Unequal probability sampling; without-replacement case */

Eigen::VectorXi ProbSampleNoReplace( int nans, Eigen::VectorXd & p_in)
{
  Eigen::VectorXd p = p_in;
  int n = p.size();
  Eigen::VectorXi ans(nans);
  double rT, mass, totalmass;
  int i, j, k, n1;

  /* Record element identities */
  int *perm;
  //perm  = Calloc(n, int);
  perm = (int*)calloc(n, sizeof(int));
  for (i = 0; i < n; i++)
    perm[i] = i ;

  /* Sort probabilities into descending order */
  /* Order element identities in parallel */
  double *pdat = &p(0);
  revsort(pdat, perm, n);

  /* Compute the sample */
  totalmass = 1;
  for (i = 0, n1 = n-1; i < nans; i++, n1--) {
    rT = totalmass * unif_rand();
    mass = 0;
    for (j = 0; j < n1; j++) {
      mass += p[j];
      if (rT <= mass)
        break;
    }
    ans[i] = perm[j] ;
    totalmass -= p[j];
    for(k = j; k < n1; k++) {
      p[k] = p[k + 1];
      perm[k] = perm[k + 1];
    }
  }
  free(perm);
  return(ans);
}


void ProbSampleNoReplace( int nans, Eigen::VectorXd & p_in, std::vector<int> & ans)
{
  Eigen::VectorXd p = p_in;
  int n = p.size();
  double rT, mass, totalmass;
  int i, j, k, n1;

  /* Record element identities */
  int *perm;
  //perm  = Calloc(n, int);
  perm = (int*)calloc(n, sizeof(int));
  for (i = 0; i < n; i++)
    perm[i] = i;

  /* Sort probabilities into descending order */
  /* Order element identities in parallel */
  double *pdat = &p(0);
  revsort(pdat, perm, n);

  /* Compute the sample */
  totalmass = 1;
  for (i = 0, n1 = n-1; i < nans; i++, n1--) {
    rT = totalmass * unif_rand();
    mass = 0;
    for (j = 0; j < n1; j++) {
      mass += p[j];
      if (rT <= mass)
        break;
    }
    ans[i] = perm[j] ;
    totalmass -= p[j];
    for(k = j; k < n1; k++) {
      p[k] = p[k + 1];
      perm[k] = perm[k + 1];
    }
  }
  free(perm);
}

void ProbSampleNoReplace( int nans,
						  Eigen::VectorXd & p_in,
						  std::vector<int> & ans,
						  std::vector<int> & selected)
{
  Eigen::VectorXd p = p_in;
  int n = p.size();
  double rT, mass, totalmass;
  int i, j, k, n1;

  /* Record element identities */
  int *perm;
  //perm  = Calloc(n, int);
  perm = (int*)calloc(n, sizeof(int));
  for (i = 0; i < n; i++)
    perm[i] = i ;

  /* Sort probabilities into descending order */
  /* Order element identities in parallel */
  double *pdat = &p(0);
  revsort(pdat, perm, n);

  /* Compute the sample */
  totalmass = 1;
  for (i = 0, n1 = n-1; i < nans; i++, n1--) {
    rT = totalmass * unif_rand();
    mass = 0;
    for (j = 0; j < n1; j++) {
      mass += p[j];
      if (rT <= mass)
        break;
    }
    ans[i] = perm[j] ;
    selected[perm[j] ] = 1;
    totalmass -= p[j];
    for(k = j; k < n1; k++) {
      p[k] = p[k + 1];
      perm[k] = perm[k + 1];
    }
  }
  free(perm);
}

int poissonSampling_internal( int nans,
					  Eigen::VectorXd & p_in,
					  Eigen::VectorXd & weight,
					  std::vector<int> & ans,
					  std::vector<int > & selected
					  )
{
  Eigen::VectorXd p = p_in;
  int n = p.size();
  p.array() /= p.sum();
  int i;
  double U;
  /* Record element identities */
  int *perm;
  //perm  = Calloc(n, int);
  perm = (int*)calloc(n, sizeof(int));
  for (i = 0; i < n; i++)
    perm[i] = i ;

  /* Sort probabilities into descending order */
  /* Order element identities in parallel */
  double *pdat = &p(0);
  revsort(pdat, perm, n);
  int counter = 0;
  double rem = 0; // reminder

  for (i = 0; i < n; i++) {
    U = unif_rand();
    int j = perm[i];
    double p_temp =  p[i] * (nans + rem);
    if( U < p_temp)
    {
    	if( selected[j ] == 0){
    		ans.push_back(j );
    		counter++;
    	}
    }
    if( p_temp >= 1){ // ensure that the expected number actually is m
    		weight[j ]  = 1.;
    		rem += p_temp - 1.;
    }else{
    	if( weight[j] > 0)
    		weight[j] = 1./( 1. - (1. - 1./weight[j]) * (1 - p_temp));
    	else
    		weight[j] =  1./p_temp;
    }

  }
  /*
  if(counter < 0.5* nans)
  {
  	double p_sum = 0;
  	double p_max = 0;
  	double p_min = 1;
  	for(int i = 0; i < n; i++){
  		p_sum += p[i];
  		if(p[i]> p_max){p_max = p[i];}
  		if(p[i]< p_min){p_min = p[i];}
  	}
  		for(int i = 0; i < n; i++)
  			Rcpp::Rcout << "p[" << i << "] = " << p[i] << "\n";
  	Rcpp::Rcout << "rem = " << rem << "\n";
  	Rcpp::Rcout << "counter2 = " << counter2 << "\n";
  	Rcpp::Rcout << "counter = " << counter << "\n";
  		Rcpp::Rcout << "max,min, sum  = " << p_max << "," << p_min << "," << p_sum << "\n";
  }*/
  free(perm);
  return(counter);
}


void groupSampling_weights (double pSubsample,
                            std::vector<Eigen::VectorXd> & groups,
                            Eigen::VectorXd & free,
                            Eigen::VectorXd & weight,
                            int * nSubsample_group)
{
  
  int n_indv_in_group   = 0;
  int ngroup       = groups.size(); //number of groups
  int nfree        = free.size(); //number of free subjects
  double n_average;
  for (int i=0; i< ngroup; i++)
    n_indv_in_group += groups[i].size();
  int nindv = n_indv_in_group + nfree; //total number of subjects
  int nSubsample = ceil(pSubsample * nindv); //number of subject to subsample
  n_average = n_indv_in_group / ((double) ngroup); //average number of subjects per group
  nSubsample_group[0] = ceil(pSubsample*ngroup); 
  if(nfree == 0){
    nSubsample_group[1] = 0;
  } else {
    double tmp = (double) pSubsample*nindv;
    tmp -= (double) nSubsample_group[0]*n_average;
    int temp = round(tmp);
    temp = std::max(temp,1);
    nSubsample_group[1] = std::min(temp, nfree); //number of free subjects to sample: min(max(round(p*n - gs*ga),1),nf)
  }
  

  for (int   i  = 0; i < ngroup; i++) {
    for (int ii = 0; ii < groups[i].size(); ii++)
      weight[groups[i][ii]] = ngroup / ((double) nSubsample_group[0]);
  }
  for(int i = 0; i < free.size(); i++)
    weight[free[i]] = free.size() / ( (double) nSubsample_group[1]);
}

void groupSampling_sampling(int * nSubsample_group,
                            std::vector<Eigen::VectorXd> & groups,
                            Eigen::VectorXd & free,
                            std::vector<int> & ans,
                            std::default_random_engine & sampler)
{
  int ngroup = groups.size();
  int nfree = free.size();
  int k = 0;
  //Sample group
  std::vector<int> groupInd;
    for (int i=0; i< ngroup; i++) groupInd.push_back(i);
  if(ngroup>0){
    std::shuffle(groupInd.begin(), groupInd.end(), sampler);
    for(int i = 0; i < nSubsample_group[0]; i++){
      k = groups[groupInd[i]].size();
      for (int ii = 0; ii < k; ii++)
        ans.push_back(groups[groupInd[i]][ii]);

    }
  }

  //Sample remaining free elements

  if(nSubsample_group[1]>0){
    std::vector<int> freeInd;
    for (int i=0; i< nfree; i++)
      freeInd.push_back(i);

    std::shuffle(freeInd.begin(), freeInd.end(), sampler);
    for (int i=0; i< nSubsample_group[1]; i++)
      ans.push_back(free[freeInd[i]]);
  }

}

