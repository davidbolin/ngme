#ifndef __SUBSAMPLER__H__
#define __SUBSAMPLER__H__

#include "subsampler.h"
#include "error_check.h"
#include "sample.h"

void subsampler::initFromList(const Rcpp::List & in_list)
{
  pSubsample = Rcpp::as< double > (in_list["pSubsample"]); 
  subsample_type = Rcpp::as< int    > (in_list["subsample_type"]);
  Rcpp::List obs_list  = Rcpp::as<Rcpp::List> (in_list["obs_list"]);
  nindv = obs_list.length(); 
  Eigen::VectorXd sampling_weights(nindv);
  nSubsample = ceil(pSubsample * nindv);
  Ysize.resize(nindv,0);
  int count = 0;
  for( Rcpp::List::iterator it = obs_list.begin(); it != obs_list.end(); ++it ) {
    if(subsample_type == 1){
      sampling_weights[count] = 1.0;
    } else if(subsample_type == 2){
      Rcpp::List obs_tmp = Rcpp::as<Rcpp::List>(*it);
      Eigen::VectorXd Ys = Rcpp::as<Eigen::VectorXd>(obs_tmp["Y"]);
      Ysize[count] = (double) Ys.size();
    } else if (subsample_type == 3){ //Biased sampling
      sampling_weights[count] = 1.0;
    } else if (subsample_type == 0){ //Biased sampling
      sampling_weights[count] = 1.0;
    }
    count++;
  }
  
  sampling_weights /= sampling_weights.sum();
  
  
  // For sampling we have the following:
  /*
  weight   -> the weight each observation has given it is selected
  (note this not 1/p)
  p[i]     -> the probability of selecting indvual [i]
  m        -> expected number of samples we want for weighted sampling
  selected -> 1, 0 variable keeping track of which individuals being selected
  */
  p_inv.setZero(nindv);
  weight.setZero(nindv);
  weight.array() += nindv / ((double) nSubsample);
  p.setZero(nindv);
  p_N.setZero(nindv);
  p_N.array() += 1. / nindv;
  selected.resize(nindv,0);
  
  //  Rcpp::Rcout << "nindv = " << nindv << "\n";
  for (int i=0; i< nindv; i++) longInd.push_back(i);
  
  //Rcpp::Rcout << " longInd:\n";
  //for (int i=0; i< nindv; i++) Rcpp::Rcout << longInd[i] << "\n";
  
  
  count = 0;
  if(subsample_type == 3){
    pSubsample2 = Rcpp::as< double > (in_list["pSubsample"]); 
  } else if(subsample_type == 4){
    group_list = Rcpp::as<Rcpp::List> (in_list["group_list"]);
    ngroup = group_list.length();
    //Rcpp::Rcout << "ngroup = " << ngroup << "\n";
    groups.resize(ngroup);
    free = Rcpp::as<Eigen::VectorXd>(in_list["free_samples"]);
    
    for( Rcpp::List::iterator it = group_list.begin(); it != group_list.end(); ++it ) {
      groups[count] = Rcpp::as<Eigen::VectorXd>(*it);
      
      //Rcpp::Rcout << "group: "<< groups[count] << "\n";
      
      
      count++;
    }
    int gsum = 0;
    for(int i=0;i<groups.size();i++){
      int ngroup = groups[i].size();
      gsum += ngroup;
    }
    double gmean = 0;
    if(groups.size()>0){
      gmean = gsum/groups.size();
    }
    //compute weights and decide how many groups to subsample
    groupSampling_weights (nSubsample,
                           groups,
                           free,
                           weight,
                           nSubsample_group);
  }
  

    

}


void subsampler::sample(int iter, std::default_random_engine & subsample_generator)
{
  nSubsample_i = nSubsample;
  // subsampling
  if(subsample_type == 1){
    std::shuffle(longInd.begin(), longInd.end(), subsample_generator);
    //Rcpp::Rcout << " longInd:\n";
    //for (int i=0; i< nindv; i++) Rcpp::Rcout << longInd[i] << "\n";
    
  } else if(subsample_type == 2){
    std::discrete_distribution<int> distribution(Ysize.begin(), Ysize.end());
    for (int i=0; i<nSubsample_i; ++i) {
      longInd[i] = distribution(subsample_generator);
    }
  }else if(subsample_type == 3){
    //longInd
    longInd.resize(nSubsample, 0);
    std::fill(longInd.begin(), longInd.end(), 0);
    std::fill(selected.begin(), selected.end(), 0);
    int m = int(pSubsample2 * nindv); // expected number of samples from the first part
    weight.setZero(nindv);
    if(iter <= 10){
      ProbSampleNoReplace(m+nSubsample, p_N, longInd, selected);
      weight.array() += nindv / ((double) (nSubsample + m));
      nSubsample_i = nSubsample + m;
    }else{
      ProbSampleNoReplace(nSubsample, p_N, longInd, selected);
      weight.array() += nindv / ((double) nSubsample);
      nSubsample_i = nSubsample;
    }
    p = p_N;
    if(iter > 10){
      nSubsample_i += poissonSampling_internal( m,
                                                p_inv,
                                                weight,
                                                longInd,
                                                selected);
    }
    double w_sum = 0;
    for(int ilong = 0; ilong < nSubsample_i; ilong++ )
      w_sum += weight[longInd[ilong]];
    
  } else if(subsample_type ==4){
    nSubsample_i = nSubsample;
    longInd.clear();
    std::fill(longInd.begin(), longInd.end(), 0);
    
    groupSampling_sampling(nSubsample_group,
                           groups,
                           free,
                           longInd,
                           subsample_generator);
    nSubsample_i = longInd.size();
  }else if(subsample_type == 0){
    nSubsample_i = longInd.size();
  }
}

void subsampler::groupSampling_weights (int nSubsample,
                            std::vector<Eigen::VectorXd> & groups,
                            Eigen::VectorXd & free,
                            Eigen::VectorXd & weight,
                            int * nSubsample_group)
{
  int n_indv_in_group   = 0;
  int ngroup       = groups.size();
  int nfree        = free.size();
  double n_average;
  for (int i=0; i< ngroup; i++)
    n_indv_in_group += groups[i].size();
  
  n_average = n_indv_in_group / ((double) ngroup); //average number of individuals per group
  
  double prop_group = ((double) n_indv_in_group) / ((double) nfree); //proportion of grouped samples
  
  //prop_group*nSubsample = how
  //sample ceil(4*proportion_in_groups)
  int n_sample_group = ceil(4 * (prop_group * nSubsample)/ ((double) n_average));
  n_sample_group = std::min(n_sample_group, ngroup); 
  if(n_indv_in_group < nSubsample)
    n_sample_group = ngroup;
  nSubsample_group[0] = n_sample_group;
  int temp = nSubsample - ceil(n_sample_group * n_average);
  nSubsample_group[1] = std::max(temp, 0);
  
  for (int   i  = 0; i < ngroup; i++) {
    for (int ii = 0; ii < groups[i].size(); ii++)
      weight[groups[i][ii]] = ngroup / ((double) n_sample_group);
  }
  for(int i = 0; i < free.size(); i++)
    weight[free[i]] = free.size() / ( (double) nSubsample_group[1]);
}

void subsampler::groupSampling_sampling(int * nSubsample_group,
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

Eigen::VectorXi subsampler::ProbSampleNoReplace( int nans, Eigen::VectorXd & p_in)
{
  Eigen::VectorXd p = p_in;
  int n = p.size();
  Eigen::VectorXi ans(nans);
  double rT, mass, totalmass;
  int i, j, k, n1;
  
  /* Record element identities */
  int *perm;
  perm  = Calloc(n, int);
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


void subsampler::ProbSampleNoReplace( int nans, Eigen::VectorXd & p_in, std::vector<int> & ans)
{
  Eigen::VectorXd p = p_in;
  int n = p.size();
  double rT, mass, totalmass;
  int i, j, k, n1;
  
  /* Record element identities */
  int *perm;
  perm  = Calloc(n, int);
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

void subsampler::ProbSampleNoReplace( int nans,
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
  perm  = Calloc(n, int);
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

int subsampler::poissonSampling_internal( int nans,
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
  perm  = Calloc(n, int);
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
  free(perm);
  return(counter);
}



#endif
