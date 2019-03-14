#ifndef __SUBSAMPLER__H__
#define __SUBSAMPLER__H__

#include "subsampler.h"
#include "error_check.h"
#include "sample.h"

void subsampler::initFromList(const Rcpp::List & in_list)
{
  pSubsample = Rcpp::as< double > (in_list["pSubsample"]); 
  subsample_type = Rcpp::as< int    > (in_list["subsample_type"]);
  int silent     = Rcpp::as< int    > (in_list["silent"]);
  if(silent == 0){
    Rcpp::Rcout << "Susample type: ";
    if(subsample_type == 4)
      Rcpp::Rcout << "grouped sampler, ";
    else if(subsample_type == 3)
      Rcpp::Rcout << "Poisson sampler, ";
    else if(subsample_type == 1)
      Rcpp::Rcout << "Uniform, ";
    else if(subsample_type == 2)
      Rcpp::Rcout << "Weighted, ";
    Rcpp::Rcout << "procent subsampled: " << pSubsample << "\n";
  }
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
  weight = sampling_weights;
  //weight.array() += nindv / ((double) nSubsample);
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
    groupSampling_weights (pSubsample,
                           groups,
                           free,
                           weight,
                           nSubsample_group);
    if(silent == 0){
      Rcpp::Rcout << "Groups to sample:" <<  nSubsample_group[0] << ", free to sample: "<< nSubsample_group[1] << "\n";
    }
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





#endif
