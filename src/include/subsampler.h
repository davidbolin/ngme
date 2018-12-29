#ifndef __SUBSAMP__
#define __SUBSAMP__
#include <vector>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <random>
#include <chrono>
#include "MatrixAlgebra.h"

//the subsampler has four different subsample types:
// 1. Uniform sampling without replacement from 1:nrep
// 2. Sapling weighted by number of samples for each patient
// 3. Poisson sampler
// 4. The grouped subsampler
class subsampler {
  protected:
    
    public:
      
      double pSubsample; //percentage to subsample each iteration
      std::vector<int> longInd;
      int subsample_type;
      int nindv;
      int nSubsample, nSubsample_i;
      int ngroup;
      double pSubsample2 = 0;
      int nSubsample_group[2];
      std::vector<double> Ysize;
      Eigen::VectorXd free;
      Rcpp::List group_list;
      Eigen::VectorXd weight, p_inv, p, p_N;
      std::vector<int> selected;
      std::vector<Eigen::VectorXd> groups;
      
      void initFromList(const Rcpp::List &);
      void sample(int, std::default_random_engine &);
};


#endif
