#ifndef LDMOD_SAMPLE_H_
#define LDMOD_SAMPLE_H_

#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector>
#include <random>
// sampling without repleacement
// nans - number of samples
// p    - (n x 1) probabilility of sampling
Eigen::VectorXi ProbSampleNoReplace( int nans, Eigen::VectorXd & p_in);
// same as above but result stored in out
void ProbSampleNoReplace( int nans, Eigen::VectorXd & p_in, std::vector<int> & ans);
// same as above but internal
void ProbSampleNoReplace( int nans,
						  Eigen::VectorXd & p_in,
						  std::vector<int> & ans,
						  std::vector<int>& selected);

//
/* Internal poisson sampling
* we want to sample nans samples with probability p_in
* We use regular poisson_sampling but we use probabililites nans * p_in
* if the probability is larger > 1 we move it equally among the remaning probabilites
* nans     - (int) number of samples
* p_in     - (n x 1) the probabilites of sampling each observations
* weight   - (n x 1) the actual probabililites of being selected
* ans      - (vector<int>)
* selected - (n x 1) is the data selected before 1 else 0
*/
int poissonSampling_internal( int nans,
					  Eigen::VectorXd &,
					  Eigen::VectorXd &,
					  std::vector<int> &,
					  std::vector<int> &
					  );

void groupSampling_internal(std::vector<Eigen::VectorXd> &,
                            Eigen::VectorXd &,
                            std::vector<int> &,
                            std::default_random_engine &);
//
/*
	computing the weights from the group sampling scheme
	
	pSubsample - (double ) percentage of subjects to sample
	groups     - (n_group) [i] - index of the indivuals in the group i
	free       - (n_free)  index of the free indivuals
	weight     - (n)       the weights put for all indiv
	int        - (2)       [0] - how many groups to sample
						   [1] - how many to sample from the rest	


*/
void groupSampling_weights(  double ,
							std::vector<Eigen::VectorXd> & ,
                            Eigen::VectorXd & ,
                            Eigen::VectorXd & ,
                            int * );
//
/*
	Sampling from groups and rest
	
	nSubsample_group        - (2)  [0] - how many groups to sample
								   [1] - how many to sample from the rest	

    groups     			    - (n_group) [i] - index of the indivuals in the group i
	free       				- (n_free)  index of the free indivuals
    ans                     - (Vec)  the ouput sample
    sampler                 - (rand enigne) for sampling

*/
void groupSampling_sampling(int * ,
                            std::vector<Eigen::VectorXd> & ,
                            Eigen::VectorXd & ,
                            std::vector<int> & ,
                            std::default_random_engine & );

#endif
