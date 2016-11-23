#ifndef LDMOD_SAMPLE_H_
#define LDMOD_SAMPLE_H_

#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector>
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
#endif 