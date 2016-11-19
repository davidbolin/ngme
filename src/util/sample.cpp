#include "sample.h"


/* Unequal probability sampling; without-replacement case */

Eigen::VectorXi ProbSampleNoReplace( int nans, Eigen::VectorXd p)
{                                   
  int n = p.size();
  Eigen::VectorXi ans(nans);  
  double rT, mass, totalmass;
  int i, j, k, n1;
  
  /* Record element identities */
  int *perm;
  perm  = Calloc(n, int);
  for (i = 0; i < n; i++)
    perm[i] = i + 1;
  
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
    ans[i] = perm[j] - 1;
    totalmass -= p[j];
    for(k = j; k < n1; k++) {
      p[k] = p[k + 1];
      perm[k] = perm[k + 1];
    }
  }
  free(perm);
  return(ans);
}

