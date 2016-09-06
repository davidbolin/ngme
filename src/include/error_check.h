#ifndef __ERRRORCHECKER__H__
#define __ERRRORCHECKER__H__


#include <RcppEigen.h>
#include <Rcpp.h>
#include <string.h>
#include <vector>

void check_Rcpplist(Rcpp::List const & , std::vector<std::string> ,const char* );



#endif