#ifndef __OPERATOR_HELP__
#define __OPERATOR_HELP__
#include <Rcpp.h>
#include <string>
#include "operatorMatrix.h"

// starting up operator
void operator_select(std::string type_operator, operatorMatrix **Kobj);

#endif