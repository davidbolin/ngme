#include "operator_helper.h"

void operator_select(std::string type_operator, operatorMatrix **Kobj)
{
  if(type_operator == "Matern" || type_operator == "matern"){
    *Kobj = new MaternOperator;
  }else if(type_operator == "fd2"){
    *Kobj = new fd2Operator;
  } else {
    Rcpp::Rcout << "Operator type not supported.\n";
  }
}