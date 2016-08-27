#include "error_check.h"


void check_Rcpplist(Rcpp::List const & list_in, std::vector<std::string> names,const char* object_name)
{
  std::vector<std::string>::iterator it;
  for ( it = names.begin(); it != names.end(); it++ )
   {
     if(list_in.containsElementNamed( (*it).c_str()) == 0)
     {
       Rcpp::Rcout << object_name << " must contain " << (*it) << "\n";
       throw("missing list error\n");
        //std::vector<std::string> name = Rcpp::as<std::vector<std::string> >(DF["name"]);
        
     }
     
   }
}