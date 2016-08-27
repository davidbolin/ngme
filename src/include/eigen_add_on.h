


#ifndef eigen_lanng_add_om_h
#define eigen_lanng_add_om_h

#define eigen_matrix_from_list(outname, listname, varname)   Eigen::MatrixXd langaddontemp_##outname = Rcpp::as< Eigen::MatrixXd >(listname[ varname ]); \
                                                                      outname.resize(langaddontemp_##outname.rows(),langaddontemp_##outname.cols());\
                                                                      outname = langaddontemp_##outname ;
#define eigen_vector_from_list(outname, listname, varname)   Eigen::VectorXd langaddontemp_##outname = Rcpp::as< Eigen::VectorXd >(listname[ varname ]); \
                                                                      outname.resize(langaddontemp_##outname.size());\
                                                                      outname = langaddontemp_##outname ;

#endif
