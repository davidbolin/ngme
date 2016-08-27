#ifndef __GIG__H
#define __GIG__H
#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>
#include <Eigen/OrderingMethods>
#include <string>
#include "rgig.h"
double dlambda_V(const double ,
                 const Eigen::VectorXd &, 
                 const Eigen::VectorXd &,
                 const int);
                 
Eigen::VectorXd sampleV_post(gig &,
                        const Eigen::VectorXd &, 
                        const Eigen::VectorXd &,
                        const double,
                        const double,
                        const double,
                        const int);

                 
Eigen::VectorXd sampleV_post(gig &,
                        const Eigen::VectorXd &, 
                        const Eigen::VectorXd &,
                        const double,
                        const double,
                        const double,
                        const std::string);
Eigen::VectorXd sampleV_pre(gig &,
                            const Eigen::VectorXd &, 
                            const double,
                            const std::string);


#endif 