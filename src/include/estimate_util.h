#ifndef __ESTIMATEUTIL_H
#define __ESTIMATEUTIL_H


#include <Rcpp.h>
#include <RcppEigen.h>
#include <random>
#include <chrono>
#include <vector>
#include <string>
#include <math.h>
#include "operatorMatrix.h"
#include "operator_helper.h"
#include "solver.h"
#include "GIG.h"
#include "MixedEffect.h"
#include "measError.h"
#include "latentprocess.h"
#include "subSampleDiagnostic.h"
#include "subsampler.h"
#include "MatrixAlgebra.h"

double estDigamma(double);

Eigen::VectorXd GibbsSampling(int,
                              Eigen::VectorXd&,
                              Eigen::SparseMatrix<double,0,int>&,
                              int,
                              int,
                              int ,
                              MixedEffect&,
                              operatorMatrix&,
                              MeasurementError&,
                              Process&,
                              int,
                              Eigen::VectorXd&,
                              gig&,
                              std::vector<  solver* >&);

void grad_caculations(int,
                      Eigen::VectorXd&,
                      Eigen::SparseMatrix<double,0,int>&,
                      double,
                      int,
                      MixedEffect&,
                      operatorMatrix&,
                      MeasurementError&,
                      Process&,
                      const int,
                      Eigen::MatrixXd& ,
                      int);

#endif 