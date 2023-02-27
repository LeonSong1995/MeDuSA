//Head file

#ifndef _MLM_H
#define _MLM_H

#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
#include <math.h>


using namespace std;
using namespace Eigen;
using namespace Rcpp;

typedef MatrixXd eigenMatrix;
typedef VectorXd eigenVector;

//REML paramater
vector<eigenMatrix> _A;
eigenMatrix Z; 
eigenMatrix X;
eigenVector y;
eigenMatrix Vi;
eigenMatrix P;
eigenVector varcomp_init;
eigenVector Py;
bool flag_EM = false;
vector<double>L_history;
eigenMatrix mat;


//REML function
vector<eigenMatrix> calcu_A(vector<eigenMatrix> &Z, int n, int rindx, eigenMatrix &S);
bool calcu_Vi(eigenMatrix &Vi, eigenVector &prev_varcmp, double &logdet,int n, int rindx);
bool calcu_P(eigenMatrix &X,eigenMatrix &Vi, eigenMatrix &P, double &logdet_Xt_Vi_X);
void calcu_tr_PA(eigenMatrix &P, eigenVector &tr_PA,int n,int rindx);
bool inverse_H(eigenMatrix &H);
bool ai_reml(eigenVector &y,eigenMatrix &P, eigenVector &Py,eigenVector &prev_varcmp, eigenVector &varcmp, int n, int rindx);
void em_reml(eigenVector &y,eigenMatrix &P, eigenVector &Py, eigenVector &prev_varcmp, eigenVector &varcmp, int n, int rindx);
double y_center_square(eigenVector &y);
void constrain_varcmp(double y_squared,eigenVector &varcmp,int rindx);
VectorXd reml_iteration(Eigen::VectorXd start, eigenMatrix &X,eigenVector &y, vector<eigenMatrix> &Z, eigenVector &varcmp, int n, int rindx, int maxiter, eigenMatrix &S, double &lgL, double y_squared);

#endif
