//READEME
//---------------------------------------------------------------------------------------------------
//The R-C++ interface function, it uses the AI-REML to optimze the model of LMM-CAR
//The return includes:
//1. inverse matrix of V (sigma_cell * (ZSZ') + sigma_e * (I))'
//2. the estimation of varcmp (sigma_cell & sigma_e)
//3. the value of Log-likelihood. 
//---------------------------------------------------------------------------------------------------

#include "INV.hpp"
#include "MLM.h"

using namespace std;
using namespace Eigen;


//[[Rcpp::export]]
std::vector<Eigen::MatrixXd>  reml(Eigen::VectorXd start, Eigen::MatrixXd &X, Eigen::VectorXd &y, std::vector<Eigen::MatrixXd> &Z, int maxiter,Eigen::MatrixXd &S)
{
	// Construct the variance component (varcmp)
	int num_random_effects = Z.size() + 1;
	int n = X.rows();
	Eigen::VectorXd varcmp(num_random_effects);
	varcmp.setZero();

	// Perform REML iteration until convergence or maximum number of iterations is reached.
    double lgL = 1e-20;
    double y_squared = y_center_square(y);
	varcmp = reml_iteration(start, X, y, Z, varcmp,n, num_random_effects,maxiter,S,lgL,y_squared);
	
	// Output the results
	eigenMatrix LogL(1,1);
	LogL(0,0) = lgL;
	std::vector<Eigen::MatrixXd> results(3);
	results[0] = Vi;
	results[1] = varcmp;
	results[2] = LogL;

	return(results);
}

// Construct the ZSZ' matrix
// The first one is the rancmp of indiviudal cells of the focal cell type. 
// The second one is the digonal matrix of residuals.
vector<eigenMatrix> calcu_A(vector<eigenMatrix> &Z, int n, int rindx,Eigen::MatrixXd &S)
{
	vector<eigenMatrix> _A;
	_A.resize(rindx);
	for (int i=0; i<(rindx-1);i++)
	{
		// Note: To compute the vraince explained by rancmp, we divide the number of cells. 
		_A[i] = (Z[i] * S * Z[i].transpose())/(Z[i].cols());
	}
	_A[rindx-1] = eigenMatrix::Identity(n, n);
	return _A;
}

// Compute the inverse of V matrix (Vi)
// V = sigma_cell * (ZSZ') + sigma_e * (I)
bool calcu_Vi(eigenMatrix &Vi, eigenVector &prev_varcmp, double &logdet,int n, int rindx)
{

	Vi = eigenMatrix::Zero(n, n);
	for (int i = 0; i < rindx; i++) {
		Vi += (_A[i]) * prev_varcmp[i];
	}
	// use the LLT method to decompose the V matrix
	INVmethod method_try;
	method_try = INV_LLT;
	int rank = 0;
	bool ret = true;
	if (!SquareMatrixInverse(Vi, logdet, rank, method_try))
	{
		//if V is invertible, a small positive value is added to the diagonals. 
		double d_buf = Vi.diagonal().mean() * 1e-4;
		for (int j = 0; j < n; j++) Vi(j, j) += d_buf;
		//return false when the V matrix is stll can not be inversed. 
		if (!SquareMatrixInverse(Vi, logdet, rank, method_try)) {
			std::cout << "Still can't be inversed" << endl;
			ret = false;
		}
	}
	return ret;
}

// Compute the P matrix (temp matrix during reml iteration)
bool calcu_P(eigenMatrix &X,eigenMatrix &Vi, eigenMatrix &P, double &logdet_Xt_Vi_X)
{
	eigenMatrix Vi_X;
	eigenMatrix Xt_Vi_X_i;
	Vi_X = Vi * X;
	Xt_Vi_X_i = X.transpose() * Vi_X;
	int rank = 0;
	INVmethod method_try;
	method_try = INV_LLT;
	if (!SquareMatrixInverse(Xt_Vi_X_i, logdet_Xt_Vi_X, rank, method_try))return false;

	P = Vi - Vi_X * Xt_Vi_X_i * Vi_X.transpose();
	return true;
}


// Compute the trace of P matrix
void calcu_tr_PA(eigenMatrix &P, eigenVector &tr_PA,int n,int rindx)
{
	tr_PA.resize(rindx);
	double *d_bufs = new double[n];
	for (int i = 0; i < rindx; i++)
	{
		double d_buf = 0.0;
		memset(d_bufs, 0, n * sizeof(double));
		#pragma omp parallel for
		for (int k = 0; k < n; k++) {
			for (int l = 0; l < n; l++) d_bufs[k] += P(k, l)*(_A[i])(k, l);
		}
		for (int k = 0; k < n; k++) {
			d_buf += d_bufs[k];
		}
		tr_PA(i) = d_buf;
	}
	delete[] d_bufs;
}

// Inverse the AI (H) matrix
bool inverse_H(eigenMatrix &H)
{
	double d_buf = 0.0;
	INVmethod method_try;
	method_try = INV_LLT;
	int rank = 0;
	if (!SquareMatrixInverse(H, d_buf, rank, method_try)) return false;
	else return true;
}

// AI-REML interation
bool ai_reml(eigenVector &y,eigenMatrix &P, eigenVector &Py,eigenVector &prev_varcmp, eigenVector &varcmp, int n, int rindx,double step)
{
	// calculate R matrix (first-order derivative) and AI matrix (second-order derivative）
	//--------------------------------------------------------------------------
	eigenVector R(rindx);
	eigenMatrix AI(rindx, rindx);

	Py = P * y;
	eigenMatrix APy(n, rindx);
	#pragma omp parallel for
	for (int i = 0; i < rindx; i++) {
		(APy.col(i)) = (_A[i]) * Py;
	}

	#pragma omp parallel for
	for (int i = 0; i < rindx; i++) {
		R(i) = (Py.transpose()*(APy.col(i)))(0, 0);
		eigenVector cvec = P * (APy.col(i));
		AI(i, i) = ((APy.col(i)).transpose() * cvec)(0, 0);
		for (int j = 0; j < i; j++) AI(j, i) = AI(i, j) = ((APy.col(j)).transpose() * cvec)(0, 0);
	}
	AI = 0.5 * AI;
	eigenVector tr_PA;
	calcu_tr_PA(P, tr_PA, n, rindx);
	R = -0.5 * (tr_PA - R);
	//--------------------------------------------------------------------------

	// update the varcmp (sigma_cell and sigma_e)
	if (!inverse_H(AI)) return false;
	eigenVector delta(rindx);
	delta = AI * R;
	varcmp = prev_varcmp + step *delta;
	return true;
}

// EM-REML interation
// We switch to EM-REML when the AI matrix can not be inversed (rare case). 
void em_reml(eigenVector &y,eigenMatrix &P, eigenVector &Py, eigenVector &prev_varcmp, eigenVector &varcmp, int n, int rindx)
{
	eigenVector tr_PA;
	calcu_tr_PA(P, tr_PA,n, rindx);
	Py = P * y;

	eigenVector R(rindx);
	#pragma omp parallel for
	for (int i = 0; i < rindx; i++)
	{
		R(i) = (Py.transpose()*(_A[i]) * Py)(0, 0);
		varcmp(i) = prev_varcmp(i) - prev_varcmp(i) * prev_varcmp(i) * (tr_PA(i) - R(i)) / n;
	}
}

// Compute the squared norm of the centered y. 
double y_center_square(eigenVector &y)
{
	double y_squared;
	eigenVector y_centered = y.array() - y.mean();
	y_squared = y_centered.squaredNorm();
	return(y_squared);

}

// Constrain the negative varcmp at a small positive value (1e-6*y_squared, rare case) 
void constrain_varcmp(double y_squared,eigenVector &varcmp,int rindx)
{
	double constr_scale = 1e-6;
	for (int i = 0; i < rindx; i++) {
		if (varcmp[i] < 0) varcmp[i] = constr_scale * y_squared;
	}
}

// REML iteration
VectorXd reml_iteration(Eigen::VectorXd start, eigenMatrix &X,eigenVector &y, vector<eigenMatrix> &Z, eigenVector &varcmp, int n, int rindx,int maxiter, eigenMatrix &S, double &lgL, double y_squared)
{
	//set inital values for iteration
	double logdet = 0.0;
	double logdet_Xt_Vi_X = 0.0;
	double prev_lgL = -1e20;
	double dlogL = 1000.0;
	double step = 0.316;
	int L_size;
	eigenVector prev_varcmp(rindx);
	prev_varcmp = start;

	//construct rancmp (ZSZ')
	_A = calcu_A(Z,n, rindx, S);
	//record the value of Log-likelihood. 
	L_history.push_back(lgL);

	//start iteration.
	for (int iter = 0; iter <= maxiter; iter++)
	{
		if (!calcu_Vi(Vi, prev_varcmp, logdet,n, rindx))
		{
			cout << "ERROR!:V matrix is not positive" << endl;
			return varcmp;
		}
		if (!calcu_P(X,Vi, P, logdet_Xt_Vi_X))
		{
			cout << "ERROR!: the X^t * V^-1 * X matrix is not invertible, please check your scRNA-seq data" << endl;
			return varcmp;
		}

		//defult with AI-REML.
		if (flag_EM) em_reml(y,P, Py, prev_varcmp, varcmp, n, rindx);
		//if AI-REML fails to converge, switch to the EM-REML. 
		else if (!ai_reml(y,P, Py, prev_varcmp, varcmp,n, rindx,step)) flag_EM = true;
		
		//constrain the negative varcmp at a small positive value.
		constrain_varcmp(y_squared,varcmp,rindx);
		//compute the Log-likelihood.
		lgL = -0.5 * (logdet_Xt_Vi_X + logdet + (y.transpose() * Py)(0, 0));

		//reduce the step size if REML fails to converge after 300 iterations (jump out of the oscillation range).
		if (iter > 300){
			L_size = L_history.size();
			for(int i = 0; i < L_size - 1; i++)
			{
				if(lgL - L_history[i]<1e-30 || lgL - L_history[i]< -1e-30)
				{
					step *= 0.8;
					break;
				}
			}
			L_history.push_back(lgL);
		}

		//compute the difference in Log-likelihood. 
		dlogL = lgL - prev_lgL;

		//check for convergence of the REML. 
		if ((varcmp - prev_varcmp).squaredNorm() / varcmp.squaredNorm() < 1e-6 && 
			(fabs(dlogL) < 1e-3 || (fabs(dlogL) < 1e-2 && dlogL < 0)))
		{
			prev_varcmp = varcmp;
			calcu_Vi(Vi, prev_varcmp, logdet,n, rindx);
			return varcmp;
		}

		//REML fails to converge when the maximum number of iterations is reached. 
		if(iter == maxiter)
		{
			calcu_Vi(Vi, prev_varcmp, logdet,n, rindx);
			cout << "Warning: Log-likelihood not converged. Results are not reliable.\nYou can specify the parameter max_iter to allow for more iterations." << endl;

		}

		//update the value. 
		prev_varcmp = varcmp;
		prev_lgL = lgL;
	}

	return varcmp;
}
