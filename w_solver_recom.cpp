#include <cstdio>
#include <cstring>
#include <iostream>
#include <stdlib.h>  
#include <vector>
#include <algorithm>
#include <math.h> 
#include <set>
#include <omp.h>
#include <string.h>
#include <unistd.h>
#include "mex.h"
#include <matrix.h>
#include <time.h>
#include "tron.h"
using namespace std;
void usage()
{
	mexErrMsgTxt("Usage:[W] = w_solver_recom(R_train,Zc,Zc^t,lambda,eps,maxiter)");
	
}
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]){

	if(nrhs!=6)
		usage();
	int N = mxGetM(prhs[0]);
	int D = mxGetN(prhs[0]);
	double* R = mxGetPr(prhs[0]);
	size_t* R_ir =  mxGetIr(prhs[0]);
	size_t* R_jc =  mxGetJc(prhs[0]);
	double* Zc = mxGetPr(prhs[1]);
	double* Zct = mxGetPr(prhs[2]);
	int K = mxGetN(prhs[1]);
	double lambda = mxGetScalar(prhs[3]);
	double eps = mxGetScalar(prhs[4]);
	double maxiter = mxGetScalar(prhs[5]);
	plhs[0]=mxCreateDoubleMatrix(K,D,mxREAL);	
	double* W_opt = mxGetPr(plhs[0]);
	//#pragma omp parallel for
    for(int i=0;i<D;i++){
        // parallezed appended ? 
    	function* Ab = new function(R,R_ir,R_jc,Zc,Zct,lambda,N,K,D,i);
    	TRON* solver = new TRON(Ab,eps,maxiter); 
    	solver->tron_quad(W_opt+K*i);

    	delete Ab;
	    delete solver;

    }
}