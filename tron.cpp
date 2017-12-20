#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include "tron.h"
#include <unistd.h>
#include <algorithm> 
#ifndef min
template <class T> static inline T min(T x,T y) { return (x<y)?x:y; }
#endif

#ifndef max
template <class T> static inline T max(T x,T y) { return (x>y)?x:y; }
#endif

#ifdef __cplusplus
extern "C" {
#endif

extern double dnrm2_(int *, double *, int *);
extern double ddot_(int *, double *, int *, double *, int *);
extern int daxpy_(int *, double *, double *, int *, double *, int *);
extern int dscal_(int *, double *, double *, int *);

#ifdef __cplusplus
}
#endif

static void default_print(const char *buf)
{
	fputs(buf,stdout);
	fflush(stdout);
}
function::function(double* _R,size_t* _R_ir,size_t* _R_jc,double* _Zc,double* _Zct,double _lambda,int _N,int _K,int _D,int i){
		R=_R;
		R_ir=_R_ir;
		R_jc=_R_jc;
		Zc = _Zc;
		Zct = _Zct;
		lambda = _lambda;
		N = _N;
		K = _K;
		D = _D;
		target = i;
	}
void function::grad(double *w, double *g){
		
		
		for(int i=0;i<K;i++){
			g[i]=0.0;
		}

	//calculate -b  neglect w 
	
		for(int j=R_jc[target];j<R_jc[target+1];j++){
			for(int i=0;i<K;i++){
				g[i] = g[i] - Zct[R_ir[j]*K+i]*R[j];					
			}
		}
	
}
void function::Hv(double *s, double *Hs){
		//sclar product lambda O(K)  lambda*I*w
		double* c = new double[K];
		int incx=1;
	    memcpy(c,s,sizeof(double)*K);
		dscal_(&K,&lambda,c,&incx);
		//Matrix-vector product Zr O(\omega_d*K)
		int nnz = R_jc[target+1]-R_jc[target];
		double* o = new double[nnz];
		for(int j=R_jc[target];j<R_jc[target+1];j++)
			o[j-R_jc[target]] = ddot_(&K,Zct+K*R_ir[j],&incx,s,&incx);		

        //Matrix-vector product Zr^t O(\omega_d*K)

        for(int i=R_jc[target];i<R_jc[target+1];i++){
        	for(int j=0;j<K;j++){
        		c[j] = c[j]+ Zct[K*R_ir[i]+j]*o[i-R_jc[target]];
        	}
        }	
		// for(int i=0;i<K;i++){
		// 	Hs[i]=c[i];
		// }
		 memcpy(Hs,c,sizeof(double)*K);
		delete[] c; 
		delete[] o;

	} 





void TRON::info(const char *fmt,...)
{
	char buf[BUFSIZ];
	va_list ap;
	va_start(ap,fmt);
	vsprintf(buf,fmt,ap);
	va_end(ap);
	(*tron_print_string)(buf);
}

TRON::TRON(const function *fun_obj, double eps, int max_iter)
{
	this->fun_obj=const_cast<function *>(fun_obj);
	this->eps=eps;
	this->max_iter=max_iter;
	tron_print_string = default_print;
}

TRON::~TRON()
{
}


void TRON::tron_quad(double *w)
{
	int n = fun_obj->get_nr_variable();
	int i, cg_iter;
	double one=1.0; 
	int inc = 1;
	double *s = new double[n];
	double *r = new double[n];
	double *g = new double[n];
	

	//double f = fun_obj->fun(w);	
	
	//w is X0,g is b(residual),
	for (i=0; i<n; i++)
		w[i] = 0;
	fun_obj->grad(w, g);
	// for (i=0; i<n; i++)
	// 	printf("%f\n",g[i]);
	//sleep(2);

	double cgtol = eps*dnrm2_(&n, g, &inc);

	//double delta = dnrm2_(&n, g, &inc);
	//double delta = 1e300;
	cg_iter = trcg(g, s, r, cgtol);
	daxpy_(&n, &one, s, &inc, w, &inc);
	
	//double fnew = fun_obj->fun(w);
	
	//info("f %5.3e, f_new %5.3e, CG_iter %3d\n", f, fnew, cg_iter);
	
	delete[] g;
	delete[] r;
	delete[] s;
}


int TRON::trcg(double *g, double *s, double *r, double cgtol)
{
	int i, inc = 1;
	int n = fun_obj->get_nr_variable();
	double one = 1.0;
	double *d = new double[n];
	double *Hd = new double[n];
	double rTr, rnewTrnew, alpha, beta;
   //  r: residual  d: conjugate vector H: A  alpha: stepsize  beta: projection
	// s: total iteration update
	for (i=0; i<n; i++)
	{
		s[i] = 0.0;
		r[i] = -g[i];
		d[i] = r[i];
	}
	int cg_iter = 0;
	rTr = ddot_(&n, r, &inc, r, &inc);
	while (cg_iter < max_iter)
	{
		double res = dnrm2_(&n, r, &inc);
		//info("target=%d cg_iter=%d  residual=%f, cgtol=%f\n",fun_obj->get_target(),cg_iter,res,cgtol);
		if (res <= cgtol)
			break;
		//sleep(2);
		cg_iter++;
		//info("Hv enter\n");
		fun_obj->Hv(d, Hd);
		//info("Hv out\n");
		//info("Hd[0]=%f\n",Hd[0]);
		alpha = rTr/ddot_(&n, d, &inc, Hd, &inc);
		daxpy_(&n, &alpha, d, &inc, s, &inc);
		/*if (dnrm2_(&n, s, &inc) > delta)
		{
			info("cg reaches trust region boundary\n");
			alpha = -alpha;
			daxpy_(&n, &alpha, d, &inc, s, &inc);

			double std = ddot_(&n, s, &inc, d, &inc);
			double sts = ddot_(&n, s, &inc, s, &inc);
			double dtd = ddot_(&n, d, &inc, d, &inc);
			double dsq = delta*delta;
			double rad = sqrt(std*std + dtd*(dsq-sts));
			if (std >= 0)
				alpha = (dsq - sts)/(std + rad);
			else
				alpha = (rad - std)/dtd;
			daxpy_(&n, &alpha, d, &inc, s, &inc);
			alpha = -alpha;
			daxpy_(&n, &alpha, Hd, &inc, r, &inc);
			break;
		}*/
		alpha = -alpha;
		daxpy_(&n, &alpha, Hd, &inc, r, &inc);
		rnewTrnew = ddot_(&n, r, &inc, r, &inc);
		beta = rnewTrnew/rTr;
		dscal_(&n, &beta, d, &inc);
		daxpy_(&n, &one, r, &inc, d, &inc);
		rTr = rnewTrnew;
	}

	delete[] d;
	delete[] Hd;

	return(cg_iter);
}

double TRON::norm_inf(int n, double *x)
{
	double dmax = fabs(x[0]);
	for (int i=1; i<n; i++)
		if (fabs(x[i]) >= dmax)
			dmax = fabs(x[i]);
	return(dmax);
}

void TRON::set_print_string(void (*print_string) (const char *buf))
{
	tron_print_string = print_string;
}


