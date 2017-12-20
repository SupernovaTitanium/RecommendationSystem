#ifndef _TRON_H
#define _TRON_H

class function
{
public:
	//virtual double fun(double *w) = 0 ;
	//f(x)=1/2x^tHx+g^tx+b  
    function(double* _R,size_t* _R_ir,size_t* R_jc,double* _Zc,double* _Zct,double _lambda,int _N,int _K,int _D,int i);

	virtual void grad(double *w, double *g);

	//
	virtual void Hv(double *s, double *Hs);
	
	virtual int get_nr_variable(void){
		return K;		
	} 
	virtual int get_target(void){
		return target;		
	} 
	virtual ~function(void){}
private:
		double* R;
		size_t* R_ir;
		size_t* R_jc;
		double* Zc;
		double* Zct;
		double lambda;
		int N;
		int K;
		int D;
		int target;	
};

class TRON
{
public:
	TRON(const function *fun_obj, double eps = 0.1, int max_iter = 100);
	~TRON();

	void tron_quad(double* w);



	void set_print_string(void (*i_print) (const char *buf));
	void set_precision(double _eps){
		eps = _eps;
	}

private:
	int trcg(double *g, double *s, double *r, double cgtol);
	double norm_inf(int n, double *x);

	double eps;
	int max_iter;
	function *fun_obj;
	void info(const char *fmt,...);
	void (*tron_print_string)(const char *buf);
};
#endif
