#include "../include/utils.h"
#define MIN(x,y) (x<y)?x:y

//numerically evaluate a jacobian. this creates a dense jacobian, so may not be
//suitable for very large functions
cs* CheckDerivativeN(int m, int n, Fn function, N_Vector curr_guess, void *user_data){

	N_Vector fxh=N_VNew_Serial(m);
	N_Vector fx=N_VNew_Serial(m);
	double h=1e-6, drop_tol=1e-6, *u0=NV_DATA_S(curr_guess), *fxh0=NV_DATA_S(fxh), *fx0=NV_DATA_S(fx), df, *Jx;
	int j, i, nz=0, *Jp, *Ji;
	
	//allocate the jacobian matrix and check memory
	cs *jac=cs_spalloc(m, n, 20*n, 1, 0);
	if(!jac) return NULL;
	Jp=jac->p; Ji=jac->i; Jx=jac->x; Jp[0]=0;
	
	//for derivatives, we need only evaluate f(x) in (f(x+h)-f(x))/h
	//once
	function(curr_guess, fx, user_data);
	
	//build the jacobian one column at a time
	for(j=0;j<n;j++){
		u0[j]+=h;
		function(curr_guess, fxh, user_data);
		for(i=0;i<m;i++){
			df=(fxh0[i]-fx0[i])/h;
			if(abs(df)>drop_tol){
				Jx[nz]=df;
				Ji[nz++]=i;
			}
		}
		Jp[j+1]=nz;
		u0[j]-=h;
	}
	//free up unused jacobian memory
	cs_sprealloc(jac, nz);
	N_VDestroy_Serial(fxh);
	N_VDestroy_Serial(fx);
	return jac;
}

double Norm2(int n, double *v){
	int i;
	double norm=0.0;
	for(i=0;i<n;i++){
		norm+=sqrt(v[i]*v[i]);
	}
	return norm;
}

void PrintVec(int n, double *v){
	int i;
	for(i=0;i<n;i++)printf("%lg\n", v[i]);
	printf("\n");
	return;
}