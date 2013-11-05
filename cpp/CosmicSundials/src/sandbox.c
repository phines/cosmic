#include <stdio.h>
#include <stdlib.h>
#include <cs.h>
#include "../include/psc.h"
#include "../include/utils.h"

//residual function we want to test
static int ComputeF(N_Vector curr_guess, N_Vector residual, void* user_data);

//main
int main(int argc, char **argv){

	int i, j, n, pv=0;
	char* file_name;
	cs *jac;
	N_Vector u;
	Gen *gens;
	double *u0;
	
	//load ps file, populate ps, and grab some data
	if(argc > 1)
		file_name = argv[1];
	else{
		printf("Please specify PS data file.");
		return 1;
	}
	
	PS ps=GetPS(file_name);
	if(ps==NULL){
		printf("Unable to build PS data.\n");
		return 1;
	}	
	n=(ps->num_buses-1)<<1;
	gens=ps->gens;
	u=N_VNew_Serial(n);
	u0=NV_DATA_S(u);
	
	//create an initial condition
	if(gens[pv].index==ps->swing_bus)pv++;
	for(i=0;i<ps->ybus->n;i++){
		if(i==ps->swing_bus) continue;
		j=(i>ps->swing_bus?(i-1)<<1:i<<1);
		if(gens[pv].index==i){
			u0[j]=0.1;
			u0[j+1]=0.1;
			pv++;
			if(pv<ps->num_gens&&gens[pv].index==ps->swing_bus)pv++;
			continue;
		}
		u0[j]=1.0; u0[j+1]=0.1;
	}

	//numerically evaluate the jacobian at initial condition
	jac=CheckDerivativeN(n, n, ComputeF, u, (void*)ps);
	cs_di_print(jac, 0);
	
	//free memory and return
	N_VDestroy_Serial(u);
	return 0;
}

int ComputeF(N_Vector curr_guess, N_Vector residual, void *user_data){
	
	//function variables
	PS ps=(PS)user_data;
	cs_ci *ybus=ps->ybus;
	int j,jp,p,i,n=ybus->n,*Yp=ybus->p,*Yi=ybus->i,swing_bus=ps->swing_bus, pv=0,sh=0,num_gens=ps->num_gens,num_shunts=ps->num_shunts;
	double *v=NV_DATA_S(curr_guess),*r=NV_DATA_S(residual),vm,vr,vg,vk2,P,Q,g,b;
	double _Complex *Yx=ybus->x;
	Shunt *shunts=ps->shunts;
	Gen *gens=ps->gens;
	
	//increment shunt and gen pointers if their first buses are the swing bus
	if(shunts[sh].index==swing_bus)sh++;
	if(gens[pv].index==swing_bus)pv++;
	for(j=0;j<((n-1)<<1);j++) r[j]=0.0;
	
	//loop over the columns of ybus to calculate I=Ybus*V+I_l
	for(j=0;j<n;j++){
	
		//skip swing bus
		if(j==swing_bus) continue;	
		jp=(j>swing_bus?(j-1):j)<<1;
		
		//add generator currents
		vr=v[jp];vm=v[jp+1];vk2=vr*vr+vm*vm;
		
		if(pv<num_gens&&j==gens[pv].index){
			vg=gens[pv].vsp;vm=vr;vk2=vg*vg;vr=sqrt(vk2-vm*vm);
			P=gens[pv].pg;
			Q=gens[pv].qg;
			r[jp]-=(P*vm-Q*vr)/vk2;r[jp+1]-=(P*vr+Q*vm)/vk2;
			pv++;
			if(gens[pv].index==swing_bus&&pv<num_gens)pv++;
		}
		//add load currents
		if(sh<num_shunts&&j==shunts[sh].index){
			P=shunts[sh].p;
			Q=shunts[sh].q;
			r[jp]+=(P*vm-Q*vr)/vk2;r[jp+1]+=(P*vr+Q*vm)/vk2;
			sh++;
			if(sh<num_shunts&&shunts[sh].index==swing_bus)sh++;
		}		
		//add ybus currents
		for(p=Yp[j];p<Yp[j+1];p++){
			if(Yi[p]==swing_bus) continue;
			g=creal(Yx[p]);b=cimag(Yx[p]);
			i=(Yi[p]>swing_bus?(Yi[p]-1):Yi[p])<<1;
			r[i]+=(b*vr+g*vm);
			r[i+1]+=(g*vr-b*vm);
		}
	}	
	//correct pv buses: P=VrIr+VmIm.
	for(pv=0;pv<num_gens;pv++){
		if(gens[pv].index==swing_bus)continue;		jp=(gens[pv].index>swing_bus?((gens[pv].index-1)<<1):(gens[pv].index<<1));
		vm=v[jp];vg=gens[pv].vsp;vk2=vg*vg;vr=sqrt(vk2-vm*vm);
		P=vr*r[jp]+vm*r[jp+1];
		r[jp]=vm*P/vk2;r[jp+1]=vr*P/vk2;
	}
	//return
	return 0;
}