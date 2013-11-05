#include "../include/ic_solver_NR.h"

//solve the initial condition using simple newton iterations
int SolveIC(PS ps, double *xy0){

	int n=(ps->num_buses-1)<<1, pv=0, i, j, iters=0;
	klu_common comm;
	klu_symbolic *symb;
	klu_numeric *numeric=NULL;
	klu_defaults(&comm);
	Gen *gens=ps->gens;
	double *v, *r, tol=1e-6;
	cs_di *jac, *jac2;
	
	//allocate memory for vectors
	v=malloc(n*sizeof(double));
	if(!v) return 1;
	r=malloc(n*sizeof(double));
	if(!r){free(v);return 1;}
	
	//initial condition
	if(gens[pv].index==ps->swing_bus)pv++;
	for(i=0;i<ps->ybus->n;i++){
		if(i==ps->swing_bus) continue;
		j=(i>ps->swing_bus?(i-1)<<1:i<<1);
		if(gens[pv].index==i){
			v[j]=0.1;v[j+1]=0.1;pv++;
			if(pv<ps->num_gens&&gens[pv].index==ps->swing_bus)pv++;
			continue;
		}
		v[j]=1.0; v[j+1]=0.1;
	}

	//create the jacobian and prep KLU
	jac=InitJac(ps);
	if(!jac){free(v);free(r);return 1;}
	ComputeJac(v,jac,ps);
	jac2=CheckDerivativeN()
	symb=klu_analyze(n, jac->p, jac->i, &comm);
	
	//main newton iterations loop
	ComputeF(v, r, ps);
	while(Norm2(n,r)>tol&&iters<20){	
		if(!numeric) numeric=klu_factor(jac->p, jac->i, jac->x, symb, &comm);
		else klu_refactor(jac->p, jac->i, jac->x, symb, numeric, &comm);
		klu_solve(symb, numeric, n, 1, r, &comm);
		for(i=0;i<n;i++) v[i]=r[i]+v[i];
		ComputeF(v,r,ps);
		cs_di_print(jac,0);
		ComputeJac(v,jac,ps);
		iters++;
	}
	printf("Newton search finished in %d iterations with residual norm %lg.\n", iters, Norm2(n,r));
	
	//free resources and return
	free(r); free(v); jac=cs_spfree(jac);
	klu_free_symbolic(&symb,&comm); klu_free_numeric(&numeric,&comm);
	return 0;
}

//initialize the Jacobian nonzero pattern
cs_di* InitJac(PS ps){

	//check input
	if(!ps||!ps->ybus||!ps->shunts||!ps->gens) return NULL;

	//set up appropriate function variables
	cs_di *jac;
	cs_ci *ybus = ps->ybus;
	int j,ij,p,pv=0,is_pv=0,col=0,nz=0,swing_bus=ps->swing_bus,n=ybus->n,*Yp=ybus->p,*Yi=ybus->i,*Jp,*Ji,num_gens=ps->num_gens;
	Gen* gens=ps->gens;
	double *Jx;
	double _Complex *Yx=ybus->x;
		
	//allocate jacobian memory
	jac=cs_di_spalloc((n-1)<<1, (n-1)<<1, ybus->nzmax<<2, 1, 0);
	if(!jac) return NULL;
	Jp=jac->p;Ji=jac->i;Jx=jac->x;
	Jp[0]=0;
	
	//loop over ybus columns
	if(gens[pv].index==swing_bus)pv++;
	for(j=0; j<n; (!col)?j++:1){
		
		//skip swing bus columns 
		if(j==swing_bus) continue;
		
		//loop through ybus rows
		for(p=Yp[j]; p<Yp[j+1]; p++){
		
			//skip swing bus rows. if we jumped to the next column, we're done.
			if(Yi[p]==swing_bus){p++;if(p==Yp[j+1]) break;}
			ij=(Yi[p]>swing_bus?((Yi[p]-1)<<1):(Yi[p]<<1));
		
			//if we're on a pv bus, the second column only needs entries for the diagonal
			if(is_pv){
				if(Yi[p]==j){				
					Ji[nz++]=ij;
					Jx[nz-1]=creal(Yx[p]);
					Ji[nz++]=ij+1;
					Jx[nz-1]=-cimag(Yx[p]);
					break;
				}
			}
			//if not on a pv column, just add in ybus entries
			else{		
				Ji[nz++]=ij;
				Jx[nz-1]=(col?creal(Yx[p]):cimag(Yx[p]));
				Ji[nz++]=ij+1;
				Jx[nz-1]=(col?-cimag(Yx[p]):creal(Yx[p]));			
			}
		}
		//set the next column pointer
		Jp[((j>swing_bus?j-1:j)<<1)+col+1]=nz;
		
		//check whether we're on a pv column (in ybus) so that we can
		//build the next jacobian column correctly.
		if((is_pv=(j==gens[pv].index)))pv++;
		if(gens[pv].index==swing_bus&&pv<num_gens)pv++;
		col=!col;		
	}
	//resize the matrix, freeing up any unused memory and return
	cs_di_sprealloc(jac, nz);
	return jac;
}

//residual computation function
int ComputeF(double *v, double *r, PS ps){
	
	//function variables
	cs_ci *ybus=ps->ybus;
	int j,jp,p,i,n=ybus->n,*Yp=ybus->p,*Yi=ybus->i,swing_bus=ps->swing_bus, pv=0,sh=0,num_gens=ps->num_gens,num_shunts=ps->num_shunts;
	double vm,vr,vg,vk2,P,Q;
	double _Complex *Yx=ybus->x;
	Shunt *shunts=ps->shunts;
	Gen *gens=ps->gens;
		
	//increment shunt and gen pointers if their first buses are the swing bus
	if(shunts[sh].index==swing_bus)sh++;
	if(gens[pv].index==swing_bus)pv++;
	for(j=n;j<(n-1)<<1;j++)r[j]=0.0;
	
	//loop over the columns of ybus to calculate I=Ybus*V+I_l
	for(j=0;j<n;j++){
	
		//skip swing bus
		if(j==swing_bus) continue;	
		jp=(j>swing_bus?((j-1)<<1):(j<<1));
		r[jp]=r[jp+1]=0.0;
		
		//add generator currents
		vr=v[jp];vm=v[jp+1];vk2=vr*vr+vm*vm;
		if(pv<num_gens&&j==gens[pv].index){
			vg=gens[pv].vsp; vm=vr; vk2=vg*vg; vr=sqrt(vg*vg-vm*vm);
			P=gens[pv].pg;
			Q=gens[pv].qg;
			r[jp]+=(P*vm-Q*vr)/vk2;r[jp+1]+=(P*vr+Q*vm)/vk2;
			pv++;
			if(gens[pv].index==swing_bus&&pv<num_gens)pv++;
		}		
		//add load currents
		if(sh<num_shunts&&j==shunts[sh].index){
			P=shunts[sh].p;
			Q=shunts[sh].q;
			r[jp]-=(P*vm-Q*vr)/vk2;r[jp+1]-=(P*vr+Q*vm)/vk2;
			sh++;
			if(sh<num_shunts&&shunts[sh].index==swing_bus)sh++;
		}
		//add ybus currents
		for(p=Yp[j];p<Yp[j+1];p++){
			if(Yi[p]==swing_bus) continue;
			i=(Yi[p]>swing_bus?((Yi[p]-1)<<1):(Yi[p]<<1));
			r[i]-=(creal(Yx[p])*vm+cimag(Yx[p])*vr);
			r[i+1]-=(creal(Yx[p])*vr-cimag(Yx[p])*vm);
		}
	}
	//correct pv buses: P=VrIr+VmIm.
	for(pv=0;pv<num_gens;pv++){
		if(gens[pv].index==swing_bus)continue;		jp=(gens[pv].index>swing_bus?((gens[pv].index-1)<<1):(gens[pv].index<<1));
		vm=v[jp];vg=gens[pv].vsp;vk2=vg*vg; vr=sqrt(vg*vg-vm*vm);
		P=vr*r[jp]+vm*r[jp+1];
		r[jp]=vm*P/vk2;r[jp+1]=vr*P/vk2;
	}	
	//return
	return 0;
}

//jacobian computation function
int ComputeJac(double *v, cs_di* jac, PS ps){
	
	//input check
	if(!jac) return 1;
	
	//set up appropriate function variables
	cs_ci *ybus=ps->ybus;
	int p, pj, j, jp, i, swing_bus=ps->swing_bus, *Jp=jac->p, *Yp=ybus->p, *Yi=ybus->i, *src_buses=ps->src_buses, num_src=ps->num_src, num_shunts=ps->num_shunts, num_gens=ps->num_gens, pv=0, sh=0, is_pv; 
	Gen *gens=ps->gens;
	Shunt *shunts=ps->shunts;
	double *Jx=jac->x, vr, vm, vg, vk2, vk6, P, Q, b1, b2, g1, g2, g, b;
	double _Complex *Yx=ybus->x;
	
	//increment shunt and gen pointers if their first buses are the swing bus
	if(shunts[sh].index==swing_bus)sh++;
	if(gens[pv].index==swing_bus)pv++;
	
	//loop over load/generation buses
	for(i=0;i<num_src;i++){
	
		//skip the swing bus
		if((j=src_buses[i])==swing_bus)continue;
		jp=(j>swing_bus?(j-1):j)<<1;
		
		//determine load/generation on this bus
		is_pv=0; P=0.0; Q=0.0;
		if(sh<num_shunts&&shunts[sh].index==j){
			P+=shunts[sh].p;
			Q+=shunts[sh].q;
			sh++;
			if(sh<num_shunts&&shunts[sh].index==swing_bus)sh++;
		}
		if(pv<num_gens&&gens[pv].index==j){
			is_pv=1;
			P-=gens[pv].pg;
			Q-=gens[pv].qg;
			vg=gens[pv].vsp;
			pv++;
			if(pv<num_gens&&gens[pv].index==swing_bus)pv++;
		}
		
		//get the current guess for bus voltage
		vr=v[jp]; vm=v[jp+1];
		if(is_pv){
			vm=vr; vr=sqrt(vg*vg-vm*vm);
		}
		vk2=vr*vr+vm*vm;vk6=vk2*vk2*vk2;
		
		//update diagonals and pv columns as necessary
		for(p=Yp[j],pj=Jp[jp]; p<Yp[j+1]; p++,pj+=2){
			
			//move on if we're on the swing bus row
			if(Yi[p]==swing_bus){p++;if(p==Yp[j+1])break;}
			
			//non pv bus columns only need diagonals updated
			if(!(is_pv||Yi[p]==j))continue;
			pj=Jp[jp];g=creal(Yx[p]);b=cimag(Yx[p]);
			
			//update diagonals
			if(Yi[p]==j){
				g1=g2=g;
				b1=-(b2=b);
				g1+=(P-2*vr*vk2*(P*vr+Q*vm))/vk6;
				b1+=(Q-2*vm*vk2*(P*vr+Q*vm))/vk6;
				b2-=(Q+2*vr*vk2*(P*vm-Q*vr))/vk6;
				g2+=(P-2*vm*vk2*(P*vm-Q*vr))/vk6;
				if(is_pv){
					Jx[pj]=g2-b2*vm/vr;
					Jx[pj+1]=b1-g1*vm/vr;
					Jx[Jp[jp+1]]=vr/vk2;
					Jx[Jp[jp+1]+1]=-vm/vk2;
					continue;
				}
				Jx[pj]=b2;
				Jx[pj+1]=g1;
				Jx[pj+(Jp[jp+1]-Jp[jp])]=g2;
				Jx[pj+(Jp[jp+1]-Jp[jp])+1]=b1;
				continue;
			}
			//otherwise, update first pv column
			Jx[pj]=-b-g*vm/vr;
			Jx[pj+1]=g-b*vm/vr;
		}
	}	
	//return
	return 0;
}