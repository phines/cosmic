/*
 file: psc.c
 author: Ben O'Hara
 email: buohara@gmail.com, buohara@uvm.edu
 description: Implementation of power system data-related routines.
 */
#include "psc.h"

//comparison functions for sorting items by their bus indices
COMPARE(Shunt)
COMPARE(Gen)
COMPARE(Bus)

//some readability macros
#define BUS(x) &(ps->buses[i].x)
#define BR(x) &(ps->branches[i].x)
#define GEN(x) &(ps->gens[i].x)
#define SH(x) &(ps->shunts[i].x)

/**
 * Create the system admittance matrix from power system data. Note that this  function is called transparently from GetPSFromFile once power system data has been populated (if the file read is successful).
 *
 * @param ps Power system data that will be used to construct YBus.
 */
static int CreateYBus(PS ps, int *ids);

/**
 * Bus comparison function used to sort buses by ID in a ps struct.
 * 
 * @param bus_1 First bus to be compared.
 * @param bus_2 Second bus to be compared.
 * 
 * @return A negative, zero, or positive integer if bus_1 is less than,
 equal, or greater than bus_2, respectively.
 */
static int CompareBuses(const void *bus_1, const void *bus_2);

/*
 * Read power system data from a file and load it into a PS struct.
 */
PS GetPS(const char* file_name){
	
	//numbers of various items in the ps file.
	int i, num_buses, num_branches, num_gens, num_shunts, *ids, pv=0, sh=0, src=0;
	double base_mva;
	PS ps;
	
	//check input
	if(file_name == NULL) return NULL;
	FILE *ps_file = fopen(file_name, "r");
	if(ps_file == NULL) return NULL;
	
	//read off ps data sizes
	if(fscanf(ps_file, "BASE_MVA %lg\n", &base_mva) != 1) return NULL;
	if(fscanf(ps_file, "BUS %d\n", &num_buses) != 1) return NULL;
	if(fscanf(ps_file, "BRANCH %d\n", &num_branches) != 1) return NULL;
	if(fscanf(ps_file, "GEN %d\n", &num_gens) != 1) return NULL;
	if(fscanf(ps_file, "SHUNT %d\n", &num_shunts) != 1) return NULL;
	
	ps = malloc(sizeof(_PS));
	if(!ps) return NULL;
	ps->base_mva = base_mva;
	ps->num_buses = num_buses;
	ps->num_branches = num_branches;
	ps->num_gens = num_gens;
	ps->num_shunts = num_shunts;
	ps->num_macs = num_gens;
	
	//allocate memory for buses, branches, etc.
	ps->buses = malloc(num_buses*sizeof(Bus));	
	ps->branches = malloc(num_branches*sizeof(Branch));
	ps->gens = malloc(num_gens*sizeof(Gen));
	ps->shunts = malloc(num_shunts*sizeof(Shunt));
	ids = malloc(num_buses*sizeof(int));
	ps->src_buses = malloc((num_shunts+num_gens)*sizeof(int));
	
	//check out memory allocations
	if(!ps->buses || !ps->branches || !ps->gens || !ps->shunts || !ids || !ps->src_buses){
		FreePS(&ps);
		if(ids) free(ids);
		return NULL;
	}
	
	//read data from file and populate ps 
	fscanf(ps_file, "\n");
	for(i = 0; i < num_buses; i++){
		fscanf(ps_file, "%d %d %lg %lg %lg %lg %d %lg %lg %lg %d %lg %lg\n", BUS(id), BUS(type), BUS(pd), BUS(qd), BUS(gs), BUS(bs), BUS(area), BUS(vmag), BUS(vang), BUS(base_kv), BUS(zone), BUS(vmax), BUS(vmin));
		ids[i] = ps->buses[i].id;	
	}
	 
	//sort the bus array, get local indices, and locate the swing bus index
	qsort(ps->buses, num_buses, sizeof(Bus), CompareBuses);
	for(i = 0; i < num_buses; i++){ 
		ps->buses[i].index = i;
		ids[i] = ps->buses[i].id;
		if(ps->buses[i].type == 3) ps->swing_bus = i;	
	}
	
	//branches
	int swap_temp;
	fscanf(ps_file, "\n");
	for(i = 0; i < num_branches; i++){
		fscanf(ps_file, "%d %d %lg %lg %lg %lg %lg %lg\n", BR(from), BR(to), BR(r), BR(x), BR(b), BR(rate_a), BR(rate_b), BR(rate_c));
		
		//check whether from is greater than to, and if so, swap. this
		//helps identify parallel branches.
		if(ps->branches[i].from > ps->branches[i].to){
			swap_temp = ps->branches[i].from;
			ps->branches[i].from = ps->branches[i].to;
			ps->branches[i].to = swap_temp;
		}
		ps->branches[i].from_index = (int*)bsearch(&(ps->branches[i].from), ids, num_buses, sizeof(int), CompareBuses) - ids;
		ps->branches[i].to_index = (int*)bsearch(&(ps->branches[i].to), ids, num_buses, sizeof(int), CompareBuses) - ids;
	}
	
	//generators
	fscanf(ps_file, "\n");
	for(i = 0; i < num_gens; i++){
		fscanf(ps_file, "%d %lg %lg %lg %lg %lg %lg %d %lg %lg\n", GEN(bus), GEN(pg), GEN(qg), GEN(qmax), GEN(qmin), GEN(vsp), GEN(m_base), GEN(status), GEN(pmax), GEN(pmin));
		ps->gens[i].index = (int*)bsearch(&(ps->gens[i].bus), ids, num_buses, sizeof(int), CompareBuses) - ids;
		ps->gens[i].pg/=base_mva;
		ps->gens[i].qg/=base_mva;
	}
	qsort(ps->gens, ps->num_gens, sizeof(Gen), CompareGen);
	
	//shunts
	fscanf(ps_file, "\n");
	for(i = 0; i < num_shunts; i++){
		fscanf(ps_file, "%d %lg %lg\n", SH(bus), SH(p), SH(q));
		ps->shunts[i].index = (int*)bsearch(&(ps->shunts[i].bus), ids, num_buses, sizeof(int), CompareBuses) - ids;
		ps->shunts[i].p/=base_mva;
		ps->shunts[i].q/=base_mva;
	}
	qsort(ps->shunts, ps->num_shunts, sizeof(Shunt), CompareShunt);
	
	//find the set union of load and generator buses
	while(sh < num_shunts || pv < num_gens){
		if(pv==num_gens){ps->src_buses[src++] = ps->shunts[sh++].index; continue;}
		if(sh==num_shunts){ps->src_buses[src++] = ps->gens[pv++].index; continue;}
		if(ps->shunts[sh].index==ps->gens[pv].index){ 
			ps->src_buses[src++] = ps->shunts[sh++].index;
			pv++;
			continue;
		}
		ps->src_buses[src++] = (ps->shunts[sh].index < ps->gens[pv].index ? ps->shunts[sh++].index : ps->gens[pv++].index);
	}
	
	//free up any extra memory in the source bus array
	realloc(ps->src_buses, src*sizeof(int));
	ps->num_src = src;
	
	//create ybus. if there are any errors, free allocated memory and
	//return NULL
	if(CreateYBus(ps, ids)){
		FreePS(&ps);
		free(ids);
		return NULL;
	}
	
	//release the bus map
	free(ids); 
	 
	//return success
	return ps;
};

/*
 Create the system admittance matrix from power system data. Note that this function is called transparently from GetPSFromFile once power system data has been populated (if the file read is successful).
 */
int CreateYBus(PS ps, int *ids){

	int i, from, to, num_buses = ps->num_buses, num_branches = ps->num_branches, nnz = num_buses + 2*num_branches;
	double r, x, b;
	double _Complex y, *diags = malloc(num_buses*sizeof(double _Complex));
	
	//allocate a ybus
	cs_ci *ybus = cs_ci_spalloc(num_buses, num_buses, nnz, 1, 1), *ybus_comp, *ybus_comp_t;
	
	//check memory
	if(!ybus) return 1;
		
	//loop through ybus branches
	for(i=0; i<num_branches; i++)
	{
		//get branch from/to bus indices
		from = ps->branches[i].from_index;
		to = ps->branches[i].to_index;
	
		//determine line admittance
		r = ps->branches[i].r;
		x = ps->branches[i].x;
		b = ps->branches[i].b;		
		y = (1.0+0.0*I)/(r+x*I);

		//add off-diagonals to ybus and increment diagonal entries
		cs_ci_entry(ybus, from, to, -y);
		cs_ci_entry(ybus, to, from, -y);
		diags[from]+=(y+I*b/2.0);
		diags[to]+=(y+I*b/2.0);
	}
		
	//add diogonal entries to ybus.
	for(i=0; i<num_buses; i++) cs_ci_entry(ybus, i, i, diags[i]);
		
	//compress ybus
	ybus_comp=cs_ci_compress(ybus);	
	
	//check that memory allocation was okay
	if(!ybus_comp){cs_ci_spfree(ybus); return 1;}
		
	//free the old ybus and sum duplicates in the compressed matrix (i.e., parallel branches)
	cs_ci_spfree(ybus);
	cs_ci_dupl(ybus_comp);
	
	//sort the columns of ybus using the transpose-transpose trick
	ybus_comp_t = cs_ci_transpose(ybus_comp, 1);
	if(!ybus_comp_t){free(ybus_comp); return 1;}
	free(ybus_comp);
	ybus_comp = cs_ci_transpose(ybus_comp_t, 1);
	if(!ybus_comp){free(ybus_comp_t); return 1;}
	free(ybus_comp_t);
	
	//assign the ps ybus pointer to the compressed, summed ybus matrix	
	ps->ybus = ybus_comp;
		
	//free the diagonal workspace and return
	free(diags);
	return 0;
}

//Trip/restore a power system transimission line.
int ModifyBranch(PS ps, int to, int from, int line_num){
	return 0;
};

/*
 * Free memory allocated for power system data.
 */
int FreePS(PS *ps){
	if((*ps)->buses){free((*ps)->buses); (*ps)->buses = NULL;}
	if((*ps)->branches){free((*ps)->branches); (*ps)->branches = NULL;}
	if((*ps)->gens){free((*ps)->gens); (*ps)->gens = NULL;}
	if((*ps)->shunts){free((*ps)->shunts);(*ps)->shunts = NULL;}
	if((*ps)->machines){free((*ps)->machines);(*ps)->machines = NULL;}
	if((*ps)->src_buses){free((*ps)->src_buses); (*ps)->src_buses=NULL;}
	if((*ps)->ybus) cs_ci_spfree((*ps)->ybus);
	*ps = NULL;
	return 0;
};

//compare buses by ID
int CompareBuses(const void *bus_1, const void *bus_2){
	return ((Bus*)bus_1)->id - ((Bus*)bus_2)->id;
}

#ifdef _VERBOSE

//some more readability constants. now i'm just getting lazy.
#define BUSp(x) ps->buses[i].x
#define BRp(x) ps->branches[i].x
#define GENp(x) ps->gens[i].x
#define SHp(x) ps->shunts[i].x

/*
 Print power system data.
 */
void PrintPS(PS ps){

	//title
	printf("PS Data:\n\n");
	
	//base MVA
	printf("Base MVA: %lg\n\n", ps->base_mva);

	int i;
	printf("Buses:\n");
	for(i = 0; i < ps->num_buses; i++){
		printf("%d %d %d %lg %lg %lg %lg %d %lg %lg %lg %d %lg %lg\n", BUSp(id), BUSp(index), BUSp(type), BUSp(pd), BUSp(qd), BUSp(gs), BUSp(bs), BUSp(area), BUSp(vmag), BUSp(vang), BUSp(base_kv), BUSp(zone), BUSp(vmax), BUSp(vmin));
	}
	 
	//branches
	printf("\nBranches:\n");
	for(i = 0; i < ps->num_branches; i++){
		printf("%d %d %d %d %lg %lg %lg %lg %lg %lg\n", BRp(from), BRp(to), BRp(from_index), BRp(to_index), BRp(r), BRp(x), BRp(b), BRp(rate_a), BRp(rate_b), BRp(rate_c));
	}
	
	//generators
	printf("\nGenerators:\n");
	for(i = 0; i < ps->num_gens; i++){
		printf("%d %d %lg %lg %lg %lg %lg %lg %d %lg %lg\n", GENp(bus), GENp(index), GENp(pg), GENp(qg), GENp(qmax), GENp(qmin), GENp(vsp), GENp(m_base), GENp(status), GENp(pmax), GENp(pmin));
	}
	
	//shunts
	printf("\nShunts:\n");
	for(i = 0; i < ps->num_shunts; i++){
		printf("%d %d %lg %lg\n", SHp(bus), SHp(index), SHp(p), SHp(q));
	}
	printf("\n");

	//return
	return;
}

#endif