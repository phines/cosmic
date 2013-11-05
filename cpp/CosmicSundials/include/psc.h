/*
 file: psc.c
 author: Ben O'Hara
 email: buohara@gmail.com, buohara@uvm.edu
 description: Data structures that contain various power system information, including
 collections of bus data, branch data, load data, and machine data.
 */

#ifndef PSC_H
#define PSC_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <cs.h>

#define COMPARE(x) int Compare##x(const void *item1, const void *item2){ \
return ((x*)item1)->index - ((x*)item2)->index; \
}

/// Bus data structure.
typedef struct{
	int id;
	int index;
	int type;
	double pd;
	double qd;
	double gs;
	double bs;
	int area;
	double vmag;
	double vang;
	double base_kv;
	int zone;
	double vmax;
	double vmin;
	
	/*
	double lam_p;
	double lam_q;
	double mu_vmax;
	double mu_vmin;
	int loc_x;
	int loc_y;
	*/
} Bus;

/// Branch data structure.
typedef struct{
	int from;
	int to;
	int from_index;
	int to_index;
	double r;
	double x;
	double b;
	double rate_a;
	double rate_b;
	double rate_c;
	
	/*
	int tap;
	double shift;
	int status;
	double pf;
	double qf;
	double pt;
	double qt;
	double mu_f;
	double mu_t;
	double imag_f;
	double imag_t;
	int switchable;
	double prob_fail;
	int type;
	*/
} Branch;

/// Generator data structure.
typedef struct{
	int bus;
	int index;
	double pg;
	double qg;
	double qmax;
	double qmin;
	double vsp;
	double m_base;
	int status;
	double pmax;
	double pmin;
	/*
	double mu_pmax;
	double mu_pmin;
	double mu_qmax;
	double mu_qmin;
	int type;
	double cost;
	double part_fact;
	double ramp_rate_up;
	double ramp_rate_down;
	*/
} Gen;

/// Machine data structure.
typedef struct{
	int gen;
	double r;
	double xd;
	double xdp;
	double xdpp;
	double xq;
	double xqp;
	double xqpp;
	double d;
	double m;
	double ea;
	double eap;
	double pm;
	double pm0;
	double delta_m;
	double omega;
	double td0;
	double td0p;
}Mac;

/*
///Exciter data structure.
typedef struct{
	int gen;
	int type;
	double ka;
	double ta;
	double tb;
	double ke;
	double te;
	double kf;
	double tf;
	double aex;
	double bex;
	double urmin;
	double urmax;
	double vref;
	double efd;
	double e1;
}Exc;

///Governor data structure.
typedef struct{
	int gen;
	int type;
	double r;
	double tsr;
	double tsm;
	double tt;
	double lcmax;
	double lcmin;
	double pmax;
	double pmin;
	double pref;
	double p1;
	double p2;
}Gov;
*/

/// Load/shunt data structure.
typedef struct{
	int bus;
	int index;
	double p;
	double q;
	/*
	double frac_s;
	double frac_z;
	double frac_y;
	int status;
	int type;
	double value;
	double frac_e;
	double gamma;
	*/
}Shunt;

/**
 * Power system data structure. This contains all data relevant to a power
 * system, including all bus, branch, load, and machine data.
 */
typedef struct{
	double base_mva;
	int num_buses;
	int num_branches;
	int num_gens;
	int num_shunts;
	int num_macs;
	//int num_govs;
	int swing_bus;
	int num_src;
	
	Bus *buses;
	Branch *branches;
	Gen *gens;
	Shunt *shunts;
	Mac *machines;
	//Gov *governors;
	int *src_buses;
	cs_ci *ybus;
} _PS, *PS;

/**
 * Read power system data from a file and load it into a PS struct.
 *
 * @param file_name The path to the power system data file.
 * @param ps Power system struct to be filled with data.
 *
 * @return Status code indicating result of file read.
 */
PS GetPS(const char* file_name);

/**
 * Trip/restore a power system transimission line.
 *
 * @param ps Power system whose branch data is to be modified.
 * @param to Destination side of the modified transmission line.
 * @param from Source side of the modified transmission line.
 * @param line_num If there are parallel lines, which line is to be modified.
 *
 * @return Status code indicating result of line modification.
 */
int ModifyBranch(PS ps, int to, int from, int line_num);

/**
 * Free memory allocated for power system data.
 *
 * @param ps The power system data to be freed.
 *
 * @return Status code indicating result of memory deallocation.
 */
int FreePS(PS *ps);

#ifdef _VERBOSE

/**
 * Print power system data.
 * 
 * @param ps Power system to print.
 */
void PrintPS(PS ps);

#endif
	
#ifdef __cplusplus
}
#endif

#endif