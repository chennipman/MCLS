#pragma once
#include<cstdlib>
#include<iostream>
#include<algorithm>
#include<math.h>
#include<stdio.h>
#include<string.h>
// #include<stdlib.h>
#include <string>
#include <sstream>
#include <fstream>
#include<math.h>
#include <sys/time.h>
#include "mmio.h"
#include "omp.h"
#include "paralution.hpp"
#include "driver_types.h"
using namespace std;
using namespace paralution;
// #define GPURUN	10000
#define BUBFLO	20000
extern void  wrap_A_intoDIA(Array2<double>, double *, int, int, int, int *);
extern void  cnvrtDIA_to_CSR(double *, double **, int **, int **, int, int, int, int);


extern int  call_paralution_dpcg(double *, int *, int *, double*, double *, int, double,
				 int, int, int *, double *, int xdim, int ydim, int zdim, int setlssd_in,
			  int defvex_perdirec_in, int lvst_offst_in, int phisize, double *lvlstptr, double *xout);


extern int call_to_cg_wrapper(Array2<double> A,Array1<double> x,Array1<double> b, int max_iter, double tolerance,
			int xdim, int ydim, int zdim, int *iter_cnt, double *l2nrm_residual, Array3<double>);

extern void get_memory_status();
