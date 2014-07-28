#include "dpcg_rohit_paralution.hpp"


int call_to_cg_wrapper(Array2<double> A,Array1<double> x,Array1<double> b, int max_iter, double tolerance,
			int xdim, int ydim, int zdim, int *iter_cnt, double *l2nrm_residual,
		       Array3<double> level_set)
{
  FILE *fp;
  double *Aptr, *bin, *xin, *acsrvals, *levelsetptr, *xout, *levref, *xref, *bref;
  int *acsrcols, *acsrrows, status;
  int dim=xdim*ydim*zdim, nnzA, i;
  int setlssd_in=1, defvex_perdirec_in=8, lvst_offst_in=1;
  int phisize=(xdim+2*lvst_offst_in)*(ydim+2*lvst_offst_in)*(zdim+2*lvst_offst_in);
  struct timeval now;
  double tick, tack;
  Aptr=(double*)calloc(4*dim,sizeof(double));
  xout=(double*)calloc(dim,sizeof(double));
  xin = (double*)calloc(dim,sizeof(double));
  bin = (double*)calloc(dim,sizeof(double));

  bref=b.get_pointer();  xref=x.get_pointer();
// #pragma omp parallel
//     {  
// #pragma omp for private(i)     
  for(i=0;i<dim;i++)
  {
    xin[i]=xref[i];	bin[i]=bref[i];
  }
//     }
//   memcpy(xin,x.get_pointer(),sizeof(double)*dim);
//   memcpy(bin,b.get_pointer(),sizeof(double)*dim);
  if(setlssd_in){
    levelsetptr=(double*)calloc(phisize,sizeof(double));
    //std::copy(level_set.get_pointer(), level_set.get_pointer()+phisize,levelsetptr);
//     memcpy(levelsetptr,level_set.get_pointer(),sizeof(double)*phisize);
    levref=level_set.get_pointer();
// #pragma omp parallel
//     {
// #pragma omp for 
    for(i=0;i<phisize;i++)
      levelsetptr[i]=(levref)[i]*(-1.0f);

    //fp=fopen("raw_lvlset_irreg.rec","wt");
//     ofstream levelsetfile;
//     levelsetfile.open("raw_lvlset_irreg.rec",ios::out);
//     for(i=0;i<phisize;i++);
//       levelsetfile<<levelsetptr[i]<<std::endl;
//     levelsetfile.close();
//     for(i=0;i<phisize;i++);
//     {
//       fprintf(fp,"%0.9f\n",levelsetptr[i]);
//       
//     }
//     fclose(fp);
//     }
  }
  
  gettimeofday(&now, NULL);
  tick = now.tv_sec*1000000.0+(now.tv_usec);
  
  wrap_A_intoDIA(A, Aptr, xdim, ydim, zdim, &nnzA);
  cnvrtDIA_to_CSR(Aptr, &acsrvals, &acsrcols, &acsrrows, xdim, ydim, zdim, nnzA);
  
  gettimeofday(&now, NULL);
  tack = now.tv_sec*1000000.0+(now.tv_usec);
  std::cout << "A conversion:" << (tack-tick)/1000000 << " sec" << std::endl;
  printf("\n GPU memory status on entry is\n");
  get_memory_status();
  
  status=call_paralution_dpcg(acsrvals, acsrrows, acsrcols, bin ,xin, max_iter,
			      tolerance, nnzA, dim, iter_cnt, l2nrm_residual, xdim,
			      ydim, zdim, setlssd_in, defvex_perdirec_in, lvst_offst_in,
			      phisize, levelsetptr, xout);
//   printf("\n Solver status has been reported as %d",status);
  //std::copy(xout, xout+dim,x.get_pointer());
// #pragma omp parallel
//   {
// #pragma omp for  
  for(i=0;i<dim;i++)
    x[i]=xout[i];
//   }
//   printf("\n going to free stuff now");
  free(Aptr);	free(xout);	
  if(setlssd_in) 
    free(levelsetptr);
  printf("\n GPU memory status on exit is\n");
  get_memory_status();
  
  if(status)
    return 0;//all non-zero status from paralution means success
  else
    return 1;
  
}