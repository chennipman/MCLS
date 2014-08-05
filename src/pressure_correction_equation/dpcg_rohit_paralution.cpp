#include "dpcg_rohit_paralution.hpp"


int call_to_cg_wrapper(Array2<double> A,Array1<double> x,Array1<double> b, int max_iter, double tolerance,
			int xdim, int ydim, int zdim, int *iter_cnt, double *l2nrm_residual,
		       Array3<double> level_set)
{
  FILE *fp;
  double *Aptr=NULL, *bin=NULL, *xin=NULL, *acsrvals=NULL;
  double *levelsetptr=NULL, *xout=NULL, *levref=NULL, *xref=NULL, *bref=NULL;
  int *acsrcols=NULL, *acsrrows=NULL, status;
  int dim=xdim*ydim*zdim, nnzA, i;
  int setlssd_in=1, defvex_perdirec_in=8, lvst_offst_in=1;
  int phisize=(xdim+2*lvst_offst_in)*(ydim+2*lvst_offst_in)*(zdim+2*lvst_offst_in);
  struct timeval now;
  double tick, tack;
  allocate_host(dim*4,&Aptr);
  allocate_host(dim,&xout);
  allocate_host(dim,&xin);
  allocate_host(dim,&bin);
  
  bref=b.get_pointer();  xref=x.get_pointer();
  printf("\n GPU memory status on entry is\n");
  get_memory_status();
  
  for(i=0;i<dim;i++)
  {
    xin[i]=xref[i];	bin[i]=bref[i];
  }
//     }
  if(setlssd_in){
    //levelsetptr=(double*)malloc(phisize*sizeof(double));
//     levelsetptr = new double [phisize];
    allocate_host(phisize, &levelsetptr);
    //std::copy(level_set.get_pointer(), level_set.get_pointer()+phisize,levelsetptr);
//     memcpy(levelsetptr,level_set.get_pointer(),sizeof(double)*phisize);
    levref=level_set.get_pointer();
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
    
  status=call_paralution_dpcg(acsrvals, acsrrows, acsrcols, bin ,xin, max_iter,
			      tolerance, nnzA, dim, iter_cnt, l2nrm_residual, xdim,
			      ydim, zdim, setlssd_in, defvex_perdirec_in, lvst_offst_in,
			      phisize, levelsetptr, xout);
//   printf("\n Solver status has been reported as %d",status);
  //std::copy(xout, xout+dim,x.get_pointer());
  for(i=0;i<dim;i++)
    x[i]=xout[i];
//   }
//   printf("\n going to free stuff now");

  free_host(&Aptr);	free_host(&xout);
  if(setlssd_in) 
    free_host(&levelsetptr);
  printf("\n GPU memory status on exit is\n");
  get_memory_status();
  
  if(status)
    return 0;//all non-zero status from paralution means success
  else
    return 1;
  
}