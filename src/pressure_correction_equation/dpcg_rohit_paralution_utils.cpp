#include "dpcg_rohit_paralution.hpp"
using namespace std;


// void  wrap_A_intoDIA(Array2<double>, double *, int, int, int, int *);
// void  cnvrtDIA_to_CSR(double *, double **, int **, int **, int, int, int, int);
// int  call_paralution_dpcg(double *, int *, int *, double*,  double *,
// 			  int, double, int, int, int *, double *);

extern int mm_write_banner(FILE *f, MM_typecode matcode);
extern int mm_write_mtx_crd_size(FILE *f, int M, int N, int nz);
extern int mm_write_mtx_array_size(FILE *f, int M, int N);

void initialize_paralution_library(void)
{
  init_paralution();
}
void stop_paralution_library(void)
{
  stop_paralution();
}
void wrap_A_intoDIA(Array2<double> Adblptr, double *Aptr, int xdim, int ydim, int zdim, int *nnzA)
{
  int dim=xdim*ydim*zdim, nsqr=xdim*ydim, i,j,k,l, lclnnzmain, lclnnz1, lclnnz2, lclnnz3;
  int max_threads, num_procs;
  MM_typecode matcode;   
  FILE *fp;
  lclnnz1=0; lclnnz2=0; lclnnz3=0; lclnnzmain=0;
  
//   max_threads=omp_get_max_threads();
//     num_procs= omp_get_num_procs();
//     printf("\n Number of thread for openMP inside makebubmap are %d. max_threads are %d. num_procs %d", omp_get_num_threads(), max_threads, num_procs);
//     omp_set_dynamic(0);
//     omp_set_num_threads(8);
//     printf("\n Number of thread for openMP inside makebubmap are now %d", omp_get_num_threads());
//  #pragma omp for   
  for(i=0;i<dim-(nsqr);i++)
  {
    Aptr[i+nsqr]=Adblptr[3][i];//storing diagonal with offset -n*n.
    if(fabs(Aptr[i+nsqr])>0.0f)lclnnz1++;
  }
//  #pragma omp for   
  for(j=0;j<dim-xdim;j++)
  {
    Aptr[dim+j+xdim]=Adblptr[2][j];//storing diagonal with offset -n.
    if(fabs(Aptr[dim+j+xdim])>0.0f)lclnnz2++; //it is stored the opposite way since Jok stores upper triangular part
  }
//  #pragma omp for   
  for(k=0;k<dim-1;k++)
  {
    Aptr[2*dim+k+1]=Adblptr[1][k];//storing diagonal with offset -1.
    if(fabs(Aptr[2*dim+k+1])>0.0f)lclnnz3++;//it is stored the opposite way since Jok stores upper triangular part
  }
//   printf("\n Num openmp threads is %d", omp_get_num_threads());
//  #pragma omp for
  for(l=0;l<dim;l++)
  {
    Aptr[3*dim+l]=Adblptr[0][l];//storing main diagonal
    if(fabs(Aptr[3*dim+l])>0.0f)lclnnzmain++;
  }
      
//   for(i=0;i<dim;i++)
//   {
//     Aptr[4*dim+i]=Adblptr[1][i];//storing diagonal with offset +1.
//     if(fabs(Aptr[2*dim+i])>0.0f)lclnnz++; //it is stored the expected way since Jok stores upper triangular part
//   }
//   for(i=0;i<dim;i++)
//       Aptr[5*dim+i]=Adblptr[2][i];//storing diagonal with offset +n.
// 				    //it is stored the expected way since Jok stores upper triangular part
//   for(i=0;i<dim;i++)
//       Aptr[6*dim+i]=Adblptr[3][i];//storing diagonal with offset +n.
// 				    //it is stored the expected way since Jok stores upper triangular part
 *nnzA=2*lclnnz1 + 2*lclnnz2 + 2*lclnnz3 + lclnnzmain;
 /*  fp=fopen("AC16_duncan.mtx","wt");
    mm_initialize_typecode(&matcode);
    mm_set_matrix(&matcode);
    mm_set_coordinate(&matcode);
    mm_set_real(&matcode);
    mm_write_banner(fp, matcode); 
	       //mm_write_mtx_crd_size(fp, M, N, nz);
//     nnz=0;
//     for(i=0;i<dim;i++)
//       for(j=0;j<ndiags;j++)
// 	if(fabs(Aptr[j*dim+i])!=0.0)
// 	  nnz++; 
	
    mm_write_mtx_crd_size(fp, dim, dim, *nnzA);
    for(i=0;i<dim;i++)
      {
	if(fabs(Aptr[i])!=0.0 && i-xdim*xdim>=0)
	      fprintf(fp,"%d %d %f\n",i-xdim*xdim+1,i+1,Aptr[i]);//diagonal with offset -n*n
	      
	if(fabs(Aptr[dim+i])!=0.0 && i-xdim>=0)
	      fprintf(fp,"%d %d %f\n",i-xdim+1,i+1,Aptr[dim+i]);//diagonal with offset -n
      if(fabs(Aptr[2*dim+i])!=0.0 && i-1>=0)
	      fprintf(fp,"%d %d %f\n",i,i+1,Aptr[2*dim+i]);//diagonal with offset -1
      if(fabs(Aptr[3*dim+i])!=0.0)
		      fprintf(fp,"%d %d %f\n",i+1,i+1,Aptr[3*dim+i]);//main diagonal
      if(i+1<=dim-1)
	if(fabs(Aptr[2*dim+i+1])!=0.0)
		      fprintf(fp,"%d %d %f\n",i+2,i+1,Aptr[2*dim+i+1]);//diagonal with offset +1
      if(i+xdim<=dim-1)
	if(fabs(Aptr[dim+i+xdim])!=0.0)
		      fprintf(fp,"%d %d %f\n",i+xdim+1,i+1,Aptr[dim+i+xdim]);//diagonal with offset +n
      if(i+xdim*xdim<=dim-1)
	if(fabs(Aptr[i+xdim*xdim])!=0.0)
		      fprintf(fp,"%d %d %f\n",i+xdim*xdim+1,i+1,Aptr[i+xdim*xdim]);//diagonal with offset +n*n

    }
   fclose(fp);//*/
  
  
  
}

void  cnvrtDIA_to_CSR(double *Aptr, double **acsrvals, int **acsrcols, int **acsrrows, int xdim,
			int ydim, int zdim, int nnzA)
{
  int dim=xdim*ydim*zdim, i, nsqr, lclctr, csridx;
  nsqr=xdim*ydim;
  *acsrvals=(double*)calloc(nnzA,sizeof(double));
  *acsrcols=(int*)calloc(nnzA,sizeof(int));
  *acsrrows=(int*)calloc(dim+1,sizeof(int));
  
  //convert DIA to CSR now
  (*acsrrows)[0]=0;
  omp_set_num_threads(8);
  lclctr=0; csridx=0;
#pragma omp for   
  for(i=0;i<dim;i++)
  {
    if((i-nsqr)>=0)// diagonal with offset -n*n avlbl
    {
      if(fabs(Aptr[i])>0.0f)
      {
	(*acsrvals)[csridx]=Aptr[i];
	(*acsrcols)[csridx]=i-nsqr;
	lclctr++;	csridx++;	
      }
    }
    if((i-xdim)>=0)// diagonal with offset -n avlbl
    {
      if(fabs(Aptr[dim+i])>0.0f)
      {
	(*acsrvals)[csridx]=Aptr[dim+i];
	(*acsrcols)[csridx]=i-xdim;
	lclctr++;	csridx++;
      }
    }
    if((i-1)>=0)// diagonal with offset -1 avlbl
    {
      if(fabs(Aptr[2*dim+i])>0.0f)
      {
	(*acsrvals)[csridx]=Aptr[2*dim+i];
	(*acsrcols)[csridx]=i-1;
	lclctr++;	csridx++;
      }
    }
    if(fabs(Aptr[3*dim+i])>0.0f)
      {
	(*acsrvals)[csridx]=Aptr[3*dim+i];
	(*acsrcols)[csridx]=i;
	lclctr++;	csridx++;
      }
    if((i+1)<=dim-1)// diagonal with offset +1 avlbl
    {
      if(fabs(Aptr[2*dim+i+1])>0.0f)
      {
	(*acsrvals)[csridx]=Aptr[2*dim+i+1];
	(*acsrcols)[csridx]=i+1;
	lclctr++;	csridx++;
      }
    }
    if((i+xdim)<=dim-1)// diagonal with offset n avlbl
    {
      if(fabs(Aptr[dim+i+xdim])>0.0f)
      {
	(*acsrvals)[csridx]=Aptr[dim+i+xdim];
	(*acsrcols)[csridx]=i+xdim;
	lclctr++;	csridx++;
      }
    }
    if((i+nsqr)<=dim-1)// diagonal with offset n*n avlbl
    {
      if(fabs(Aptr[i+nsqr])>0.0f)
      {
	(*acsrvals)[csridx]=Aptr[i+nsqr];
	(*acsrcols)[csridx]=i+nsqr;
	lclctr++;	csridx++;	
      }
    }
    (*acsrrows)[i+1]=lclctr;
  }
}


int  call_paralution_dpcg(double *acsrvals, int *acsrrows, int *acsrcols, double *bin, double* xin, 
			  int max_iter, double tolerance, int nnzA, int dim, int *iter_count,
			  double *res_end, int xdim, int ydim, int zdim, int setlssd_in,
			  int defvex_perdirec_in, int lvst_offst_in, int phisize_in,
			  double* level_set, double *xout)
{
    struct timeval now;
  double tick, tack;
  int status, i;
  int *bubmap_ptr=NULL, maxbmap;
  //wrap A into paralution CSR object
  
  LocalVector<double> x;
  LocalVector<double> rhs;
  LocalMatrix<double> mat;
  mat.SetDataPtrCSR(&acsrrows, &acsrcols, &acsrvals, "A", nnzA, dim,dim);
  rhs.SetDataPtr(&bin, "rhs", dim);
  x.SetDataPtr(&xin, "x", dim);
  //wrap b and x into paralution vectors
  //setup DPCG solver
  // Linear Solver
  DPCG<LocalMatrix<double>, LocalVector<double>, double > ls;

    
#ifdef GPURUN  
  mat.MoveToAccelerator();
  x.MoveToAccelerator();
  rhs.MoveToAccelerator();
  ls.MoveToHost();
#else
  mat.MoveToHost();
  x.MoveToHost();
  rhs.MoveToHost();
  ls.MoveToHost();
#endif
  
gettimeofday(&now, NULL);
tick = now.tv_sec*1000000.0+(now.tv_usec);
  
#ifdef BUBFLO
  ls.SetNVectors_eachdirec(defvex_perdirec_in, defvex_perdirec_in, defvex_perdirec_in);
  ls.Set_alldims(xdim, ydim, zdim);
  //ls.Setxdim(xdim);
  ls.Setlvst_offst(lvst_offst_in);
  ls.SetNVectors(defvex_perdirec_in);
  ls.SetZlssd(setlssd_in);
  if(setlssd_in){
      LocalVector<int> bubmap;
      bubmap.Allocate("bubmap",mat.get_nrow());
      bubmap.LeaveDataPtr(&bubmap_ptr);
      bubmap_create(level_set, bubmap_ptr, xdim, ydim, zdim, mat.get_nrow(), &maxbmap, lvst_offst_in);
  }
#endif  

  ls.SetOperator(mat);
  ls.Init(0.0, tolerance, 1e8, max_iter);

#ifdef BUBFLO
  ls.MakeZ_CSR(); // requires xdim_ and novecni_ and zlssd_ to be set
  if(setlssd_in)
    ls.MakeZLSSD(bubmap_ptr, maxbmap); // bubmap must be ready and maxbmap available
#endif 
 
  ls.Build();
  mat.ConvertToDIA();
  
gettimeofday(&now, NULL);
tack = now.tv_sec*1000000.0+(now.tv_usec);
std::cout << "Building:" << (tack-tick)/1000000 << " sec" << std::endl;
  
//   ls.Verbose(2);

  mat.info();

gettimeofday(&now, NULL);
tick = now.tv_sec*1000000.0+(now.tv_usec);
  //call DPCG solver
//   x.info();	rhs.info();
  ls.Solve(rhs, &x);

gettimeofday(&now, NULL);
tack = now.tv_sec*1000000.0+(now.tv_usec);
std::cout << "Solver execution:" << (tack-tick)/1000000 << " sec" << std::endl;

  status=ls.GetSolverStatus();
  *iter_count=ls.GetIterationCount();
  *res_end=ls.GetCurrentResidual();
  //cout<<"\n iteration count is "<<*iter_count<<" and residual is" <<std::scientific<<*res_end<<endl;
  
gettimeofday(&now, NULL);
tick = now.tv_sec*1000000.0+(now.tv_usec);
  
  x.MoveToHost();
  mat.MoveToHost();
  rhs.MoveToHost();
  mat.ConvertToCSR();
  free(acsrcols);	free(acsrrows);	free(acsrvals);
  acsrcols=NULL;	acsrrows=NULL;	acsrvals=NULL;
  mat.LeaveDataPtrCSR(&acsrrows, &acsrcols, &acsrvals);

  free(xin);	free(bin);	xin=NULL;	bin=NULL;
  rhs.LeaveDataPtr(&bin);
  x.LeaveDataPtr(&xin);
  memcpy(xout, xin,sizeof(double)*dim);
  ls.Clear();
  
gettimeofday(&now, NULL);
tack = now.tv_sec*1000000.0+(now.tv_usec);
std::cout << "Cleanup:" << (tack-tick)/1000000 << " sec" << std::endl;
  
  
  return status;
}


