//#include "../headers/array.h"

#include <stdio.h>      /* printf, scanf, NULL */
#include <stdlib.h>     /* calloc, exit, free */
#include<iostream>

/*****************************************************************************
Prototype:    Array2<double>    double_Matrix(int row, int col) 
Description:     To allocate a  matrix of double type. 
Input value: 
    int row --- the row number of the matrix 
    int col ---  column number of the matrix 
Return value:     a pointer to the matrix 
*****************************************************************************/
Array2<double> double_Matrix(int row, int column){
    	int i;
    	Array2<double> m;



    	m=(Array2<double> ) calloc((unsigned) row,sizeof(Array1<double> ));

    	for(i=0;i<row;i++){

        	m[i]=(Array1<double> ) calloc((unsigned) column,sizeof(double));
    		if(m[i]==NULL){
			fprintf(stderr,"Allocation failed in double_Matrix with \n \
		        arguments row %i column %i \n",row,column);
			exit(1);
    		}
	}
    	if(m==NULL){
		fprintf(stderr,"Allocation failed in double_Matrix with \n \
			        arguments row %i column %i \n",row,column);
		exit(1);
    	}
   
    return m;
}



/*****************************************************************************
Prototype:    Array1<double> double_Vector(int nh) 
Description:    To allocate a vector of double type. 
Input value: 
    int nh --- length of the vector 
Return value:    pointer to the vector 
*****************************************************************************/
Array1<double> double_Vector( int nh){
    	Array1<double> v;


     	v=(Array1<double> )calloc((unsigned) nh,sizeof(double)); 
   	
	if(v==NULL){
		fprintf(stderr,"Allocation failed in double_Vector with\n  \
		        	argument nh %i \n",nh);
		exit(1);
    	}

    	return v;
}




/*****************************************************************************
Prototype:  void  free_double_Matrix( Array2<double> m, int row, int column) 
Description:     To free a matrix of double type. 
Input value: 
    Array2<double> m --- pointer to a matrix of double type. 
    int row --- row number of the matrix 
    int column --- column number of the matrix 
Return value:     none 
Note: no reference is needed to the last dimension
*****************************************************************************/
void free_double_Matrix(Array2<double> m, int row)
{
    	int i;

    	for(i=row-1;i>=0;i--){ 
		free((char*) m[i]);
	}

        free((char*) m);
}


/****************************************************************************
Prototype:  void  free_double_Vector(Array1<double> v, int nh) 
Description:     To free a vector of  double type. 
Input value: 
    Array1<double>  v --- pointer to a vector of double type 
    int nh --- length of the vector 
Return value:     none 
****************************************************************************/
void free_double_Vector(Array1<double> v)
{
    	free((char *) v);
}
