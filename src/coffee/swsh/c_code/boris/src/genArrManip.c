/*
  Library for manipulation of array structures
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>

/* Macros */
#include "macros.h"
/* Generic functions */
#include "genArrManip.h" //array manipulation / killing


/* Local function prototypes */
/* Reshape 3d tri to 1d arr [for wignerDel] */
double* reshapeTriTo1DArr(double ***delArr_ptr, float L);
/* Inverse of the above */
double*** reshape1DtoTriArr(double *delArr1D_ptr, float L);
/* Beat salm arrays into 1d [ret real / imag part] */
double* reshapealmTo1D(complex double **alm_ptr, int L, int reIm);

/* 
   Return [re/im=0/1][eqNum=0...L][Nst=0...TotSteps]
   Reshape time-evolution RK->(0)alm (use for scalar advection for instance) 
*/
double*** reshapeRKtEvTo0alm(complex double **yVec,
			     long unsigned int eqNum,
			     long unsigned int Nst);

/* Deallocation of data structures */
void killTriArr(double ***delArr_ptr, float L);
void kill1DArr(double *delArr_ptr);


/*
  Reshape the 3d triangular array to a 1d array for writing.
  For a given L, the total number of elements (white blocks fig 4.1) is
  1/6(L+1)(L+2)(L+3).
  We can malloc a contiguous block then populate each entry.
*/
double* reshapeTriTo1DArr(double ***delArr_ptr, float L){
  int l,m,n; //iterator(s)
  int Lin = (int) L;

  //Allocate 1D flattened arr
  double* delArr1D;
  double *delArr1D_ptr;
  delArr1D = (double*) malloc((int)(1+1.0/6.0*(Lin+1)*(Lin+2)*(Lin+3))*
			      sizeof(double));
  delArr1D_ptr = &delArr1D[0];

  int k=0; //index for 1d arr
  //Populate
  for(l=0; l<=Lin; l++){
    for(n=0; n<=l; n++){
      for(m=n; m<=l; m++){
	//delArr1D[l+n+(m-n)] = delArr_ptr[l][n][m-n];
	delArr1D[k] = delArr_ptr[l][n][m-n];
	//printf("%f, %d, %f\n", delArr1D[k], Lin, L);
	k++;
      }
    }
  }
  return delArr1D_ptr;
}

/*
  Inverse of the function reshapeTriTo1DarrWignerDel.
*/
double*** reshape1DtoTriArr(double *delArr1D_ptr, float L){
 
  int l,m,n; //iterator(s)
  int Lin = (int) L;

  //Allocate 3D triangular arr
  double*** delArr;
  double ***delArr_ptr;
  delArr = (double***) malloc((Lin+1)*sizeof(double*));
  delArr_ptr = &delArr[0];

  for(l=0; l<=Lin; l++){
    delArr[l] = (double**) malloc((l+1)*sizeof(double*));
    for(n=0; n<=l; n++){
      delArr[l][n] = (double*) malloc((l+1-n)*sizeof(double));      
    }
  }

  int k=0; //index for 1d arr
  //Populate
  for(l=0; l<=Lin; l++){
    for(n=0; n<=l; n++){
      for(m=n; m<=l; m++){
	//delArr_ptr[l][n][m-n] = delArr1D_ptr[l+n+(m-n)];
	delArr_ptr[l][n][m-n] = delArr1D_ptr[k];
	//printf("r:%f, %d, %f\n", delArr1D_ptr[k], Lin, L);
	k++;
      }
    }
  }

  return delArr_ptr;
}

/* Beat salm arrays into 1d [ret real / imag part] */
double* reshapealmTo1D(complex double **alm_ptr, int L, int reIm){
  int i = 0;
  int l, m;
  double* alm1d = (double*) malloc(pow(L+1,2)*sizeof(double));
  double *alm1d_ptr = &alm1d[0];


  for(l=0; l<=L; l++){
    for(m=0; m<=2*l; m++){
      if(reIm){
	alm1d_ptr[i] = creal(alm_ptr[l][m]);
      } else {
	alm1d_ptr[i] = cimag(alm_ptr[l][m]);
      } 
      i++;
    }
  }
  return alm1d_ptr;
}


/* 
   Return [re/im=0/1][eqNum=0...L][Nst=0...TotSteps]
   Reshape time-evolution RK->(0)alm (use for scalar advection for instance) 
*/
double*** reshapeRKtEvTo0alm(complex double **yVec,
			     long unsigned int eqNum,
			     long unsigned int Nst){
  return NULL;
}

/* Helper procedure for deallocation */
void killTriArr(double ***delArr_ptr, float L){
  int l, n; //iterators

  if(fmod(L,1)>0.1){
    /* Kill half-int */
    for(l=0; l<=(int) L; l++){
      for(n=0; n<=l; n++){
	free(delArr_ptr[l][n]);
      }
      free(delArr_ptr[l]);
    }
    free(delArr_ptr);

  } else {
    /* Kill int */
    for(l=0; l<=(int) L; l++){
      for(n=0; n<=l; n++){
	free(delArr_ptr[l][n]);
      }
      free(delArr_ptr[l]);
    }
    free(delArr_ptr);
  }
}

/* For symmetry... */
void kill1DArr(double *delArr_ptr){
  free(delArr_ptr);
}

