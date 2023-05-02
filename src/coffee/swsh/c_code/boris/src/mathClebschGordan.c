/*
  Library for computation of Clebsch-Gordan coefficients.
  Relation to w3j is used.
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h> //Needed by genArrManip

/* Macros */
#include "macros.h"
/* Generic functions */
#include "genArrManip.h" //array manipulation / killing
#include "genArrSorting.h" //Array sorting
/* Mathematical functions */
#include "mathWigner3j.h"
#include "mathClebschGordan.h"

/* Local function prototypes */
double clebschGordanLookup(double *lArr_ptr, double *mArr_ptr, 
			   cplCw3jStruct *cplCw3jStruct_ptr);
double* decompspSphHarmProd(double *sphA, double *sphB,
			    cplCw3jStruct *cplCw3jStruct_ptr);
/*
  Get Clebsch-Gordan coefficient from w3j
*/
double clebschGordanLookup(double *lArr_ptr, double *mArr_ptr, 
			   cplCw3jStruct *cplCw3jStruct_ptr){

  double mArrW[] = {mArr_ptr[0], mArr_ptr[1], -mArr_ptr[2]};
  double w3jVal = getW3jPrecalc(lArr_ptr, &mArrW[0], cplCw3jStruct_ptr);

  double ph = pow(-1, -lArr_ptr[0]+lArr_ptr[1]-mArr_ptr[2]);
  
  return ph*sqrt(2*lArr_ptr[2]+1)*w3jVal;
}

/*
  sphA={sA,lA,mA}, sphB={sB,lB,mB}
  return dbl [arrSize, lmin, cg1, cg2, ...]
*/
double* decompspSphHarmProd(double *sphA, double *sphB,
			    cplCw3jStruct *cplCw3jStruct_ptr){
  int i; //Iterator(s)

  //Find l index summation range
  double lMinS[] = {fabs(sphA[1]-sphB[1]), 
		    fabs(-(sphA[0]+sphB[0])),
		    fabs(sphA[2]+sphB[2])};

  double lMinTmp = lMinS[0];
  double lMaxTmp = sphA[1]+sphB[1];

  for(i=1; i<3; i++){
    if(lMinS[i]>lMinTmp)
      lMinTmp = lMinS[i];
  }
  int lNum = (int) (lMaxTmp-lMinTmp)+1;
  
  //Overall front factor
  double fF = sqrt((2*sphA[1]+1)*(2*sphB[1]+1)/(4*PI));

  double* decompArr = (double*) malloc((2+lNum)*sizeof(double));
  decompArr[0] = lNum; //Ensure that array size is returned

  int j=1;
  //  printf("->vals: sz=%.1f ", decompArr[0]);
  for(i=0; i<lNum; i++){
    double lArr[] = {sphA[1], sphB[1], lMinTmp+i};
    double mArrA[] = {-sphA[0], -sphB[0], -(sphA[0]+sphB[0])};
    double mArrB[] = {sphA[2], sphB[2], sphA[2]+sphB[2]};
    /*  
    printf("%.1f, %.1f, %.1f\n", lArr[0], lArr[1], lArr[2]);
    printf("%.1f, %.1f, %.1f\n", mArrA[0], mArrA[1], mArrA[2]);
    printf("%.1f, %.1f, %.1f\n\n", mArrB[0], mArrB[1], mArrB[2]);
    */
    if(lArr[2]-0.1<=cplCw3jStruct_ptr->L){
      if((lArr[0]>fabs(mArrA[0])-0.1) && (lArr[0]>fabs(mArrB[0])-0.1) && \
	 (lArr[1]>fabs(mArrA[1])-0.1) && (lArr[1]>fabs(mArrB[1])-0.1) && \
	 (lArr[2]>fabs(mArrA[2])-0.1) && (lArr[2]>fabs(mArrB[2])-0.1)){
	j++;

	decompArr[j] = fF*1.0/sqrt(2*(lMinTmp+i)+1)*			\
	  clebschGordanLookup(&lArr[0], &mArrA[0], cplCw3jStruct_ptr)*	\
	  clebschGordanLookup(&lArr[0], &mArrB[0], cplCw3jStruct_ptr);
	
	//printf("%.1f, %d\n", lMinTmp, i);
      }
    }
  }
  decompArr = (double*) realloc(decompArr, (j+1)*sizeof(double));
  decompArr[0] = j+1; //Array size
  decompArr[1] = lMinTmp; //minimum l value returned

  return &decompArr[0];
}
