/*
  Library for r/w of specific data structures with hdf5 interface
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h> //Needed by genArrManip
#include <time.h> //For timestamping

/* Macros */
#include "macros.h"
/* For array manipulation */
#include "genArrSorting.h"
#include "genArrManip.h"
/* File I/O functions */
#include "hdf5Interface.h" //main interface
#include "hdf5rwDataStructs.h" //for relevant data types we define

/* Define Functions */
void writeWignerDelArr(double ***delArr_ptr, float L);
int checkWignerDelArrExists(float L);
double*** readWignerDelArr(float L);
void writeCw3jStruct(cplCw3jStruct *cplCw3jStruct_ptr);
int checkCw3jStructExists(double L, double dl);
cplCw3jStruct* readcplCw3jStruct(double L, double dl);
void getCurTmStmp(char *tstampBuffer);
int writeS2ScalarAdvectionSimData(int L, double h, long unsigned int Nst,
				  long unsigned int eqNum,
				  double t0,
				  double *tArr_ptr,
				  complex double ***usalm_ptr,
				  complex double **a2dArr_ptr);


/*
  To Do:
   C.f. speed when storing by indexing params not by hash
   (This will be irrelevant for scalar advection calc)
*/


/*
  Write Wigner Del array. First reshape to 1d array then write using
  hdf5Interface
*/
void writeWignerDelArr(double ***delArr_ptr, float L){  
  int Lin = floor(L);

  /* Reshape to 1D array for writing */
  double *delArr1D_ptr;
  delArr1D_ptr = reshapeTriTo1DArr(delArr_ptr, L);

  /* Write using hdf5Interface */
  const char wignerDelAbsPath[] = "/calculationData/wigner/delFuncs";
  int dblArr_sz = ceil(1.0/6.0*(Lin+1)*(Lin+2)*(Lin+3));
  char dataName[50];
  sprintf(dataName, "(%.1f)", L);
  writeDouble1DArr(delArr1D_ptr, dblArr_sz, wignerDelAbsPath, dataName); 

  /* Free/release resources */
  kill1DArr(delArr1D_ptr);
}

/*
  Check if a precalculated Wigner Del array exists.
  ret 1 if yes, 0 if no
*/
int checkWignerDelArrExists(float L){
  char pathTmp[] = "/calculationData/wigner/delFuncs/";
  char dataName[50];
  char strTmp[140];
  
  sprintf(strTmp, "%s", pathTmp);
  sprintf(dataName, "(%.1f)", L);  
  strcat(strTmp, dataName);
  strcat(strTmp, "\0");
  
  const char *pathChk_ptr = &strTmp[0]; //function takes const
  if(checkDatasetExists(pathChk_ptr)){
    return 1;
  }  
  return 0;
}


/*
  Read Wigner Del from file and convert to tri arr
*/
double*** readWignerDelArr(float L){
  double ***delArr_ptr;
  double *delArr1D_ptr;

  char strPathTmp[] = "/calculationData/wigner/delFuncs";
  const char *absolutePath_ptr = &strPathTmp[0];

  char dataName[50];
  sprintf(dataName, "(%.1f)", L);
  const char *dataName_ptr = &dataName[0];
  delArr1D_ptr = readDouble1DArr(absolutePath_ptr, dataName_ptr);

  /* Reshape to tri array */
  delArr_ptr = reshape1DtoTriArr(delArr1D_ptr, L);
  
  /* Free resources */
  kill1DArr(delArr1D_ptr);

  return delArr_ptr;
}


/*
  Write coupled cArr/wArr struct as arrays.
*/
void writeCw3jStruct(cplCw3jStruct *cplCw3jStruct_ptr){

  const char w3jAbsPath[] = "/calculationData/wigner/w3jFuncs";
  char dataName[50];
  
  //Write ullArr (hashRegge generated indices)
  sprintf(dataName, "c(%.1f,%.1f)",
	  cplCw3jStruct_ptr->L, cplCw3jStruct_ptr->dl);
  writeUllint1DArr(cplCw3jStruct_ptr->ullArr_ptr,
		   cplCw3jStruct_ptr->sz,
		   w3jAbsPath, dataName);
  
  //Write w3j symbols
  sprintf(dataName, "w(%.1f,%.1f)",
	  cplCw3jStruct_ptr->L, cplCw3jStruct_ptr->dl);
  writeDouble1DArr(cplCw3jStruct_ptr->dblArr_ptr,
		   cplCw3jStruct_ptr->sz,
		   w3jAbsPath, dataName);

  //Write dimension information
  sprintf(dataName, "i(%.1f,%.1f)",
	  cplCw3jStruct_ptr->L, cplCw3jStruct_ptr->dl);
  ullint dimWrite = cplCw3jStruct_ptr->sz;
  writeUllint1DArr(&dimWrite,
		   1,
		   w3jAbsPath, dataName);

  /* Free/release resources */
  //Kill input?
  //killCplCw3jStruct(cplCw3jStruct_ptr);
}


/*
  Check if a precalculated Cw3j structure exists.
  (The check is based on if the i(L, dl) dataset exists)
  ret 1 if yes, 0 if no
*/
int checkCw3jStructExists(double L, double dl){
  char pathTmp[] = "/calculationData/wigner/w3jFuncs/";
  char dataName[50];
  char strTmp[140];
  
  sprintf(strTmp, "%s", pathTmp);
  sprintf(dataName, "i(%.1f,%.1f)", L, dl);  
  strcat(strTmp, dataName);
  strcat(strTmp, "\0");
  
  const char *pathChk_ptr = &strTmp[0]; //function takes const
  if(checkDatasetExists(pathChk_ptr)){
    return 1;
  }  
  return 0;
}


/*
  Construct cplCw3jStruct from saved data
*/
cplCw3jStruct* readcplCw3jStruct(double L, double dl){
  cplCw3jStruct *cplCw3jStruct_ptr;
  ullint *ullArr_ptr, *sz_ptr;
  double *dblArr_ptr;

  char strPathTmp[] = "/calculationData/wigner/w3jFuncs";
  const char *absolutePath_ptr = &strPathTmp[0];

  char dataName[50];
  const char *dataName_ptr = &dataName[0];

  //Read reggeHash index
  sprintf(dataName, "c(%.1f,%.1f)", L, dl);
  ullArr_ptr = readUllint1DArr(absolutePath_ptr, dataName_ptr);

  //Read w3jSymbols
  sprintf(dataName, "w(%.1f,%.1f)", L, dl);
  dblArr_ptr = readDouble1DArr(absolutePath_ptr, dataName_ptr);
  
  //Read dim info
  sprintf(dataName, "i(%.1f,%.1f)", L, dl);
  sz_ptr = readUllint1DArr(absolutePath_ptr, dataName_ptr);
  
  //Construct from input data
  cplCw3jStruct_ptr = constrCplCw3jStruct(*sz_ptr, ullArr_ptr, dblArr_ptr);
  cplCw3jStruct_ptr->L = L;
  cplCw3jStruct_ptr->dl = dl;
  
  //Free resources
  free(sz_ptr);
  return cplCw3jStruct_ptr;
}


/*
  Get current time-stamp (for labelling output data)
  requires an input buffer
*/
void getCurTmStmp(char *tstampBuffer){
  //char tstampBuffer[50];
  time_t rawtime;
  struct tm* timeinfo;

  time(&rawtime);
  timeinfo = localtime(&rawtime);
  
  sprintf(tstampBuffer, "(%d,%02d,%02d,%02d,%02d,%02d)",
	  1900+(timeinfo->tm_year),
	  (timeinfo->tm_mon)+1,
	  timeinfo->tm_mday,
	  timeinfo->tm_hour,
	  timeinfo->tm_min, 
	  timeinfo->tm_sec);
}


/*
  Write S2 Scalar Advection simulation data
*/
int writeS2ScalarAdvectionSimData(int L, double h, long unsigned int Nst,
				  long unsigned int eqNum,
				  double t0,
				  double *tArr_ptr,
				  complex double ***usalm_ptr,
				  complex double **a2dArr_ptr){

  /* Make directory for data with timestamp */
  char strDirBase[]="/simulationData/S2advection/scalar/";
  char tstampBuffer[200];
  char dirTmp[500];
  
  getCurTmStmp(&tstampBuffer[0]);
  sprintf(dirTmp, "%s%s", strDirBase, tstampBuffer);
  
  printf("\nOutput in: %s\n", dirTmp);
  writeSimFolder(dirTmp);


  /* Write using hdf5Interface */

  //First write usalm (vector field we advect along)
  char udTmp[500];
  sprintf(udTmp, "%s%s/ualm", strDirBase, tstampBuffer);
  writeSimFolder(udTmp); //write directory for ualm coeffs

  long unsigned int i;
  for(i=0; i<3; i++){
    double *reulm1d_ptr = reshapealmTo1D(usalm_ptr[i], L, 0);
    double *imulm1d_ptr = reshapealmTo1D(usalm_ptr[i], L, 1);
    char nameTmpA[100];
    char nameTmpB[100];
    sprintf(nameTmpA, "ualm(%d)-Re", ((int)i)-1);
    sprintf(nameTmpB, "ualm(%d)-Im", ((int)i)-1);
     
    writeDouble1DArr(reulm1d_ptr, pow(L+1,2), udTmp, nameTmpA); 
    writeDouble1DArr(imulm1d_ptr, pow(L+1,2), udTmp, nameTmpB); 
  
    free(reulm1d_ptr);
    free(imulm1d_ptr);
  }

  //Scalar field coefficients
  char almTmp[500];
  sprintf(almTmp, "%s%s/alm", strDirBase, tstampBuffer);
  writeSimFolder(almTmp); //write directory for ualm coeffs

  writeComplexDouble2DArr(a2dArr_ptr,
			  Nst+1,eqNum,
			  almTmp,
			  "alm(0)");


  //Step-length
  writeDouble1DArr(&h, 1, dirTmp, "h(RK_Step)"); 
  //Initial time
  writeDouble1DArr(&t0, 1, dirTmp, "t0(RK_Ini)"); 
  //Time-steps evaluated
  writeDouble1DArr(tArr_ptr, Nst+1, dirTmp, "t(RK_t)"); 

  //Number of steps taken
  ullint NstPass = (ullint) Nst;
  writeUllint1DArr(&NstPass, 1, dirTmp, "Nst(RK_StepNum)");
  
  //Total number of equations
  ullint eqNumPass = (ullint) eqNum;
  writeUllint1DArr(&eqNumPass, 1, dirTmp, "eqNum(salm)");
  //Band-limit
  ullint LPass = (ullint) L;
  writeUllint1DArr(&LPass, 1, dirTmp, "L");

   
  return 1;
}
