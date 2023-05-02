/*
  Library for calculation of Wigner 3j symbols via hybrid recursion method
  THIS CODE NEEDS TO BE REFACTORED/CLEANED UP:
  void calcCoeff(double *jmVals_ptr, double *retW3jArr_ptr); will be retained
  as an external function that may be called, the other internal stuff may be
  changed..
*/

#include <stdio.h>
#include <stdlib.h> //for calloc
#include <math.h>
#include <float.h> //for precision constants (DBL_EPSILON)
#include <string.h> //for memcpy

/* Mathematical functions */
#include "mathWigner3jRecursion.h"

// Macros
#ifndef SIZEOF_DBL
#define SIZEOF_DBL (sizeof(double))
#endif

#ifndef SIZEOF_LONG
#define SIZEOF_LONG (sizeof(long))
#endif

#ifndef UNDERFLOW_MULT
#define UNDERFLOW_MULT  1000
#endif

//If a ratio exceeds this number, revert to linear
#ifndef FRAC_REVERT 
#define FRAC_REVERT 1000
#endif

/* Function prototypes */

void calcCoeff(double *ptr_jmVals, double *retw3jArr); //Prototype of the main algo func.

double coeffA(const double jvar, const double *ptr_jmVals); //Prototype of the recursion factors A,B
double coeffB(const double jvar, const double *ptr_jmVals);

double leftNL(const double Jval, const double hm1, const double *ptr_jmVals, int *ptr_isSingular); //Prototype of the Nonlinear recursion functions
double rightNL(const double Jval, const double gp1, const double *ptr_jmVals, int *ptr_isSingular);

double leftLi(const double Jval, const double fm1, const double fe1, const double *ptr_jmVals); //Prototype of the Linear recrusion functions
double rightLi(const double Jval, const double fe1, const double fp1, const double *ptr_jmVals);

int phaseConv(double *ptr_jmVals); //Prototype for the phase convention

//Hybrid scheme prototypes
int leftNLfun(int j1num, double j1min, double j1max, int *ptr_classicalLeak, double *jmVals, double *ptr_arrRetLeftNL); //Prototype of the Nonlinear (left->right) recursion algorithm
int rightNLfun(int j1num, double j1min, double j1max, int *ptr_classicalLeak, double *jmVals, double *ptr_arrRetRightNL, double zZrnd, int *ptr_sEp); //Prototype of the Nonlinear (left<-right) recursion algorithm
void patchLinear(int j1num, double j1min, double j1max, double zZrnd, int j1L, int j1R, int *ptr_classicalLeak, double *ptr_jmVals, double *ptr_arrRet, double *ptr_arrRetLeftNL, double *ptr_arrRetRightNL); //Prototype of patching linear scheme

//Three term fully linear scheme
void fullLinear(int j1num, double j1min, double j1max, double zIni, double zZrnd, double *ptr_jmVals, double *ptr_arrRet); //Prototype of fully linear scheme

void normalize(double j1min, double j1max, double j1num, int phase, int signArrMax, double *ptr_arrRet); //Normalise array


void calcCoeff(double *ptr_jmVals, double *retw3jArr){
  double zIni=0.1; //In case of singularity, use the following initialisation
  double zZrnd=DBL_EPSILON*1000; //Machine round to zero
  
  int sEp=1; //record the signs in a seperate variable
  int *ptr_sEp=&sEp;

  double j1max=ptr_jmVals[0]+ptr_jmVals[1]; //j2+j3
  double j1min=fmax(fabs(ptr_jmVals[0]-ptr_jmVals[1]), fabs(ptr_jmVals[2]));//max(|j2-j3|,|m1|)
  
  int j1num=(int)(j1max-j1min+1); //Total number of valid j1 entries for this set of indices
  //printf("j1num=%d\n",j1num);
  //if(j1num>1){
  int cL=j1num<=3; //If less than 3 terms, revert to 3-term linear recursion
  int *ptr_classicalLeak=&cL;
  *ptr_classicalLeak=1; //force test of linear algo
  //*ptr_classicalLeak=0; //force test of hybrid algo
  if(*ptr_classicalLeak==0){ //now, if recursion wasnt singular, attempt right to left

    //printf("\n====]Attempting hybrid calculation..\n");
    //printf("=>Nonlinear [Left->Right]\n");
    double *ptr_arrRetLeftNL=calloc(j1num, SIZEOF_DBL);
    
    int j1L=0;
    
    if(fabs(coeffB(j1min,ptr_jmVals))>zZrnd){
      j1L=leftNLfun(j1num, j1min, j1max, ptr_classicalLeak, ptr_jmVals, ptr_arrRetLeftNL);
    } else {
      *ptr_classicalLeak=1;
    }
    
    if(*ptr_classicalLeak==0){ //now, if recursion wasnt singular, attempt right to left
      int j1R;
      //double *arrRetRightNL=calloc(j1num, SIZEOF_DBL);
      //double *ptr_arrRetRightNL=&arrRetRightNL[0];
      double *ptr_arrRetRightNL=calloc(j1num, SIZEOF_DBL);
      //printf("=>Nonlinear [Left<-Right]\n");

      
      if(fabs(coeffB(j1max,ptr_jmVals))>zZrnd){
	j1R=rightNLfun(j1num, j1min, j1max, ptr_classicalLeak, ptr_jmVals, ptr_arrRetRightNL, zZrnd, ptr_sEp);
      } else {
	*ptr_classicalLeak=1;
      }
      
      
      if(*ptr_classicalLeak==0){//if still no leak, attempt to patch solutions together
	//double arrRetpl[j1num];//The array used for holding the w3j (nontrivially 0 positions are retained)
	//double *ptr_arrRetpl=&arrRetpl[0];
	//printf("=>patchLinear\n");
	patchLinear(j1num, j1min, j1max, zZrnd, j1L, j1R, ptr_classicalLeak, ptr_jmVals, retw3jArr, ptr_arrRetLeftNL, ptr_arrRetRightNL);
	//free(arrRetRightNL);
	//free(ptr_arrRetRightNL);
	
	if(*ptr_classicalLeak==0){//if patching succeeds, normalize and apply phase
	  //printf("=>normalization and phase application\n");
	  int phase=phaseConv(ptr_jmVals); //the phase convention
	  //printf("ptr_sEp=%d\n",*ptr_sEp); //sign compensation for numerical underflow

	  //printf("Normalized::\n");
	  normalize(j1min, j1max, j1num, phase, *ptr_sEp, retw3jArr);
	  
	  //we now check if we have inf/nan
	  int i;
	  for(i=0; i<j1num; i++){
	    if(isnan(retw3jArr[i])||isinf(retw3jArr[i])){
	      *ptr_classicalLeak=1;
	      //printf("\n ====]Reverting to linear\n"); //Perform the entire recursion linearly i.e., classicalLeak=true
	    }
	  }
	  
	  //int i;
	  //for(i=0; i<j1num; i++){
	  //printf("ptr_arrRetpl:=%1.13f\n",ptr_arrRetpl[i]);
	  //}
	   
	  //} else {
	  //printf("\n ====]Reverting to linear\n"); //Perform the entire recursion linearly i.e., classicalLeak=true
	}
      }
      
      //free(arrRetRightNL);
      free(ptr_arrRetRightNL);
    }
    
    free(ptr_arrRetLeftNL);
    //free(arrRetLeftNL);
  }
  
  
  //*ptr_classicalLeak=1; //force test of linear algo
  if(*ptr_classicalLeak==1){//Revert to full linear recursion
    
    
    //double arrRet[j1num];//The array used for holding the w3j (nontrivially 0 positions are retained)
    //double *ptr_arrRet=&arrRet[0];
    
    fullLinear(j1num, j1min, j1max, zIni, zZrnd, ptr_jmVals, retw3jArr);

    int phase=phaseConv(ptr_jmVals); //the phase convention
    
    //printf("Normalized::\n");
    normalize(j1min, j1max, j1num, phase, pow(-1,(retw3jArr[j1num-1]<0)), retw3jArr);
    /*
      int i;
      for(i=0; i<j1num; i++){
      printf("ptr_arrRet:=%1.13f\n",ptr_arrRet[i]);
      }
    */
    
  }

  /*
    if(cmpNL==1){
    printf("====================\n");
    int i;
    for(i=0; i<j1num; i++){
    printf("w==r:=%d wArrTmp[i]:=%1.13f   retw3jArr[i]=%1.13f\n",(abs(wArrTmp[i]-retw3jArr[i])<zZrnd),wArrTmp[i],retw3jArr[i]);
    }
    printf("====================\n");
    }
  */
  /*
    int i;
    for(i=0; i<j1num; i++){
    if((retw3jArr[i]<-0.018486845)&(retw3jArr[i]>-0.018486847)){
    printf("====================\n");
    printf("====================\n");
    printf("j1min:=%f j1max:=%f i:=%d retw3jArr[i]:=%1.16f\n",j1min,j1max, i, retw3jArr[i]);
    printf("::%f,%f,%f,%f,%f::\n",ptr_jmVals[0],ptr_jmVals[1],ptr_jmVals[2],ptr_jmVals[3],ptr_jmVals[4]);
    printf("====================\n");
    printf("j1min:=%f j1max:=%f i:=%d retw3jArr[i]:=%1.16f\n",j1min,j1max, i+1, retw3jArr[i+1]);
    printf("::%f,%f,%f,%f,%f::\n",ptr_jmVals[0],ptr_jmVals[1],ptr_jmVals[2],ptr_jmVals[3],ptr_jmVals[4]);
    printf("====================\n");
    printf("====================\n");
    }
    }

  */
}


//=================Define reused eqns
//The phase convention employed (-)^(j2-j3-m1)
int phaseConv(double *ptr_jmVals){
  //printf("::%f,%f,%f::\n",ptr_jmVals[0],ptr_jmVals[1],ptr_jmVals[2]);
  return pow(-1,ptr_jmVals[0]-ptr_jmVals[1]-ptr_jmVals[2]);
}
//Normalization (take input pointer to array, modify that inplace)
void normalize(double j1min, double j1max, double j1num, int phase, int signArrMax, double *ptr_arrRet){
  //double s=phase*sArrMax;
  double tmp=0;
  int i;
  for(i=(2*j1min+1); i<=(2*j1max+1); i=i+2){//calculate normalization constant
    tmp=tmp+i*pow(ptr_arrRet[(int) ((i-2*j1min-1)/2)],2);
  }
  
  double sApp=pow(-1, signArrMax!=phase);//Flip sign during normalization (adopt phase convention)
  //now apply normalization to input array inplace
  for(i=0; i<j1num; i++){
    ptr_arrRet[i]=sApp*ptr_arrRet[i]/sqrt(tmp);
  }

}



//==Recursion `factors' A(j1), B(j1)
//A(j1; j2, j3,m1)  (j1 is the swept index (variable), j2,j3,m1 are fixed)
double coeffA(const double jvar, const double *ptr_jmVals){
  //printf("coeffAdbg=%1.20f\n",(pow(jvar,2)-pow(ptr_jmVals[0]-ptr_jmVals[1],2))*((pow(ptr_jmVals[0]+ptr_jmVals[1]+1,2)-pow(jvar,2)))*(pow(jvar,2)-pow(ptr_jmVals[2],2)));
  //do not allow -ve arg. for surds ==> set to zero.
  double arg=(pow(jvar,2)-pow(ptr_jmVals[0]-ptr_jmVals[1],2))*((pow(ptr_jmVals[0]+ptr_jmVals[1]+1,2)-pow(jvar,2)))*(pow(jvar,2)-pow(ptr_jmVals[2],2));
  if(arg<0){
    //printf("wtf... arg=%1.16f\n",arg); ~-0.0000000000....
    return 0;
  } else {
    return sqrt(arg);
  }
  //return sqrt(pow(jvar,2)-pow(ptr_jmVals[0]-ptr_jmVals[1],2))*sqrt(pow(ptr_jmVals[0]+ptr_jmVals[1]+1,2)-pow(jvar,2))*sqrt(pow(jvar,2)-pow(ptr_jmVals[2],2));
}
//B(j1; j2,j3,m1,m2,m3)  (j1 is a variable, the rest are params)
double coeffB(const double jvar, const double *ptr_jmVals){
  return -(2*jvar+1)*(ptr_jmVals[0]*(ptr_jmVals[0]+1)*ptr_jmVals[2]-ptr_jmVals[1]*(ptr_jmVals[1]+1)*ptr_jmVals[2]-jvar*(jvar+1)*(ptr_jmVals[4]-ptr_jmVals[3]));
}

//==Nonlinear Recursion functions
//Two-term nonlinear, left to right
double leftNL(const double Jval, const double hm1, const double *ptr_jmVals, int *ptr_isSingular){
  double denom=coeffB(Jval,ptr_jmVals)+(Jval+1)*coeffA(Jval,ptr_jmVals)*hm1;
  
  if(fabs(denom)>0){
    double retVal=-Jval*coeffA(Jval+1,ptr_jmVals)/denom;
    
    if(fabs(retVal)>FRAC_REVERT){
      *ptr_isSingular=1;
      printf("retVal Large=%f denom=%f\n", retVal, denom);
      return 0.0f;
    }
     
    return retVal;
  } else {
    *ptr_isSingular=1;
    return 0.0f;
  }
   
  return -Jval*coeffA(Jval+1,ptr_jmVals)/denom;
}
//Two-term nonlinear, right to left
double rightNL(const double Jval, const double gp1, const double *ptr_jmVals, int *ptr_isSingular){
  double denom=coeffB(Jval, ptr_jmVals)+Jval*coeffA(Jval+1,ptr_jmVals)*gp1;
  
  if(fabs(denom)>0){
    double retVal=-(Jval+1)*coeffA(Jval, ptr_jmVals)/denom;
    
    if(fabs(retVal)>FRAC_REVERT){
      *ptr_isSingular=1;
      printf("retVal Large=%f denom=%f\n", retVal, denom);
      return 0.0f;
    }
     
    return retVal;
  } else {
    *ptr_isSingular=1;
    return 0.0f;
  }
  
  return -(Jval+1)*coeffA(Jval, ptr_jmVals)/denom;
}

//==Linear Recursion functions
//Three-term linear, left to right
double leftLi(const double Jval, const double fm1, const double fe1, const double *ptr_jmVals){
  return -(coeffB(Jval,ptr_jmVals)*fe1+(Jval+1)*coeffA(Jval,ptr_jmVals)*fm1)/(Jval*coeffA(Jval+1,ptr_jmVals));
}
//Three-term linear, right to left
double rightLi(const double Jval, const double fe1, const double fp1, const double *ptr_jmVals){
  //printf("coeffA(Jval+1,ptr_jmVals)=%f\n",coeffA(Jval+1,ptr_jmVals));
  //printf("rightLidbg=%f\n",-(Jval*coeffA(Jval+1,ptr_jmVals)*fp1+coeffB(Jval,ptr_jmVals)*fe1));
  return -(Jval*coeffA(Jval+1,ptr_jmVals)*fp1+coeffB(Jval,ptr_jmVals)*fe1)/((Jval+1)*coeffA(Jval,ptr_jmVals));
}

//=================Functions for recursion
//=====Nonlinear scheme (from left)
//==First calculate min index-->1st turning point using nonlinear recursion, this
//assumes all entries to the first turning point are non-zero (in the so
//called nonclassical regime)
int leftNLfun(int j1num, double j1min, double j1max, int *ptr_classicalLeak, double *ptr_jmVals, double *ptr_arrRetLeftNL){
  //double hj1leftnonlin[j1num]; //min index--->1st tp. index recursion
  //jInd=1; //for taking care of non-integer J (which may be integral/half-integral)
  int *ptr_isSingular=ptr_classicalLeak; //couple pointers (if denom singular in A/B coeffs => we need to use 3-term algo)
  
  if(coeffB(j1min, ptr_jmVals)!=0){
    ptr_arrRetLeftNL[0]=leftNL(j1min,0,ptr_jmVals,ptr_isSingular);
  }else{
    *ptr_classicalLeak=1;
  }
  //printf("ptr_classicalLeak=%d ptr_isSingular=%d\n",*ptr_classicalLeak,*ptr_isSingular);
  
  int jInd;
  if(*ptr_classicalLeak==0){
    int tp=0; //boolean for first turning point
    //int jInd=2; //next pos
    jInd=2;
    double J;
    double mCurVal;
    double mPreVal;
    
    while(tp!=1){//iterate nonlinear rec. to 1st tp.
      J=j1min+jInd-1; //current value of J

      
      ptr_arrRetLeftNL[jInd-1]=leftNL(J, ptr_arrRetLeftNL[jInd-2],ptr_jmVals,ptr_isSingular);
      mCurVal=fabs(ptr_arrRetLeftNL[jInd-1]);
      mPreVal=fabs(ptr_arrRetLeftNL[jInd-2]);
      

      if(*ptr_isSingular==1){//Kill if denom singular
	tp=1;
      } else if((mCurVal>1)||(mCurVal<mPreVal)){
	tp=1;
      } else if(J>j1max){
	*ptr_classicalLeak=1;
	tp=1;
      } else {
	jInd=jInd+1;
	if((jInd==j1num+1)||(jInd==0)){
	  *ptr_classicalLeak=1;
	  tp=1;
	}
      }
      
    }
    
  }
  
  //int k;
  //for(k=0; k<j1num; k++){
  //printf("aaptr_arrRetLeftNL:=%1.13f\n",ptr_arrRetLeftNL[k]);
  //}
  
  //we now convert ratios to individual terms
  if(*ptr_classicalLeak==0){
    int j1L=jInd; //index that exceeds unity
    int i;
    /*
      double *fj1leftnonlin = calloc(j1num,SIZEOF_DBL);
      double tmp=1;
      
      for(i=0; i<=(j1L-1); i++){
      tmp=tmp*ptr_arrRetLeftNL[j1L-i-1];
      fj1leftnonlin[j1L-i-1]=tmp;
      //ptr_arrRetLeftNL[j1L-i-1]=ptr_arrRetLeftNL[j1L-i-1]*ptr_arrRetLeftNL[j1L-i-2];
      }
    */
    for(i=1; i<=(j1L-1); i++){
      ptr_arrRetLeftNL[j1L-i-1]=ptr_arrRetLeftNL[j1L-i-1]*ptr_arrRetLeftNL[j1L-i];
      //fj1leftnonlin[j1L-i-1]=tmp;
      //ptr_arrRetLeftNL[j1L-i-1]=ptr_arrRetLeftNL[j1L-i-1]*ptr_arrRetLeftNL[j1L-i-2];
    }
    /*
      for(i=0; i<=(j1L-1); i++){
      ptr_arrRetLeftNL[i]=fj1leftnonlin[i];
      }
    */
    /*
      printf("ttttt\n");
      int k;
      for(k=0; k<j1num; k++){
      printf("ArrRet:=%1.13f ][%d\n",fj1leftnonlin[k],j1L);
      }
      printf("ttttt\n");
    */
    //remove garbage
    //free(fj1leftnonlin);
    
    return j1L;
  }
  
  return 0;
}

//=====Nonlinear scheme (from right)
int rightNLfun(int j1num, double j1min, double j1max, int *ptr_classicalLeak, double *ptr_jmVals, double *ptr_arrRetRightNL, double zZrnd, int *ptr_sEp){
  
  int *ptr_isSingular=ptr_classicalLeak; //couple pointers (if denom singular in A/B coeffs => we need to use 3-term algo)
  int jInd=j1num;
  //int sEp=1; //sign of maximum entry (to avoid signloss on numerical underflow)
  int j1R;

  if(coeffB(j1max, ptr_jmVals)!=0){
    ptr_arrRetRightNL[j1num-1]=rightNL(j1max,0,ptr_jmVals,ptr_isSingular);
    *ptr_sEp=pow(-1,!(ptr_arrRetRightNL[j1num-1]>0)); 
  }else{
    *ptr_classicalLeak=1;
  }
  //printf("ptr_classicalLeak(%d)\n",*ptr_classicalLeak);
  
  
  if(*ptr_classicalLeak==0){
    int tp=0; //2nd turning point not reached
    jInd=jInd-1; //towards centre
    double J;
    double mCurVal;
    double mPreVal;
    
    
    
    while(tp!=1){
      J=j1max+jInd-j1num;
      ptr_arrRetRightNL[jInd-1]=rightNL(J,ptr_arrRetRightNL[jInd],ptr_jmVals,ptr_isSingular);
      
      mCurVal=fabs(ptr_arrRetRightNL[jInd-1]);
      mPreVal=fabs(ptr_arrRetRightNL[jInd]);
      
      
      if(*ptr_isSingular==1){//Kill if denom singular
	*ptr_classicalLeak=1; //redundant.
	tp=1;
      } else if((mCurVal>1)||(mCurVal<mPreVal)){
	tp=1;
      } else if(J<j1min){
	*ptr_classicalLeak=1;
	tp=1;
      } else {
	jInd=jInd-1; //move to next index
	if((jInd==j1num+1)||(jInd==0)){
	  *ptr_classicalLeak=1;
	  tp=1;
	}
      }
    }
    
    if(*ptr_classicalLeak==0){
      j1R=jInd; //index exceeding unity
      
      //convert ratios to actual terms
      
      if(fabs(ptr_arrRetRightNL[j1num-1])< zZrnd*UNDERFLOW_MULT){//check for numerical underflow
	int k;
	//int tmpS=1;
	for(k=0; k<j1num; k++){
	  *ptr_sEp=*ptr_sEp*pow(-1,(ptr_arrRetRightNL[k]<0));
	}
	
	//printf("underFlow: *ptr_sEp=%d\n",*ptr_sEp);
	
      } else {//no overflow, proceed normally.
	int i;
	for(i=1; i<=(j1num-j1R); i++){
	  ptr_arrRetRightNL[j1R+i-1]=ptr_arrRetRightNL[j1R+i-1]*ptr_arrRetRightNL[j1R+i-2];
	}
	*ptr_sEp=pow(-1,(ptr_arrRetRightNL[j1num-1]<0)); //no underflow so take sign of last entry
      }
      return j1R;
    }
    
    
  } else {
    *ptr_classicalLeak=1;
  }
  //Now perform linear recursion over central segment.
  return 0;
}


//Connect the two nonlinear regimes, by using linear recursion over the central (classical) portion
void patchLinear(int j1num, double j1min, double j1max, double zZrnd, int j1L, int j1R, int *ptr_classicalLeak, double *ptr_jmVals, double *ptr_arrRet, double *ptr_arrRetLeftNL, double *ptr_arrRetRightNL){
  int doLinear; //we only perform linear recursion if required
  int j1midInd;
  int jInd;
  double J;
  
  
  if(*ptr_classicalLeak==0){
    if(j1L>j1R){//1st turning point > 2nd turning pt => inconsistency, revert to fully linear algo
      *ptr_classicalLeak=1;
      doLinear=0;
    } else if(j1L==j1R){//we do not need the linear recurrence (only 1 tp)
      j1midInd=j1L;
      doLinear=0;
    } else {
      j1midInd=floor((j1L+j1R)/2);
      doLinear=1;
    }
  }

  //=====Linear scheme [might have been cleaner to separate this and reuse code from fullLinear]
  if(doLinear==1){

    //double *arrRetLeftLiZ=calloc(j1num*1000000, SIZEOF_DBL);
    //printf("j1L=%d\n",j1L);
    
    //double *ptr_arrRetLeftLi=calloc(j1num, SIZEOF_DBL);
    double ptr_arrRetLeftLi[j1num];
    
    //double *arrRetLeftLi=calloc(j1num, SIZEOF_DBL);
    //double *ptr_arrRetLeftLi=&arrRetLeftLi[0];

    //initialise left right recursion
    ptr_arrRetLeftLi[j1L-2]=ptr_arrRetLeftNL[j1L-2];
    ptr_arrRetLeftLi[j1L-1]=ptr_arrRetLeftNL[j1L-1];
    
    
    
    for(jInd=j1L+1; jInd<=(j1midInd+1); jInd++){
      J=jInd+j1min-1;
      ptr_arrRetLeftLi[jInd-1]=leftLi(J-1,ptr_arrRetLeftLi[jInd-3],ptr_arrRetLeftLi[jInd-2], ptr_jmVals);
    }

    //initialise right left recursion
    //double *ptr_arrRetRightLi=calloc(j1num, SIZEOF_DBL);
    double ptr_arrRetRightLi[j1num];
    
    //double *arrRetRightLi=calloc(j1num, SIZEOF_DBL);
    //double *ptr_arrRetRightLi=&arrRetRightLi[0];
    ptr_arrRetRightLi[j1R-1]=ptr_arrRetRightNL[j1R-1];
    ptr_arrRetRightLi[j1R]=ptr_arrRetRightNL[j1R];
    
    int j1m1;
    for(jInd=j1midInd; jInd<=(j1R-1); jInd++){
      J=(j1max+(j1R-j1num)-1)+j1midInd-jInd;
      j1m1=j1R-jInd+j1midInd-1;
      ptr_arrRetRightLi[j1m1-1]=rightLi(J+1,ptr_arrRetRightLi[j1m1],ptr_arrRetRightLi[j1m1+1],ptr_jmVals);
    }
    
    if(*ptr_classicalLeak==0){
      //recursion must agree at midpoint!
      double cMul;
      if(fabs(ptr_arrRetRightLi[j1midInd-1])>zZrnd){
	cMul=ptr_arrRetRightLi[j1midInd-1]/ptr_arrRetLeftLi[j1midInd-1]; //centre-point nonzero
      } else {
	cMul=ptr_arrRetRightLi[j1midInd]/ptr_arrRetLeftLi[j1midInd]; //as centre-point zero, use next overlapping point
      }
      
      //combine and rescale the four sets  output=[cMul*leftNl lefLi rightLi rightNl]
      int i;
      for(i=1; i<=(j1L-1);i++){
	ptr_arrRet[i-1]=cMul*ptr_arrRetLeftNL[i-1];
      }
      for(i=j1L; i<=j1midInd; i++){
	ptr_arrRet[i-1]=cMul*ptr_arrRetLeftLi[i-1];
      }
      for(i=j1midInd+1; i<=j1R; i++){
	ptr_arrRet[i-1]=ptr_arrRetRightLi[i-1];
      }
      for(i=j1R+1; i<=j1num; i++){
	ptr_arrRet[i-1]=ptr_arrRetRightNL[i-1];
      }
    }
    //free(ptr_arrRetLeftLi);// free(arrRetLeftLi); 
    //free(ptr_arrRetRightLi);// free(arrRetRightLi); 
    
  } else {
    if(*ptr_classicalLeak==0){//otherwise tps match
      double cMul;
      if(fabs(ptr_arrRetRightNL[j1midInd-1])>zZrnd){
	cMul=ptr_arrRetRightNL[j1midInd-1]/ptr_arrRetLeftNL[j1midInd-1]; //centre-point nonzero
      } else {
	cMul=ptr_arrRetRightNL[j1midInd]/ptr_arrRetLeftNL[j1midInd]; //centre-point nonzero
      }
      //combine and rescale the sets
      int i;
      for(i=1; i<=(j1midInd-1);i++){
	ptr_arrRet[i-1]=cMul*ptr_arrRetLeftNL[i-1];
      }
      for(i=j1midInd; i<=j1num;i++){
	ptr_arrRet[i-1]=ptr_arrRetRightNL[i-1];
      }
      //free(ptr_arrRetRightNL);
      //free(ptr_arrRetLeftNL); 
    }
  }


}








//Perform the entire recursion linearly i.e., classicalLeak=true
void fullLinear(int j1num, double j1min, double j1max, double zIni, double zZrnd, double *ptr_jmVals, double *ptr_arrRet){
  int j1midInd=(int)round(((double) j1num)/2); //midpoint where recursions from both sides are equated
  //double fj1left[j1num];
  double *fj1left=calloc(j1num, SIZEOF_DBL);
  fj1left[0]=zIni; //initialisation
  int jInd;
  double J;
  
  if(j1num>1){
    if(j1min!=0){ //j1min!=0
      fj1left[1]=leftLi(j1min,0,fj1left[0], ptr_jmVals);
    }else{
      fj1left[1]=fj1left[0]*(ptr_jmVals[3]-ptr_jmVals[4])/(2*sqrt(ptr_jmVals[1]*(1+ptr_jmVals[1])));
    }
    //three term recursion
    for(jInd=3; jInd<=j1midInd+1; jInd++){
      J=jInd+j1min-2;
      fj1left[jInd-1]=leftLi(J,fj1left[jInd-3],fj1left[jInd-2], ptr_jmVals);
    }
  }
  /*
    int z;
    for(z=0; z<j1num; z++){
    printf("zfj1left:=%f\n",fj1left[z]);
    }
  */

  //double fj1right[j1num];
  double *fj1right=calloc(j1num, SIZEOF_DBL);
  fj1right[j1num-1]=zIni;
  

  
  
  if(j1num>1){
    fj1right[j1num-2]=rightLi(j1max, fj1right[j1num-1],0,ptr_jmVals);
    //printf("rightLi(j1max, fj1right[j1num-1],0,ptr_jmVals)=%f\n",rightLi(j1max, fj1right[j1num-1],0,ptr_jmVals));
    int j1m1;
    for(jInd=(int) fmax(1,(j1midInd-1)); jInd<=(j1num-2); jInd++){
      J=(j1max-2)+(j1midInd-1)-jInd+1; //current value of J (may be half-integral)
      j1m1=j1num-jInd+(j1midInd-1)-2; //current index
      fj1right[j1m1-1]=rightLi(J,fj1right[j1m1], fj1right[j1m1+1],ptr_jmVals);
    }
  }
  //printf("j1midInd(%d):j1num(%d):(fabs(fj1right[j1midInd-1])>zZrnd)=(%d)\n",j1midInd,j1num,fabs(fj1right[j1midInd-1])>zZrnd);




  //Recursion must agree at midpoint
  double cMul;
  if(fabs(fj1right[j1midInd-1])>zZrnd){
    cMul=fj1right[j1midInd-1]/fj1left[j1midInd-1];
  } else {
    
    if(j1midInd==1){
      cMul=fj1right[j1midInd]/fj1left[j1midInd];
    } else {
      cMul=fj1right[j1midInd-2]/fj1left[j1midInd-2];
    }
    
  }
  
  /*
    printf("cMul(%1.7f)\n", cMul);
    int d;
    for(d=0; d<j1num; d++){
    printf("fj1right:=%f\n",fj1right[d]);
    }
    
    for(d=0; d<j1num; d++){
    printf("fj1left:=%f\n",fj1left[d]);
    }
  */
  //Combine & rescale the four sets
  
  int i;
  for(i=0; i<=(j1midInd-1); i++){
    ptr_arrRet[i]=cMul*fj1left[i];
  }
  
  for(i=j1midInd; i<=(j1num-1); i++){
    ptr_arrRet[i]=fj1right[i];
  }
  
  free(fj1left); free(fj1right); //Remove garbage

}
