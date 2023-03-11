/*
  Library for computation of (half) integer Wigner-3j symbols
  -Useful for Clebsch Gordan decompositions
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* Settings */
#include "settings.h" //Some "magic numbers" stored here
/* Macros */
#include "macros.h"
/* Generic functions */
#include "genArrSorting.h" //Array sorting
/* Mathematical functions */
#include "mathWigner3j.h"
#include "mathWigner3jRecursion.h"

/* Define data structures */
typedef struct{ //w3j symbol info
  double *w3jArr_ptr; //calculated symbols
  double *jmVals_ptr; //calculation params
  ullint **cIndArr_ptr; //[cIndArr, oddFlipArr, J=j1+j2+j3 ] (ptr arr)
  double j1Min; //minimal valid j1 value
  ullint j1Num; //total number of symbols
} w3jStruct;

/* Local function prototypes */
ullint* hashRegge(double *jArr_ptr, double *mArr_ptr);
//Wrap recursion calculation, indexing function & nice datastructure
w3jStruct* calculateW3jFamily(double *jmArr_ptr, int toCalcHash);
void killW3jStruct(w3jStruct *w3jStruct_ptr, int hashCalced); //Free resources
void printW3jStruct(w3jStruct *w3jStruct_ptr, int hashCalced); //Debugging
ullint countW3jSymbols(double L, double dl); //Symbol counter
cplCw3jStruct* calcw3jSymbols(double L, double dl);
double calculateW3jVal(double *lArr_ptr, double *mArr_ptr);
//Extract position of WignerIndex from array
ullint findArrInd(ullint cInd, cplCw3jStruct* cplCw3jStruct_ptr);
double getW3jPrecalc(double *lArr_ptr, double *mArr_ptr,
		     cplCw3jStruct *cplCw3jStruct_ptr);

/*
To do:
 For large L numbers, dump recursion output to temporary files, then sort
*/

void swapPtr(double* jArr_ptr, double* mArr_ptr){
  double tmp_ptr = *jArr_ptr;
  *jArr_ptr = *mArr_ptr;
  *mArr_ptr = tmp_ptr;
}

/*
 * References: Rasch, J.; Yu, A.C.H. (2003) Efficient storage scheme for 
 * precalculated Wigner 3j, 6j and Gaunt coefficients. SIAM Journal of 
 * Scientific Computing. 25 (4), 1416-1428.
 */

/*
  J=j1+j2+j3;
  Regge symbol:
  R=[(-j1+j2+j3) (j1-j2+j3) (j1+j2-j3);
     (j1-m1)     (j2-m2)    (j3-m3);
     (j1+m1)     (j2+m2)    (j3+m3);];
  Permute to the form:

  R=[ S       L       (X+B-T);
      X       B       (S+L-T);
      (L+B-T) (S+X-T)   T      ]
  L>=X>=T>=B>=S and Each row/col sums to J

  The X,B,T are determined by requiring that R_{22}<R_{32} or if
  R_{22}=R_{32} then R_{23}<=R_{33}

  The Regge symbol allows for a computation of a 1d index, the total process
  may be viewed as a collisionless hash function.
*/
ullint* hashRegge(double *jArr_ptr, double *mArr_ptr){
  int i; //Iterator(s)
  int minInd, maxInd;
  int oddFlip=0;



  //The R symbol R[1,1]->R[0], R[1,2]->R[1],...., R[3,3]->R[8]
  double R[9];
  R[0]=-jArr_ptr[0]+jArr_ptr[1]+jArr_ptr[2]; //-j1+j2+j3
  R[1]=jArr_ptr[0]-jArr_ptr[1]+jArr_ptr[2]; //j1-j2+j3
  R[2]=jArr_ptr[0]+jArr_ptr[1]-jArr_ptr[2]; //j1+j2-j3
  R[3]=jArr_ptr[0]-mArr_ptr[0]; //j1-m1
  R[4]=jArr_ptr[1]-mArr_ptr[1]; //j2-m2
  R[5]=jArr_ptr[2]-mArr_ptr[2]; //j3-m3
  R[6]=jArr_ptr[0]+mArr_ptr[0]; //j1+m1
  R[7]=jArr_ptr[1]+mArr_ptr[1]; //j2+m2
  R[8]=jArr_ptr[2]+mArr_ptr[2]; //j3+m3

  
  //Find position of minimum index
  for(i=1, minInd=0; i<9; i++){
    if(R[minInd]>R[i]){
      minInd=i;
    }
  }
  
  //Perform a circular shift to move the smallest entry into R[1,1]->R[0]
  //Corresponds to circshift(R, [mrowSh, mcolSh]) in Matlab
  int mrowSh=1-((((int)floor(minInd/3))%3)+1); //(1-minRowInd)
  int mcolSh=1-((minInd%3)+1); //(1-minColInd)

  int tmp;
  tmp=-mrowSh;
  mrowSh=-mcolSh;
  mcolSh=tmp;

  /*
    We know in advance that |mrowSh|<=3 and |mcolSh|<=3, hence transpositions
    are sufficient to achieve a circshift in this case.
    Given [0 1 2]:
    We shift left by 1 (yielding 120) by transposing [ind1 ind2 0]=[1 0 2],
    and transposing [0 ind 2 ind 3]=[1 2 0]
    We shift left by 2 (yielding 201) by transposing [ind1 0 ind3]=[2 1 0],
    and transposing [0 ind 2 ind 3]=[2 0 1];
   */
  if(mrowSh<0)//remap {-2,-1}->{1,2}
    mrowSh=(3+mrowSh)%3;
  if((mrowSh>0)&(mrowSh<3)){
    //Only the first index transposition changes (shift by 1-->2) 
    //(shift by 2-->3) [-1 as indexing starts at 0]
    int tra=(1+1%mrowSh);
    for(i=0; i<3; i++){
      swapPtr(&R[3*i],&R[3*i+tra]); //1<-> [2/3] (switch columns)
      swapPtr(&R[3*i+1],&R[3*i+2]); //2<->3
    }
  }
  
  //Columns performed in a similar manner
  if(mcolSh<0)//remap {-2,-1}->{1,2}
    mcolSh=(3+mcolSh)%3;
  if((mcolSh>0)&(mcolSh<3)){
    int tra=((1%mcolSh)*3+3);
    for(i=0; i<3; i++){
      swapPtr(&R[i],&R[i+tra]);  //1<-> [2/3] (switch rows)
      swapPtr(&R[3+i],&R[6+i]);  //2<->3
    }
  }

  //Find position of maximum index
  for(i=1, maxInd=0; i<9; i++){
    if(R[maxInd]<R[i]){
      maxInd=i;
    }
  }


  int maxRow=(((int)floor(maxInd/3))%3)+1; //maxRowInd
  int maxCol=(maxInd%3)+1; //maxColInd

  //Move S and L into position
  if(maxCol==1){
    //Transpose
    swapPtr(&R[5],&R[7]);
    
    for(i=1; i<3; i++){
      swapPtr(&R[3*i],&R[i]);
    }
    
    if(maxRow==3){
      for(i=0; i<3; i++){
	swapPtr(&R[3*i+1],&R[3*i+2]); //First column unchanged
      }
      oddFlip=1;
    }
  } else if(maxCol==3){
    //Interchange last two columns
    for(i=0; i<3; i++){
      swapPtr(&R[3*i+1],&R[3*i+2]);
    }
    oddFlip=1;
  }

  
  if(R[7]<R[4]){
    for(i=0; i<3; i++){ //Interchange last two rows
      swapPtr(&R[6+i],&R[3+i]); //First row unchanged, swap last two rows
    }
    oddFlip=1-oddFlip;
  }else if((R[7]==R[4])&(R[8]<R[5])){
    
    for(i=0; i<3; i++){ //Interchange last two rows
      swapPtr(&R[6+i],&R[3+i]);
    }
    oddFlip=1-oddFlip;
  } 

  /*
    Now R has the correct form 
    ..Are the declarations here optimised out?
  */
  double S=R[0];
  double L=R[1];
  double X=R[3];
  double B=R[4];
  double T=R[8];
  
  double J=jArr_ptr[0]+jArr_ptr[1]+jArr_ptr[2];
  /*
  printf("(S,L,X,B,T,oddFlip*J)=(%f,%f,%f,%f,%f,%f)\n", 
	 S, L, X, B, T, oddFlip*J);

  printf("(%f, %f)\n", fmod(7,2), fmod(6,2));
  */
  ullint* retVals = (ullint*) 
    malloc(3*sizeof(ullint));

  retVals[0] = (ullint) L*(24+L*(50+L*(35+L*(10+L))))/120 + \
    X*(6+X*(11+X*(6+X)))/24+T*(2+T*(3+T))/6+B*(B+1)/2+S+1;
  retVals[1] = oddFlip;
  retVals[2] = (ullint) J; //J is non-negative int
  /*
  if(round(fmod(oddFlip*J,2))==0){
    //phase 1
    retVals[2]=1;
  } else {
    //phase -1
    retVals[2]=-1;
  }
  */

  return &retVals[0];
}


/*
  Calculate a particular family of W3j (sweeping j1),
  Return pointer arr:
  [toW3jSymbols, j1min, j1num]
 */
w3jStruct* calculateW3jFamily(double *jmArr_ptr, int toCalcHash){
  ullint i, j; //Iterator(s)

  //jmArr_ptr -> {j2,j3,m1,m2,m3}

  /* Construct data structure */
  w3jStruct *w3jStruct_ptr = (w3jStruct*) malloc(sizeof(w3jStruct));
  
  
  double j1Max = jmArr_ptr[0]+jmArr_ptr[1];
  ullint j1Num = (ullint)(j1Max-fmax(fabs(jmArr_ptr[0]-jmArr_ptr[1]),
				     fabs(jmArr_ptr[2]))+1);

  double* w3jArr = calloc(j1Num, sizeof(double));
  double *w3jArr_ptr = &w3jArr[0];
  calcCoeff(jmArr_ptr, w3jArr_ptr); //Perform recursion

  ullint** cIndArr = NULL; //Placeholder
  if(toCalcHash){ //Do we care about a 1d index?
    //cIndArr[j1num]x[3] = [j1]x[cIndArr, oddFlipArr, J=j1+j2+j3 ] (ptr arr)  
    cIndArr = (ullint**) malloc(j1Num*sizeof(ullint*));
    for(i=0; i<j1Num; i++){
      cIndArr[i] = (ullint*) malloc(3*sizeof(ullint));
    }
  
    //Calculate 1d indices
    for(i=0; i<j1Num; i++){
      double jArr[] = {j1Max-j1Num+i+1, jmArr_ptr[0], jmArr_ptr[1]}; //j1,j2,j3
      double mArr[] = {jmArr_ptr[2], jmArr_ptr[3], jmArr_ptr[4]}; //m1,m2,m3
 
      ullint *cRetArr_ptr = hashRegge(&jArr[0], &mArr[0]); //Calc index

      for(j=0; j<3; j++){
	cIndArr[i][j] = cRetArr_ptr[j];
      }
      free(cRetArr_ptr); //Free mem allocated in hashRegge
    
    }
  }
  /* Populate structure */
  w3jStruct_ptr->w3jArr_ptr = w3jArr_ptr;
  //Ensure jmVals are copied to struct
  double* jmValsCpy =  (double*) malloc(5*sizeof(double));
  for(i=0; i<5; i++){
    jmValsCpy[i] = jmArr_ptr[i];
  }
  w3jStruct_ptr->jmVals_ptr = &jmValsCpy[0];
  if(toCalcHash){ //Do we care about a 1d index?
    w3jStruct_ptr->cIndArr_ptr = &cIndArr[0];
  }  
  w3jStruct_ptr->j1Num = j1Num;
  w3jStruct_ptr->j1Min = j1Max-j1Num;
  
  return w3jStruct_ptr;
}

/* Free data structure */
void killW3jStruct(w3jStruct *w3jStruct_ptr, int hashCalced){
  int i; //Iterator(s)
  
  free(w3jStruct_ptr->w3jArr_ptr);
  free(w3jStruct_ptr->jmVals_ptr);
  if(hashCalced){ //Did we calculate a 1d index?
    for(i=0; i<w3jStruct_ptr->j1Num; i++){
      free(w3jStruct_ptr->cIndArr_ptr[i]);
    }  
    free(w3jStruct_ptr->cIndArr_ptr);
  }
  free(w3jStruct_ptr);
  
}

/* Print data within container */
void printW3jStruct(w3jStruct *w3jStruct_ptr, int hashCalced){
  ullint i; //Iterator(s)
  
  //Print j, m vals
  printf("w3jStruct(");
  for(i=0; i<5; i++){
    printf("%.1f ", w3jStruct_ptr->jmVals_ptr[i]);
  }
  printf("), j1Min=%.1f,j1Num=%llu\n\n",
	 w3jStruct_ptr->j1Min,
	 w3jStruct_ptr->j1Num);
  for(i=0; i<w3jStruct_ptr->j1Num; i++){
    printf("j1=%.1f, w3j=%.3f",
	   w3jStruct_ptr->j1Min+i,
	   w3jStruct_ptr->w3jArr_ptr[i]);
    if(hashCalced){ //Did we calculate a 1d index?
      printf(", cInd=(%llu,%llu,%llu) \n",
	     w3jStruct_ptr->cIndArr_ptr[i][0],
	     w3jStruct_ptr->cIndArr_ptr[i][1],
	     w3jStruct_ptr->cIndArr_ptr[i][2]);
    } else {
      printf("\n");
    }

  }
}

/* 
   Count the total number of valid symbols that exist
   (dl={0.5, 1}, L=maximum j1 value)
   (Almost no symmetries taken into account)
*/
ullint countW3jSymbols(double L, double dl){
  double l2, l3, m1, m2, m3; //Iterator(s)
  ullint sNum = 0;
  for(l2=0; l2<=L; l2+=dl){
    for(l3=0; l3<=l2; l3+=dl){
      for(m1=-L; m1<=L; m1++){
	for(m2=-l2; m2<=l2; m2++){
	  m3=-(m1+m2); //Condition on W3j
	  
	  //Ensure that (|m3|<=l3) /\ (|m3|%1==l3%1)
	  if((fabs(m3)<=l3) && 
	     (fmod(fabs(m3),1) == (fmod(l3,1)))){ 
	    //Total number of symbols in this family
	    sNum += (ullint)(l2+l3-fmax(fabs(l2-l3), fabs(m1))+1);
	  }
	}	
      }
    }
  }
  return sNum;
}

/*
  Append new data from w3jStruct to cplCw3jStruct,
  dynamically resize as required
*/
void appendw3jStructTocplCw3jStruct(w3jStruct *w3jStruct_ptr,
				    cplCw3jStruct *cplCw3jStruct_ptr){
  ullint i; //Iterator(s)

  //Current max size of the arrays
  ullint maxCplSz = cplCw3jStruct_ptr->sz;
  //Current filled size
  ullint curCplSz = cplCw3jStruct_ptr->curSz;
  //Size of data to append
  ullint appDatSz = w3jStruct_ptr->j1Num;
  
  //Resize if required
  if(maxCplSz-curCplSz<appDatSz){
    ullint newSz = maxCplSz;
    while(newSz-curCplSz<appDatSz){
      newSz *= 2; //Keep multiplying size until we can fit new data
    }

    /*
      Before increasing size, decrease to last entry,
      Sort this and remove repetitions.
      Then increase to new size.
    */
    resizecplCw3jStruct(cplCw3jStruct_ptr, curCplSz);
    mergeSortCplCw3j(cplCw3jStruct_ptr);
    cplCw3jStruct *rtmp = rmRepCplCw3j(cplCw3jStruct_ptr);
    
    //Free old data
    free(cplCw3jStruct_ptr->ullArr_ptr);
    free(cplCw3jStruct_ptr->dblArr_ptr);
    //Point to new data
    cplCw3jStruct_ptr->sz = rtmp->sz; //Update new size
    curCplSz = rtmp->sz; //Repetitions removed
    cplCw3jStruct_ptr->curSz = curCplSz;
    cplCw3jStruct_ptr->ullArr_ptr = rtmp->ullArr_ptr;
    cplCw3jStruct_ptr->dblArr_ptr = rtmp->dblArr_ptr;

    //Free temp data
    free(rtmp);
    
    //Resize (zero new locations)
    resizecplCw3jStruct(cplCw3jStruct_ptr, newSz);
    cplCw3jStruct_ptr->sz = newSz; //New maximal size
  }
  
  //Copy data
  for(i=curCplSz; i<curCplSz+w3jStruct_ptr->j1Num; i++){
    
    //Account for phase
    double ph; //=(-1)^(oddFlip*J)
    ullint oddFlip = w3jStruct_ptr->cIndArr_ptr[i-curCplSz][1];
    ullint J = w3jStruct_ptr->cIndArr_ptr[i-curCplSz][2];

    if(oddFlip*J==0){ //Instead of pow(-1, oddFlip*J);
      ph = 1.0;
    } else if(J%2>0.1){
      ph = -1.0;
    } else {
      ph = 1.0;
    }

    //Copy data and apply phase
    cplCw3jStruct_ptr->ullArr_ptr[i] = \
      w3jStruct_ptr->cIndArr_ptr[i-curCplSz][0];
    cplCw3jStruct_ptr->dblArr_ptr[i] = \
      ph*(w3jStruct_ptr->w3jArr_ptr[i-curCplSz]);
  }
  //Increase storage index
  cplCw3jStruct_ptr->curSz += w3jStruct_ptr->j1Num;
}


/*
  Calculate symbol set up to some maximal L (some may be past L)
  dl = 1 => SO(3)
  dl = 0.5 => SO(3) /\ SU(2)
*/
cplCw3jStruct* calcw3jSymbols(double L, double dl){
  double l2, l3, m1, m2, m3; //Iterator(s)
  ullint curSz = 0;
  cplCw3jStruct *cplCw3jStruct_ptr = constrSzCplCw3jStruct(W3JINISZ);
  cplCw3jStruct_ptr->curSz = curSz;

  for(l2=0; l2<=L; l2+=dl){
    for(l3=0; l3<=l2; l3+=dl){
      for(m1=-L; m1<=L; m1++){
	for(m2=-l2; m2<=l2; m2++){
	  m3=-(m1+m2); //Condition on W3j
	  
	  //Ensure that (|m3|<=l3) /\ (|m3|%1==l3%1)
	  if((fabs(m3)<=l3) && 
	     (fmod(fabs(m3),1) == (fmod(l3,1)))){ 
	    //Total number of symbols in this family
	    //ullint l1Num = (ullint)(l2+l3-fmax(fabs(l2-l3), fabs(m1))+1);
	    w3jStruct *w3jStruct_ptr;
	    double lmArr[] = {l2, l3, m1, m2, m3};      
	    w3jStruct_ptr = calculateW3jFamily(&lmArr[0], 1);
	    
	    //Update cplCw3jStruct
	    appendw3jStructTocplCw3jStruct(w3jStruct_ptr, cplCw3jStruct_ptr);
	    
	    //Free resources
	    killW3jStruct(w3jStruct_ptr,1);

	    /*
	      Should dump to a temporary output file after each outer loop
	      iterate for large L numbers
	    */

	  }
	}	
      }
    }

  }

  //Sort and remove trailing entries, then return
  resizecplCw3jStruct(cplCw3jStruct_ptr, cplCw3jStruct_ptr->curSz);
  mergeSortCplCw3j(cplCw3jStruct_ptr);

  //Remove repetitions
  cplCw3jStruct *rtmp = rmRepCplCw3j(cplCw3jStruct_ptr);    
  
  //Free old data
  free(cplCw3jStruct_ptr->ullArr_ptr);
  free(cplCw3jStruct_ptr->dblArr_ptr);

  cplCw3jStruct_ptr->sz = rtmp->sz; //Update new size
  ullint curCplSz = rtmp->sz; //Repetitions removed
  cplCw3jStruct_ptr->curSz = curCplSz;
  cplCw3jStruct_ptr->ullArr_ptr = rtmp->ullArr_ptr;
  cplCw3jStruct_ptr->dblArr_ptr = rtmp->dblArr_ptr;

  cplCw3jStruct_ptr->L = L;
  cplCw3jStruct_ptr->dl = dl;
  //Free temp data
  free(rtmp);

  return cplCw3jStruct_ptr;
}

/*
  Calculate a particular symbol with l1,l2,l3,m1,m2,m3 given,
  We assume that the choice is valid.
  Recursion is performed and a single value is extracted
*/
double calculateW3jVal(double *lArr_ptr, double *mArr_ptr){
  
  double jmArr_ptr[] = {lArr_ptr[1], lArr_ptr[2],
			mArr_ptr[0], mArr_ptr[1], mArr_ptr[2]};

  double j1Max = jmArr_ptr[0]+jmArr_ptr[1];
  ullint j1Num = (ullint)(j1Max-fmax(fabs(jmArr_ptr[0]-jmArr_ptr[1]),
				     fabs(jmArr_ptr[2]))+1);

  double* w3jArr = calloc(j1Num, sizeof(double));
  double *w3jArr_ptr = &w3jArr[0];
  calcCoeff(jmArr_ptr, w3jArr_ptr); //Perform recursion

  double j1Min = j1Max-j1Num+1;
  double toRet = w3jArr_ptr[(ullint) (lArr_ptr[0]-j1Min)];

  //Free resources
  free(w3jArr);

  return toRet;
}



int comparator(const void *a, const void *b){
  if(*(ullint*)a >  *(ullint*)b){ return 1; }
  if(*(ullint*)a == *(ullint*)b){ return 0; }
  if(*(ullint*)a <  *(ullint*)b){
    return -1;
  } else {
    return 1; //Value not present.
  }
} //Three statements as return type must be int
/*
  Find array index of a particular cInd within a 1d array.
  Uses Binary search algorithm (stdlib implementation) with O(log(n))
  This function assumes that the index exists.
*/
ullint findArrInd(ullint cInd, cplCw3jStruct* cplCw3jStruct_ptr){
  ullint *cIndPos_ptr = (ullint*) bsearch(&cInd, //key to find
					  cplCw3jStruct_ptr->ullArr_ptr, //srch
					  cplCw3jStruct_ptr->sz, //size of arr 
					  sizeof(ullint), //size of type
					  comparator); //comparing function

  return cIndPos_ptr-cplCw3jStruct_ptr->ullArr_ptr; //Return array index
}


/*
  Given a set of value {l1,l2,l3,m1,m2,m3} return the correspond w3j symbol
  from precalculation.
*/
double getW3jPrecalc(double *lArr_ptr, double *mArr_ptr,
		     cplCw3jStruct *cplCw3jStruct_ptr){

  ullint *hash_ptr = hashRegge(lArr_ptr, mArr_ptr); //Regge index
  ullint arrInd = findArrInd(hash_ptr[0], cplCw3jStruct_ptr); //Array position

  //Account for phase
  double ph; //=(-1)^(oddFlip*J)
  ullint oddFlip = hash_ptr[1]; //w3jStruct_ptr->cIndArr_ptr[i-curCplSz][1];
  ullint J = hash_ptr[2]; //w3jStruct_ptr->cIndArr_ptr[i-curCplSz][2];

  if(oddFlip*J==0){ //Instead of pow(-1, oddFlip*J);
    ph = 1.0;
  } else if(J%2>0.1){
    ph = -1.0;
  } else {
    ph = 1.0;
  }

  //Free resources
  free(hash_ptr);

  return ph*(cplCw3jStruct_ptr->dblArr_ptr[arrInd]); //Value we care about
}
