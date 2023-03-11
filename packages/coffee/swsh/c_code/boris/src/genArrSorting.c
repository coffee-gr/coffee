/*
  Library containing various sorting algorithms
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/* Macros */
#include "macros.h"

/* Generic functions */
#include "genArrSorting.h" //Array sorting

/* Structure for sorting */
typedef struct _sortArr{
  ullint sz;
  double *Arr_ptr;
} sortArr;

/* Local macros */
#define mergeSortSwitchSz = 32; //Switch to insertion sort below this len

/*
  To Do:
   Switch merge sorted to insertion sort for sufficiently small structs
*/

/* Function prototypes */
sortArr* constrArrS(sortArr *toDup_ptr);
sortArr* duplicateArrS(sortArr *toDup_ptr);
void printArrS(sortArr *toPrint_ptr);
void insertionSortS(sortArr *toSort_ptr);

void mergeS(sortArr *toSort_ptr, sortArr *tmp_ptr,
	    ullint left, ullint right, ullint rightEnd);


void mergeSortS(sortArr *toSort_ptr, sortArr *tmp_ptr,
		ullint left, ullint right);
//Construct from input data
cplCw3jStruct* constrCplCw3jStruct(ullint sz,
				   ullint *ullArr_ptr,
				   double *dblArr_ptr);
//Construct zero filled struct of specified size
cplCw3jStruct* constrSzCplCw3jStruct(ullint sz);
//Allocate memory for tmp struct
cplCw3jStruct* qkDupCplCw3jStruct(cplCw3jStruct *toDup_ptr);
//Realloc internal arrays
void resizecplCw3jStruct(cplCw3jStruct *cplCw3jStruct_ptr, ullint newSz);
//Free resources
void killCplCw3jStruct(cplCw3jStruct *toKill_ptr);
void printCplCw3jStruct(cplCw3jStruct *toPrint_ptr);
//Sorting prototypes for coupled struct
void insertionSortCplCw3j(cplCw3jStruct *cplCw3jStruct_ptr);
void mergeSortCplCw3j(cplCw3jStruct *cplCw3jStruct_ptr);
cplCw3jStruct* rmRepCplCw3j(cplCw3jStruct *cplCw3jStruct_ptr);

/*
int main(){

  // Generate random array
  
//  for(i=0; i<toSortRand_sz; i++){
//    toSortRand[i] = toSortRand_sz-i; 
//    //toSortRand[i] = rand()%toSortRand_sz+		\
//	//  (double) toSortRand_sz/(rand()%toSortRand_sz+1);
//  }




  //sortArr* tsRand = (sortArr*) malloc(sizeof(sortArr));
  sortArr *tsRand_ptr = (sortArr*) malloc(sizeof(sortArr));
  ullint tsRand_sz = 10;
  double* tsRandArr = (double*) malloc(tsRand_sz*sizeof(double));
  //double tsRandArr[] = {1.3, 4.5, 7.4, 0.1, 1.1, 3.4, 6.2};
  double *tsRandArr_ptr = &tsRandArr[0];

  //  (*tsRand).sz = tsRand_sz;
  
  tsRand_ptr->sz = tsRand_sz;
  tsRand_ptr->Arr_ptr = tsRandArr_ptr;

  ullint j; //Iterator
  srand(time(NULL)); //Seed random number generator
  for(j=0; j<tsRand_ptr->sz; j++){
    tsRand_ptr->Arr_ptr[j] = tsRand_ptr->sz-j;
  }
  
  //  printArrS(tsRand_ptr);

  //insertionSortS(tsRand_ptr);
  //printArrS(tsRand_ptr);

  // Duplicate current struct to sort 
  sortArr *tmp_ptr = constrArrS(tsRand_ptr);
  //  printArrS(tmp_ptr);

  //mergeS(tsRand_ptr, tmp_ptr, 0, 3, (tsRand_ptr->sz)-1);
  //printArrS(tsRand_ptr);
  mergeSortS(tsRand_ptr, tmp_ptr, 0, (tsRand_ptr->sz)-1);
   printf("aa=%llu\n", (tsRand_ptr->sz)-1);

   //insertionSortSran(tsRand_ptr, 4, 2);

  printArrS(tsRand_ptr);
  return 0;

}
*/


/* Construct the dbl sorting structure */
sortArr* constrArrS(sortArr *toDup_ptr){
  sortArr *ret_ptr = (sortArr*) malloc(sizeof(sortArr));
  double* dblArr = (double*) malloc((toDup_ptr->sz)*sizeof(double));
  double *dblArr_ptr = &dblArr[0];
  
  ret_ptr->sz = toDup_ptr->sz;
  ret_ptr->Arr_ptr = dblArr_ptr;

  return ret_ptr;
}

/* Duplicate the dbl sorting structure */
sortArr* duplicateArrS(sortArr *toDup_ptr){
  sortArr *ret_ptr = (sortArr*) malloc(sizeof(sortArr));
  double* dblArr = (double*) malloc((toDup_ptr->sz)*sizeof(double));
  double *dblArr_ptr = &dblArr[0];
  
  ret_ptr->sz = toDup_ptr->sz;
  ret_ptr->Arr_ptr = dblArr_ptr;
  //Populate
  ullint i;
  for(i=0; i<toDup_ptr->sz; i++){
    ret_ptr->Arr_ptr[i] = toDup_ptr->Arr_ptr[i];
  }
  return ret_ptr;
}


void printArrS(sortArr *toPrint_ptr){
  ullint i=0;
  printf("sortArr=(");
  for(i=0; i<toPrint_ptr->sz; i++){
    printf("%.1f ", toPrint_ptr->Arr_ptr[i]);
  }
  printf(")\n");
}


/* Insertion sorting O(n^2) */
void insertionSortS(sortArr *toSort_ptr){
  ullint i, j; //Iterator(s)
  double tmp;
  
  for(i=1; i<toSort_ptr->sz; i++){
    tmp = toSort_ptr->Arr_ptr[i];
    for(j=i-1; j>=0; j--){
      if(toSort_ptr->Arr_ptr[j]<=tmp)
	break;
      toSort_ptr->Arr_ptr[j+1]=toSort_ptr->Arr_ptr[j];
    }
    toSort_ptr->Arr_ptr[j+1]=tmp;
  }
}


/* Merge sorted sublists contained within struct */
void mergeS(sortArr *toSort_ptr, sortArr *tmp_ptr,
	    ullint left, ullint right, ullint rightEnd){
  ullint i, num, temp, leftEnd=right-1;
  temp = left;
  num = rightEnd-left+1;

  while((left<=leftEnd)&&(right<=rightEnd)){
    if(toSort_ptr->Arr_ptr[left]<=toSort_ptr->Arr_ptr[right]){
      tmp_ptr->Arr_ptr[temp++] = toSort_ptr->Arr_ptr[left++];
    } else {
      tmp_ptr->Arr_ptr[temp++] = toSort_ptr->Arr_ptr[right++];
    }
  }

  while(left<=leftEnd){
    tmp_ptr->Arr_ptr[temp++] = toSort_ptr->Arr_ptr[left++];
  }
  while(right<=rightEnd){
    tmp_ptr->Arr_ptr[temp++] = toSort_ptr->Arr_ptr[right++];
  }

  for(i=1; i<=num; i++, rightEnd--){
    toSort_ptr->Arr_ptr[rightEnd] = tmp_ptr->Arr_ptr[rightEnd];
  }

}

/* Merge sort the data structure */
void mergeSortS(sortArr *toSort_ptr, sortArr *tmp_ptr,
		ullint left, ullint right){
  ullint centre;
  if(left<right){
    centre=(left+right)/2;
    //printf("%llu, %llu, %llu, %llu\n", left, right, centre, right-left);
    mergeSortS(toSort_ptr, tmp_ptr, left, centre);
    mergeSortS(toSort_ptr, tmp_ptr, centre+1, right);
    mergeS(toSort_ptr, tmp_ptr, left, centre+1, right);
  }
}


/*
  Create a cplCw3jStruct from input data
*/
cplCw3jStruct* constrCplCw3jStruct(ullint sz, 
				   ullint *ullArr_ptr,
				   double *dblArr_ptr){
  cplCw3jStruct *ret_ptr = (cplCw3jStruct*) malloc(sizeof(cplCw3jStruct));
  
  //Populate returned struct
  ret_ptr->sz = sz;
  ret_ptr->dblArr_ptr = dblArr_ptr;
  ret_ptr->ullArr_ptr = ullArr_ptr;
  
  return ret_ptr;
}

/*
  Create a cplCw3jStruct of specified sz (zero filled)
*/
cplCw3jStruct* constrSzCplCw3jStruct(ullint sz){
  ullint i; //Iterator(s)
  cplCw3jStruct *ret_ptr = (cplCw3jStruct*) malloc(sizeof(cplCw3jStruct));

  double* dblArr = (double*) malloc(sz*sizeof(double));
  double *dblArr_ptr = &dblArr[0];
  ullint* ullArr = (ullint*) malloc(sz*sizeof(ullint));
  ullint *ullArr_ptr = &ullArr[0];
  
  //Zero fill
  for(i=0; i<sz; i++){
    dblArr_ptr[i] = (double) 0;
    ullArr_ptr[i] = (ullint) 0;
  }
  //Populate returned struct
  ret_ptr->sz = sz;
  ret_ptr->dblArr_ptr = dblArr_ptr;
  ret_ptr->ullArr_ptr = ullArr_ptr;
  
  return ret_ptr;
}


/*
  Quickly duplicate cplCw3jStruct (malloc but do not dup. elements)
*/
cplCw3jStruct* qkDupCplCw3jStruct(cplCw3jStruct *toDup_ptr){
  cplCw3jStruct *ret_ptr = (cplCw3jStruct*) malloc(sizeof(cplCw3jStruct));
  double* dblArr = (double*) malloc((toDup_ptr->sz)*sizeof(double));
  double *dblArr_ptr = &dblArr[0];
  ullint* ullArr = (ullint*) malloc((toDup_ptr->sz)*sizeof(ullint));
  ullint *ullArr_ptr = &ullArr[0];

  //Populate returned struct
  ret_ptr->sz = toDup_ptr->sz;
  ret_ptr->dblArr_ptr = dblArr_ptr;
  ret_ptr->ullArr_ptr = ullArr_ptr;
  
  return ret_ptr;
}

/*
  Resize arrays within cplCw3jStruct, zero fill
*/
void resizecplCw3jStruct(cplCw3jStruct *cplCw3jStruct_ptr, ullint newSz){
  ullint i; //Iterator(s)
  
  ullint* ullTmp = (ullint*) realloc(cplCw3jStruct_ptr->ullArr_ptr,
				     newSz*sizeof(ullint));
  
  double* dblTmp = (double*) realloc(cplCw3jStruct_ptr->dblArr_ptr,
				     newSz*sizeof(double));
  cplCw3jStruct_ptr->ullArr_ptr = &ullTmp[0];
  cplCw3jStruct_ptr->dblArr_ptr = &dblTmp[0];
  
  //Zero fill any additional mem locations
  for(i=cplCw3jStruct_ptr->sz; i<newSz; i++){
    cplCw3jStruct_ptr->ullArr_ptr[i] = (ullint) 0;
    cplCw3jStruct_ptr->dblArr_ptr[i] = (double) 0;
  }
  
  cplCw3jStruct_ptr->sz = newSz;
  
}


/* Free resources */
void killCplCw3jStruct(cplCw3jStruct *toKill_ptr){
  free(toKill_ptr->ullArr_ptr);
  free(toKill_ptr->dblArr_ptr);
  free(toKill_ptr);
}

void printCplCw3jStruct(cplCw3jStruct *toPrint_ptr){
  ullint i; //Iterator(s)
  
  printf("cplCw3jStruct->sz=%llu\n", toPrint_ptr->sz);
  for(i=0; i<toPrint_ptr->sz; i++){
    printf("(cInd[%llu], w3j[%llu])=(%llu, %.3f)\n",
	   i, i,
	   toPrint_ptr->ullArr_ptr[i],
	   toPrint_ptr->dblArr_ptr[i]);
  }
}

/* 
   Insertion sorting O(n^2) of coupled cInd / w3jArr.
   Sort by cInd, by keep corresponding w3j val aligned.
*/
void insertionSortCplCw3j(cplCw3jStruct *cplCw3jStruct_ptr){
  ullint i, j; //Iterator(s)
  double dbltmp;
  ullint ulltmp;

  for(i=1; i<cplCw3jStruct_ptr->sz; i++){
    //Copy current element
    dbltmp = cplCw3jStruct_ptr->dblArr_ptr[i];
    ulltmp = cplCw3jStruct_ptr->ullArr_ptr[i];

    //for(j=i-1; j>=0; j--){
    //ullint => dec iter var at beginning of loop, not end
    for(j=i-1+1; j-->0; ){
      //Sort by ull index
      if(cplCw3jStruct_ptr->ullArr_ptr[j]<=ulltmp)
	break;
	
      //Copy if required
      cplCw3jStruct_ptr->ullArr_ptr[j+1] = cplCw3jStruct_ptr->ullArr_ptr[j];
      cplCw3jStruct_ptr->dblArr_ptr[j+1] = cplCw3jStruct_ptr->dblArr_ptr[j];
    }
    cplCw3jStruct_ptr->ullArr_ptr[j+1] = ulltmp;
    cplCw3jStruct_ptr->dblArr_ptr[j+1] = dbltmp;
    
  }

}


/* 
   Merge sorted sublists contained within struct.
   ullArr is sorted, position of dblArr is aligned
*/
void mergeCplCw3jS(cplCw3jStruct *cplCw3jStruct_ptr, 
		   cplCw3jStruct *cplCw3jStructTmp_ptr,
		   ullint left, ullint right, ullint rightEnd){
  ullint i, num, temp, leftEnd=right-1;
  temp = left;
  num = rightEnd-left+1;

  while((left<=leftEnd)&&(right<=rightEnd)){
    //toSort--cplCw3jStruct
    //tmp->cplCw3jStructTmp
    if(cplCw3jStruct_ptr->ullArr_ptr[left] <=	\
       cplCw3jStruct_ptr->ullArr_ptr[right]){
      cplCw3jStructTmp_ptr->ullArr_ptr[temp] =	\
	cplCw3jStruct_ptr->ullArr_ptr[left];
      cplCw3jStructTmp_ptr->dblArr_ptr[temp++] =	\
	cplCw3jStruct_ptr->dblArr_ptr[left++];
    } else {
      cplCw3jStructTmp_ptr->ullArr_ptr[temp] =	\
	cplCw3jStruct_ptr->ullArr_ptr[right];
      cplCw3jStructTmp_ptr->dblArr_ptr[temp++] =	\
      	cplCw3jStruct_ptr->dblArr_ptr[right++];
    }
  }

  while(left<=leftEnd){
    cplCw3jStructTmp_ptr->ullArr_ptr[temp] =	\
      cplCw3jStruct_ptr->ullArr_ptr[left];
    cplCw3jStructTmp_ptr->dblArr_ptr[temp++] =	\
      cplCw3jStruct_ptr->dblArr_ptr[left++];
  }
  while(right<=rightEnd){
    cplCw3jStructTmp_ptr->ullArr_ptr[temp] =	\
      cplCw3jStruct_ptr->ullArr_ptr[right];
    cplCw3jStructTmp_ptr->dblArr_ptr[temp++] =		\
      cplCw3jStruct_ptr->dblArr_ptr[right++];
  }
  
  for(i=1; i<=num; i++, rightEnd--){
    cplCw3jStruct_ptr->ullArr_ptr[rightEnd] =		\
      cplCw3jStructTmp_ptr->ullArr_ptr[rightEnd];
    cplCw3jStruct_ptr->dblArr_ptr[rightEnd] =	\
      cplCw3jStructTmp_ptr->dblArr_ptr[rightEnd];
  }

}


/*
 Merge sort the data structure.
 ullArr is sorted, dblArr is coupled
*/
void mergeSortCplCw3jS(cplCw3jStruct *cplCw3jStruct_ptr, 
		       cplCw3jStruct *cplCw3jStructTmp_ptr,
		       ullint left, ullint right){
  ullint centre;
  if(left<right){
    centre=(left+right)/2;
    //printf("%llu, %llu, %llu, %llu\n", left, right, centre, right-left);
    mergeSortCplCw3jS(cplCw3jStruct_ptr, 
		      cplCw3jStructTmp_ptr,
		      left, centre);
    mergeSortCplCw3jS(cplCw3jStruct_ptr, 
		      cplCw3jStructTmp_ptr,
		      centre+1, right);
    mergeCplCw3jS(cplCw3jStruct_ptr, 
		  cplCw3jStructTmp_ptr,
		  left, centre+1, right);
  }
}


/*
  Merge sort the cplCw3jStruct, 
  ullArr_ptr is sorted with dblArr_ptr kept aligned
*/
void mergeSortCplCw3j(cplCw3jStruct *cplCw3jStruct_ptr){
  //Construct a temporary array
  cplCw3jStruct *tmpCw3jStruct_ptr = qkDupCplCw3jStruct(cplCw3jStruct_ptr);
  
  mergeSortCplCw3jS(cplCw3jStruct_ptr, 
		    tmpCw3jStruct_ptr,
		    0, (cplCw3jStruct_ptr->sz)-1);

  killCplCw3jStruct(tmpCw3jStruct_ptr); //Free resources
}

/*
  Remove repeated ullArr entries (and resp. aligned dblArr)
  Input must be sorted
*/
cplCw3jStruct* rmRepCplCw3j(cplCw3jStruct *cplCw3jStruct_ptr){
  ullint i, j=0; //Iterator(s)

  //Construct a temporary array
  cplCw3jStruct *tmpCw3jStruct_ptr = qkDupCplCw3jStruct(cplCw3jStruct_ptr);
  if(cplCw3jStruct_ptr->sz>0){
    tmpCw3jStruct_ptr->ullArr_ptr[0] = cplCw3jStruct_ptr->ullArr_ptr[0];
    tmpCw3jStruct_ptr->dblArr_ptr[0] = cplCw3jStruct_ptr->dblArr_ptr[0];
    j++;
  }
  for(i=1; i<cplCw3jStruct_ptr->sz; i++){
    if(cplCw3jStruct_ptr->ullArr_ptr[i-1] != \
       cplCw3jStruct_ptr->ullArr_ptr[i]){
      
      tmpCw3jStruct_ptr->ullArr_ptr[j] = cplCw3jStruct_ptr->ullArr_ptr[i];
      tmpCw3jStruct_ptr->dblArr_ptr[j] = cplCw3jStruct_ptr->dblArr_ptr[i];
      j++;

    }
  }

  tmpCw3jStruct_ptr->sz = j;
  return tmpCw3jStruct_ptr;
}

