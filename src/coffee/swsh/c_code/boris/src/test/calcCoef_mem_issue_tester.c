#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "macros.h"
#include "mathWigner3jRecursion.h"


int main(){
  double *jms = calloc(5, sizeof(double));
  jms[0] = 2;
  jms[1] = 2;
  jms[2] = 0;
  jms[3] = 0;
  jms[4] = 0;
  double *jmArr_ptr = &jms[0];
  
  double j1Max = jmArr_ptr[0]+jmArr_ptr[1];
  int j1Num = (ullint)(j1Max-fmax(fabs(jmArr_ptr[0]-jmArr_ptr[1]),
				     fabs(jmArr_ptr[2]))+1);
  
  double* w3jArr = calloc(j1Num, sizeof(double));
  double *w3jArr_ptr = &w3jArr[0];
  
  //calcCoeff(jmArr_ptr, w3jArr_ptr);
  
  free(jms);
  free(w3jArr);

  return 0;
}
