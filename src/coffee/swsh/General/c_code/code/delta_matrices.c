/////////////////////////////////////////////////////
// Dr. Leon Escobar
// Department of Mathematics,
// Universidad del Valle,
// Colombia
// 2019
//////////////////////////////////////////////////////

#include <delta_matrices.h>
#include <omp.h>


//-----------------------------------------------------
void Compute_Delta_Matrices( double *** restrict D ,  int lmax ){ 
  // the counters
  int m1,m2,l;  
  //initializiing Delta matrices
  D[0][0][0]=1.0;

  for( l=1; l <= lmax ; l++ ){ //computing planes lmax-planes of Delta matrices      
    // moving in dwn in l
    for( m2=0; m2 <=l ; m2++ ){ 
      
      if(m2==0){ // moving down for m2=0  
      D[l][l][0] = - sqrt( (2.*l-1)/(2.*l) )* D[l-1][l-1][0];
      }
      else{ // moving in diagonal in l and m2
      D[l][l][m2] = sqrt( ((l/2.)*(2.*l-1)) / ((l+m2)*(l+m2-1.)) ) * D[l-1][l-1][m2-1]; 
      }        
      // moving in m1 (backwards sence!)
      for( m1 = l - 1 ; m1>= 0  ; m1-- ){      
          if(m1==l-1){
            D[l][l-1][m2] = ( (2.*m2)/sqrt((l-m1)*(l+m1+1.)) )* D[l][l][m2];
          }
          else{ 
            D[l][m1][m2] = ( (2.*m2)/sqrt((l-m1)*(l+m1+1.)) )* D[l][m1+1][m2] - sqrt( ( (l-m1-1.)*(l+m1+2.) )/( (l-m1)*(l+m1+1.) ) ) * D[l][m1+2][m2];       
        }             
      }//loop m1 
    }//loop m2
  }//loop l 
}//function


//-----------------------------------------------------
//prcompute delta matrices
double ***Precompute_Delta_Matrices( int lmax ){
  //Creating the 3 dimensional array
  double ***D = malloc( (lmax+1) * sizeof(double **) );
  Initialize_Delta_Matrices( D, lmax ); 

  //Computing Deltas
  Compute_Delta_Matrices(  D,  lmax ); // otra diferencia con respecto a la general
  //Show_Delta_Matrices( Deltas,  abs(s) + 1  ,  lmax );  
   return D;
}

///////////////////////////////////VERSION 2.0///////////////////////////////

//----------------------Functions for fillling the DD amtrices---------------
double Dlmn(double *** restrict D, int l, int m, int n){ // mirror properties!!!

  if (0<=m && 0<=n){
    return D[l][m][n];
  }
  else if(0<=m && n<0){
    return pow(-1,l+m) * D[l][m][abs(n)];
  }

  else if (m<0 && 0<=n){    
    return pow(-1,l+n) * D[l][abs(m)][n];
  }
  else{
    return pow(-1,2*l-m+n) * D[l][abs(m)][abs(n)];
  }

}


//------------------ factors for speed up!!!-------------------
fftw_complex FACTOR1(int s,int l,int m,int n){
return pow(-1.0 , s + l + m) * cpow(I,s+n) * sqrt( ( 2.*l + 1 )/(4. * M_PI) );
}

fftw_complex FACTOR2(int s, int l, int m, int n){
return pow(-1.0, 2*l - 2*m + s+n);
}


//------------------chosing the correct pyramidal index----------
int pindex(int l,int i){
 
 if(0<=i){
  return i;
 }
 else{
 return 2*l+1+i ;
 }

}


//-----------------------------------------------------
fftw_complex *****Precompute_Delta_Delta_Matrices( int lmax , int smax ){ 

//(speed up: this is the big ingredient for speeding up the code!!!) 

double ***D = malloc( (lmax+1) * sizeof(double**) );
Initialize_Delta_Matrices( D, lmax ); 
Compute_Delta_Matrices(  D,  lmax );


//allocating memory for the matrices 
// DD forward, DD backward
fftw_complex *****DD = calloc( 2 , sizeof(fftw_complex****) );


//////////////////////////// parallel seccion ///////////////////////

//Allocating memory for the DD matrices in pyramidal form!!!!
int p,i,j,k;
for( p=0; p <= 2 ; p++ ){ 
  // DD[0] for the forward transform
  // DD[1] for the backward transform
   DD[p] = calloc( (lmax+1) , sizeof(fftw_complex***)); 
   for( i=0; i < lmax+1 ; i++ ){
     DD[p][i] = calloc( 2* smax + 1 , sizeof(fftw_complex**));  
     for( j=0; j < 2*smax+1  ; j++ ){
        DD[p][i][j] = calloc( 2*i+1 , sizeof(fftw_complex*));
        for( k=0; k < 2*i+1 ; k++ ){   
          DD[p][i][j][k] = calloc( 2*i+1 , sizeof(fftw_complex)); 
          //printf("p=%i,i=%i,j=%i, k=%i, \n ",p,i,j,k );
        }
      } 
   }
} 
//////////////////////////////////////////////////////////////////////

// the counters
int l,s,m,n;


double Tolerance = pow(10,-14);
//#pragma omp parallel 
//    {
//    #pragma omp for

//indices are DD[l][s][m][n] (for speed up!)
int Inf(int a,int b){
 if(a<=b){
 return a;
 } 
 else{
 return b;
 }
}

int m_counter=0;
int sm_counter=0;

//int ***Index = calloc( lmax, sizeof(int***) );


for( l=0; l<=lmax ; l++  ){  

 //Index[l] = calloc( Inf(smax,l) , sizeof(int**) );

 for( s=0; s <= Inf(smax,l) ; s++){ //only computes up to the necesary spin weight...this is for saving memory!!!
 
  //Index[l][s] = calloc( l, sizeof(int*) );

  for( n=0; n<=l ; n++ ){ 


    for( m=0; m <=l ;  m++ ){ 

            sm_counter +=1;

           // We split the weight possibilities for speed up the code!!!  
            if ( Tolerance <  fabs( Dlmn(D, l,m,n) * Dlmn(D,l,m,s) ) ){
            // fabs ....the absolute value of an double number   
            
            //m_counter += 1; 

            //if ( Tolerance <  Dlmn(D, l,m,n) * Dlmn(D,l,m,s) ||  Dlmn(D, l,m,n) * Dlmn(D,l,m,s) < -1.0* Tolerance  ){

            //(s+, m+, n+)       
            DD[0][l][pindex(smax,s)][pindex(l,m)][pindex(l,n)] =  FACTOR1(s,l,m,n) * Dlmn(D, l,m,n) * Dlmn(D,l,m,s) ; 
            DD[1][l][pindex(smax,s)][pindex(l,m)][pindex(l,n)] =  FACTOR2(s,l,m,n) * DD[0][l][pindex(smax,s)][pindex(l,m)][pindex(l,n)] ; 
           
            //(s+, m+, n-) 
            if(0<n){
            DD[0][l][pindex(smax,s)][pindex(l,m)][pindex(l,-n)] =  FACTOR1(s,l,m,-n) * Dlmn(D, l,m,-n) * Dlmn(D,l,m,s) ; 
            DD[1][l][pindex(smax,s)][pindex(l,m)][pindex(l,-n)] =  FACTOR2(s,l,m,-n) * DD[0][l][pindex(smax,s)][pindex(l,m)][pindex(l,-n)] ; 
            }

            //(s+, m-, n+) 
            if(0<m){
            DD[0][l][pindex(smax,s)][pindex(l,-m)][pindex(l,n)] =  FACTOR1(s,l,-m,n) * Dlmn(D, l,-m,n) * Dlmn(D,l,-m,s) ; 
            DD[1][l][pindex(smax,s)][pindex(l,-m)][pindex(l,n)] =  FACTOR2(s,l,-m,n) * DD[0][l][pindex(smax,s)][pindex(l,-m)][pindex(l,n)] ; 
            }

            //(s+, m-, n-)
            if(0<m && 0<n){
            DD[0][l][pindex(smax,s)][pindex(l,-m)][pindex(l,-n)] =  FACTOR1(s,l,-m,-n) * Dlmn(D, l,-m,-n) * Dlmn(D,l,-m,s) ; 
            DD[1][l][pindex(smax,s)][pindex(l,-m)][pindex(l,-n)] =  FACTOR2(s,l,-m,-n) * DD[0][l][pindex(smax,s)][pindex(l,-m)][pindex(l,-n)] ; 
            }

            //(s-, m+, n+)
            if(0<s){
            DD[0][l][pindex(smax,-s)][pindex(l,m)][pindex(l,n)] =  FACTOR1(-s,l,m,n) * Dlmn(D, l,m,n) * Dlmn(D,l,m,-s) ; 
            DD[1][l][pindex(smax,-s)][pindex(l,m)][pindex(l,n)] =  FACTOR2(-s,l,m,n) * DD[0][l][pindex(smax,-s)][pindex(l,m)][pindex(l,n)] ; 
            }

            //(s-, m+, n-)
            if(0<s && 0<n){
            DD[0][l][pindex(smax,-s)][pindex(l,m)][pindex(l,-n)] =  FACTOR1(-s,l,m,-n) * Dlmn(D, l,m,-n) * Dlmn(D,l,m,-s) ; 
            DD[1][l][pindex(smax,-s)][pindex(l,m)][pindex(l,-n)] =  FACTOR2(-s,l,m,-n) * DD[0][l][pindex(smax,-s)][pindex(l,m)][pindex(l,-n)] ; 
            }

            //(s-, m-, n+)
            if(0<s && 0<m){
            DD[0][l][pindex(smax,-s)][pindex(l,-m)][pindex(l,n)] =  FACTOR1(-s,l,-m,n) * Dlmn(D, l,-m,n) * Dlmn(D,l,-m,-s) ; 
            DD[1][l][pindex(smax,-s)][pindex(l,-m)][pindex(l,n)] =  FACTOR2(-s,l,-m,n) * DD[0][l][pindex(smax,-s)][pindex(l,-m)][pindex(l,n)] ; 
            }

            //(s-, m-, n-)
            if(0<s && m && 0<n){
            DD[0][l][pindex(smax,-s)][pindex(l,-m)][pindex(l,-n)] =  FACTOR1(-s,l,-m,-n) * Dlmn(D, l,-m,-n) * Dlmn(D,l,-m,-s) ; 
            DD[1][l][pindex(smax,-s)][pindex(l,-m)][pindex(l,-n)] =  FACTOR2(-s,l,-m,-n) * DD[0][l][pindex(smax,-s)][pindex(l,-m)][pindex(l,-n)] ; 
            }

           //}
           }

         else{   
        
          m_counter +=1;
         // printf("s=%i, l=%i,,m=%i, n=%i, abs=%f \n",s,l,m,n, fabs( Dlmn(D, l,m,n) * Dlmn(D,l,m,s) ) );  
         }

           // printf("--------------------------\n" );  
        
         }//loop m       

         
         //Index[l][s][n] = calloc( m_counter , sizeof(int*) );
 

      }//loop n

    }//loop l
  }//loop s 
//}


// printf("sm_counter  %i, counter %i\n", 8*sm_counter,8* m_counter);
// printf("reason %e \n",  sm_counter/m_counter );

printf("TEST: %.16f + %.16fi", creal(D[0][0][0][0][0]),imag(D[0][0][0][0][0]));

return DD;

}//function
