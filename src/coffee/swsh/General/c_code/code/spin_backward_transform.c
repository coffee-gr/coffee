
/////////////////////////////////////////////////////////////////////
// Dr. Leon Escobar
// Departmento de Matematicas,
// Universidad del Valle,
// Colombia
// 2019
/////////////////////////////////////////////////////////////////////
//
// Note: 1. Here for all the 2D-arrays we will ue row-major ordering 
//       2. My m is equivalent to mp from HW-code 
//
/////////////////////////////////////////////////////////////////////


#include <spin_backward_transform.h>
#include <time.h>
#include <stdlib.h>


//-----------------------------------------------------
void Compute_Gmn( fftw_complex * restrict aln, fftw_complex * restrict Gmn, fftw_complex ***** restrict DD,  int s, int lmax, int smax, int Sampling_on_torus_theta, int Sampling_on_torus_phi ){ 
//Crating and computing Gmo ( m can takes negatives values )
/////// for speeding up the code we split the fillign of Gmn by parts /////

int l,n,m;
int s_index = pindex(smax,s);
int l_centro;//(speed up)
int abss = abs(s);
double Tolerance = pow(10,-13);


//case m=n=0
for(l=abss; l<=lmax; l+=2){ // DD is zero for l=odd!!!

    if ( Tolerance <  cabs( DD[1][l][ s_index ][0][0]  ) ){ 

   // printf("l=%i,s=%i,m=%i, n=%i, \n ",l,s, 0,0 );
   Gmn[0] += DD[1][ l ][ s_index ][0][0] * aln[ l*(l+1) ];  
 }

      // printf("l,m,n, index =%i,%i,%i,%i\n",loo,0,0,0); 
    // printf("!DD[1][%i][%i][%i][%i]=%f \n ",index_s,loo,0, 0, creal( DD[1][ index_s ][loo][0][0] ) ); 
 }

//printf("Gmn[%i] = %f + %f i \n",0, creal( Gmn[ 0] ) , cimag( Gmn[ 0] ) ); 

//------------------for the index-------------------
int Sup(int a,int b){
 if(a<=b){
 return b;
 } 
 else{
 return a;
 }
}
//-------------------------------------------


//printf("==========================case m=0=========================\n");

// case m=0
//initializing n and l!
n=0;
while( n <= lmax ){ 
   n+=1;// raising the n!     
   l= Sup( abss ,n); 
   //l=max(abss,n);
   while (l<=lmax){
        l_centro= l*(l+1); 
       
        //printf("!!!l=%i,m=%i, n=%i, \n ",l, 0,n );
        // printf("DD[1][%i][%i][%i][%i]=%f+%fi \n ", s_index,l,0, n, creal( DD[1][l][ s_index ][0][n] ) , cimag( DD[1][l][  s_index ][0][n] ) ); 
    
         if ( Tolerance <  cabs( DD[1][l][ s_index ][0][n]  ) ){

        // (n+)
         Gmn[ n ] += DD[1][ l ][ s_index ][0][ n ] * aln[ l_centro + n ];  
        //Gmn[ n ] += aln[ l_centro + n ];  

        // (n-)
         Gmn[ Sampling_on_torus_phi - n ] += DD[1][l][ s_index ][0][ 2*l+1 - n ] * aln[ l_centro - n ];         
        //Gmn[ Sampling_on_torus_phi - n ] += aln[ l_centro - n ];         
 
        }
       // printf("l,m,n, index =%i,%i,%i,%i\n",lo,0,no,   no ); 
       // printf("DD[1][%i][%i][%i][%i]=%f \n ",index_s,lo,0, no, cimag( DD[1][ index_s ][lo][0][no] ) ); 
       // printf("DD[1][%i][%i][%i][%i]=%f \n ",index_s,lo,0, number_of_n-no, cimag( DD[1][ index_s ][lo][0][number_of_n-no] ) ); 
       // printf("l=%i\n",l);
      
         l+=2; // DD is zero for l+n=odd!!!  
   }  
   
  //printf("Gmn[%i] = %f + %f i \n",n , creal( Gmn[ n ] ) , cimag( Gmn[   n ] ) );   
  //printf("Gmn[%i] = %f + %f i \n",(Sampling_on_torus_phi - n)  , creal( Gmn[ (Sampling_on_torus_phi - n)  ]) , cimag( Gmn[ (Sampling_on_torus_phi - n)  ]  ) );   

  //printf("fuck\n");

}

//printf("==========================case n=0=========================\n");     

// case n=0
//initializing m and l!
m=0;
int Power1= pow(-1,-s);
while( m<=lmax ){ // DD is zero for l+m=odd!!!
    m+=1; // raising the n!  
    // l= Sup(abss ,m);  
    l=m;
    while(l<=lmax){
        l_centro= l*(l+1); 

        
       // printf("l=%i,m=%i, n=%i, \n ",l,m,0 );
       // printf("DD[1][%i][%i][%i][%i]=%f+%fi \n ", s_index,l,m, 0, creal( DD[1][l][ s_index ][m][0] ) , cimag( DD[1][l][  s_index ][m][0] ) ); 

        if ( Tolerance <  cabs( DD[1][l][ s_index ][m][0]  ) ){  

        // (m+)
        Gmn[ m*Sampling_on_torus_phi ] += DD[1][l][ s_index ][ m ][0] * aln[ l_centro ];  
        //Gmn[ m*Sampling_on_torus_phi ] +=  aln[ l_centro ];                        
                             
        // (m-)
        // Gmn[ (Sampling_on_torus_theta-m)* Sampling_on_torus_phi  ] = Power1 * Gmn[ m*Sampling_on_torus_phi ];
        Gmn[ (Sampling_on_torus_theta-m)* Sampling_on_torus_phi  ] = Power1 * Gmn[ m*Sampling_on_torus_phi ];

      }
       // printf("l,m,n, index =%i,%i,%i,%i\n",lo2,mo,0, mo*Sampling_on_torus_phi  ); 
       // printf("DD[1][%i][%i][%i][%i]=%f \n ",index_s,lo2,mo, 0, creal( DD[1][ index_s ][lo2][mo][0] ) ); 
       // printf("----------\n");
    


       l+=2; // DD is zero for l+n=odd!!!
       }    
//printf("Gmn[%i] = %f + %f i \n", m*Sampling_on_torus_phi , creal( Gmn[m*Sampling_on_torus_phi] ) , cimag( Gmn[ m*Sampling_on_torus_phi ] ) );   
//printf("Gmn[%i] = %f + %f i \n", (Sampling_on_torus_theta-m)* Sampling_on_torus_phi  , creal( Gmn[ (Sampling_on_torus_theta-m)* Sampling_on_torus_phi ]) , cimag( Gmn[ (Sampling_on_torus_theta-m)* Sampling_on_torus_phi ]  ) );   

 }

//printf("==========================others=========================\n");

//initializing m,l!
m=0;

while(m <= lmax){ 
   m+=1; // raising the m!  
   n=0; //initializing n
   while(n<= lmax){
    n+=1; // raising the n!  
    
    l= Sup( Sup( abss ,n ) , Sup( abss , m ) );

    while(l<=lmax){ //l start from value of n. It cannot be lower than n,m or s!!! 
        l_centro= l*(l+1); 
       //////////////////////////////////////
   
        if ( Tolerance <  cabs( DD[1][l][ s_index ][m][n]  ) ){ 

        // (m+,n+)
        Gmn[ m*Sampling_on_torus_phi +  n ] += DD[1][l][ s_index ][m][n] * aln[ l_centro + n ];
        //Gmn[ m*Sampling_on_torus_phi +  n ] +=  aln[ l_centro + n ];

        // (m+,n-)
        Gmn[ m*Sampling_on_torus_phi + (Sampling_on_torus_phi - n) ] += DD[1][l][ s_index ][m][ 2*l+1-n ] * aln[ l_centro - n ];         
        //Gmn[ m*Sampling_on_torus_phi + (Sampling_on_torus_phi - n) ] +=   aln[ l_centro - n ];         

 
       // printf("l,m,n, index =%i,%i,%i,%i\n",l,m,n, m*Sampling_on_torus_phi +  n ); 
       // printf("DD[1][%i][%i][%i][%i]=%f \n ",s_index,l,m, n, creal( DD[1][l][ s_index ][m][n] ) ); 
       // printf("DD[1][%i][%i][%i][%i]=%f \n ",s_index,l,m, 2*l+1-n, creal( DD[1][l][ s_index ][m][2*l+1-n] ) ); 
       // printf("----------\n");
       // }

        }

        l+=1; //raising the l!!!        
        }
       
        // Using mirror for properties for m-
        //(m-,n+)
        Gmn[ (Sampling_on_torus_theta-m)* Sampling_on_torus_phi + n ] = pow(-1,-s-n) * Gmn[ m*Sampling_on_torus_phi + n ];
       //(m-,n-)
        Gmn[ (Sampling_on_torus_theta-m)* Sampling_on_torus_phi + (Sampling_on_torus_phi-n) ] = pow(-1,-s+n) * Gmn[ m*Sampling_on_torus_phi + (Sampling_on_torus_phi -n) ];
   
       // printf("DD[1][%i][%i][%i][%i]=%f \n ",index_s,l,m, n, cimag( DD[1][ index_s ][l][m][n] ) );  
       //printf("DD[1][%i][%i][%i][%i]=%f \n ",index_s,l,m,index_s - n, cimag(DD[1][ index_s ][l][m][ index_s-n ] ) );      
       // printf("!Gmn[%i] = %f + %f i \n", m*Sampling_on_torus_phi +  n , creal( Gmn[ m*Sampling_on_torus_phi +  n ] ) , cimag( Gmn[ m*Sampling_on_torus_phi +  n ] ) );   
       // printf("!Gmn[%i] = %f + %f i \n",m*Sampling_on_torus_phi + (Sampling_on_torus_phi - n)  , creal( Gmn[ m*Sampling_on_torus_phi + (Sampling_on_torus_phi - n) ]) , cimag( Gmn[ m*Sampling_on_torus_phi + (Sampling_on_torus_phi - n) ]  ) );   
       // printf("!Gmn[%i] = %f + %f i \n", (Sampling_on_torus_theta-m)* Sampling_on_torus_phi + n ,  creal( Gmn[ (Sampling_on_torus_theta-m)* Sampling_on_torus_phi + n ] ), cimag( Gmn[ (Sampling_on_torus_theta-m)* Sampling_on_torus_phi + n ] ) );   
       // printf("!Gmn[%i] = %f + %f i \n",(Sampling_on_torus_theta-m)* Sampling_on_torus_phi + (Sampling_on_torus_phi-n)  , creal( Gmn[ (Sampling_on_torus_theta-m)* Sampling_on_torus_phi + (Sampling_on_torus_phi-n)  ] ),cimag( Gmn[ (Sampling_on_torus_theta-m)* Sampling_on_torus_phi + (Sampling_on_torus_phi-n)  ] ) );   
    
     }
   }  
 }  

//-----------------------------------------------------
void spin_backward_transform(fftw_complex * restrict f, fftw_complex * restrict aln, int s, int lmax , int smax, int Ntheta, int Nphi, fftw_complex ***** restrict DD ){
  //Sampling of the extended function : HW
  
  int Sampling_on_torus_theta= 2*(Ntheta-1);
  int Sampling_on_torus_phi = Nphi; //sim incluir el 2pi!!!!
  
  if ( 2*lmax+1 > Sampling_on_torus_phi ) {
        printf("Number of points in phi are less than 2*lmax+1 (read the manual!)\n");
        printf("Abort swsh-backward-transform!\n");
        return ;
  }
  if ( 2*lmax+1 > Sampling_on_torus_theta ){
        printf("Number of points in theta are less than lmax+3/2 (read the manual!)\n");
        printf("Abort swsh-backward-transform!\n");
        return ;
  }

  //Creating pointer for Gmo
  fftw_complex *Gmn = calloc( Sampling_on_torus_theta * Sampling_on_torus_phi , sizeof(fftw_complex) );
  //fftw_complex *Gmn =  fftw_malloc( Sampling_on_torus_theta*Sampling_on_torus_phi*sizeof(fftw_complex));
 

//  clock_t begin = clock();
  //Computing Gmo
  Compute_Gmn( aln, Gmn, DD,  s, lmax, smax, Sampling_on_torus_theta, Sampling_on_torus_phi );

// clock_t end = clock();
//  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
 


 // clock_t begin2 = clock();

  //Function on the torus
  //fftw_complex *f_torus = calloc( Sampling_on_torus_theta*Sampling_on_torus_phi  , sizeof(fftw_complex) );
  fftw_complex *f_torus = fftw_malloc( Sampling_on_torus_theta * Sampling_on_torus_phi *sizeof(fftw_complex)); 


  //Computing the backward fourier transform
  fftw_plan fftplan = fftw_plan_dft_2d( Sampling_on_torus_theta, Sampling_on_torus_phi , Gmn, f_torus, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(fftplan);
  fftw_destroy_plan(fftplan);

   //Passing the function from the torus to the sphere
   int itheta;
   int iphi;

   for (itheta = 0; itheta < Ntheta; itheta++){ 
    for (iphi = 0; iphi < Nphi; iphi++){

      f[ itheta * Sampling_on_torus_phi + iphi ] =  f_torus[ itheta * Sampling_on_torus_phi + iphi ]  ;     
      //printf("f[%i] = %f fuck \n", itheta * Sampling_on_torus_phi + iphi , creal(f[itheta * Sampling_on_torus_phi + iphi]) );  
    }
    //printf("--------------------------------------------------------\n");
  }

free(f_torus);
free(Gmn);   

//clock_t end2 = clock();
//double time_spent2 = (double)(end2 - begin2) / CLOCKS_PER_SEC;
//printf("time_back2!!!! =%f \n ", time_spent2 ); 

}  




