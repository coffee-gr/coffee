
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

#include <spin_forward_transform.h>
#include <time.h>


void Compute_Imn_by_extending_f_to_torus( fftw_complex * restrict f, fftw_complex * restrict Imn, fftw_complex * restrict W, int s,  int lmax, int Sampling_on_torus_theta, int Sampling_on_torus_phi ){
  // F is complex version of f, extended to the whole sphere Works for pixelization where first ring is north pole and the last ring is the south pole
  // i.e., theta = itheta*pi/(Ntheta-1). The number of theta rows in F, the extended version of f, is 2*(Ntheta-1) 
  // Extend f with F(theta > pi) = (-1)^s f(2pi - theta). Finally we return the fft of F
 
  int itheta;
  int iphi;

  double Norm = 2 * M_PI/ (Sampling_on_torus_theta * Sampling_on_torus_phi); // =  2pi/Ntheta_extended 
  // Building the extended function F and multipliying by the real function W and the Norm. 
  fftw_complex *F = fftw_malloc(Sampling_on_torus_theta*Sampling_on_torus_phi*sizeof(fftw_complex)); 

  // Extending the function f in S2 to the T2
  int Ntheta =  Sampling_on_torus_theta/2 + 1;
  int opp_iphi;

  double power_s= pow(-1.0,s) ;//(speed up)

  for (itheta = 0; itheta < Ntheta; itheta++) { 

     for (iphi = 0; iphi < Sampling_on_torus_phi; iphi++) {          

        opp_iphi=(iphi + Sampling_on_torus_phi/2) % Sampling_on_torus_phi;
        
        F[itheta * Sampling_on_torus_phi  + iphi] =  creal( W[itheta] ) * f[itheta * Sampling_on_torus_phi + iphi] * Norm ;  
                
        // Trick for teh extention of HW
        // int opp_iphi = (iphi + Sampling_on_torus_phi/2) % Sampling_on_torus_phi;
        // printf(" opp_iphi =%i , iphi=%i ",opp_iphi,iphi);             
        // Trick for the extention of McEwen & Wiaux
        //int m = (iphi <= Sampling_on_torus_phi/2) ? iphi : (iphi - Sampling_on_torus_phi); // convert index of m to m value
      
        if( (itheta > 0) && (itheta < Ntheta) ){
      
        F[ (Sampling_on_torus_theta - itheta)* Sampling_on_torus_phi + opp_iphi ] = power_s * creal(W[Sampling_on_torus_theta - itheta]) * f[ itheta*Sampling_on_torus_phi + iphi ] *  Norm;          
       }
    }
 }  
  //Applying the 2D-transform 
  fftw_plan fftplan = fftw_plan_dft_2d( Sampling_on_torus_theta, Sampling_on_torus_phi,  F , Imn , FFTW_FORWARD, FFTW_ESTIMATE );   
  fftw_execute(fftplan);
  fftw_destroy_plan(fftplan); 
  free(F);
}

//-----------------------------------------------------
void Compute_Jmn( fftw_complex * restrict f, fftw_complex * restrict Jmn, fftw_complex * restrict W, int s,int lmax,  int Sampling_on_torus_theta, int Sampling_on_torus_phi ) { 
  //fftw_complex *Imn = calloc( Sampling_on_torus_theta * Sampling_on_torus_phi , sizeof(fftw_complex*) );  
  fftw_complex *Imn =  fftw_malloc(sizeof(fftw_complex) * Sampling_on_torus_theta * Sampling_on_torus_phi); 

  //Compute_Imn_by_extending_f_to_torus_convolution(f,Jmn,W,s,lmax,Ntheta,Nphi,Sampling_on_torus_theta,Sampling_on_torus_phi);
  Compute_Imn_by_extending_f_to_torus( f , Imn , W , s , lmax, Sampling_on_torus_theta, Sampling_on_torus_phi ); 
 
  //Building the Jmn terms  ( m just takes positives values , So the size is just: lmax + 1  ) 
  int m,n;  
  int Nn = 2*lmax+1;
  int Nm = lmax;

 double power_s_plus_n;//(speed up)
 double power_s_minus_n;//(speed up)

 for (n=0; n <= Nm ; n++){ 

  power_s_plus_n = pow(-1.0,s+n);
  power_s_minus_n = pow(-1.0,s-n);

  for (m=0; m <= Nm ; m++){      
    //printf("m=%i, n=%i,\n", m , n );  
    // cases for n!
    if (m==0){
        Jmn[n]= Imn[n] ;    
        //printf("Jmn[%i] = %f\n", n , creal( Jmn[n] ) );   

        // filling n-
        if ( 0<n && n<lmax){
        Jmn[ Nn - n ] = Imn[ Sampling_on_torus_phi - n ];
        //printf("Jmn[%i] = %f\n",  Nn - n , creal( Jmn[ Nn - n] ) );   
        }      
     }
    else{ //Recall that from the FFWT frecuencies{0,1,.....n,-n,....-1}!!!!!!!
    // filling n+!!! 
    Jmn[ m * Nn + n ] = Imn[ m * Sampling_on_torus_phi +  n ] + power_s_plus_n * Imn[  (Sampling_on_torus_theta - m) * Sampling_on_torus_phi + n]; // Where Imo[-m] = Imo[ 2*lmax+1 + m ]   
    //printf("Jmn[%i] = %f\n", m * Nn + n , creal( Jmn[m * Nn + n] ) );     
   
     // filling n-!!!
    if (0<n && n<=lmax){
    Jmn[ m * Nn + (Nn- n) ] = Imn[ m * Sampling_on_torus_phi +  (Sampling_on_torus_phi - n ) ] +  power_s_minus_n * Imn[ ( Sampling_on_torus_theta - m ) * Sampling_on_torus_phi+  (Sampling_on_torus_phi - n  )]; // Where Imo[-m] = Imo[ 2*lmax+1 + m ]   
    //printf("Jmn[%i] = %f\n", m * Nn + (Nn- n) , creal( Jmn[m * Nn + (Nn- n)] ) );    
 
    }
   }
  }
 } 

fftw_free(Imn);  
} 
//-----------------------------------------------------
// Here we denote the aln and not alm!!!! 
// Note that f is flatten!!!!....following the row-major-ordering convention

void spin_forward_transform(fftw_complex * restrict f, fftw_complex * restrict aln, int s, int lmax, int smax, int Ntheta , int Nphi, fftw_complex ***** restrict DD ,  fftw_complex * restrict W ){
 
  //Sampling of the extended function : HW
  int Sampling_on_torus_theta = 2*(Ntheta-1);
  int Sampling_on_torus_phi = Nphi; //sim incluir el 2pi!!!!

   // Note that one condition for the transform to make sence is that both 
  // Sampling_on_phi, Sampling_on_torus > 2* lmax + 1 .
  // That is why we imporse the restriction as follows:  

  if ( 2*lmax+1 > Sampling_on_torus_phi ) {
        printf("Number of points in phi are less than 2*lmax+1 (read the manual!)\n");
        printf("Abort swsh-forward-transform\n");
        return ;
  }
  if ( 2*lmax+1 > Sampling_on_torus_theta ){
        printf("Number of points in theta are less than lmax+3/2 (read the manual!)\n");
        printf("Abort swsh-forward-transform\n");
        return ;
  }

//Creating pointer matrix for Jmn with row-major-ordering 
fftw_complex *Jmn = calloc(  (lmax+1) * (2*lmax+1) , sizeof(fftw_complex) );
//fftw_complex *Jmn = fftw_malloc((lmax+1)*(2*lmax+1)*sizeof(fftw_complex));


//computing Jmn

//clock_t begin2 = clock();
Compute_Jmn(f, Jmn, W, s, lmax, Sampling_on_torus_theta, Sampling_on_torus_phi);
//clock_t end2 = clock();
//double time_spent2 = (double)(end2 - begin2) / CLOCKS_PER_SEC;
//printf("time_spent_forward1=%f \n ", time_spent2 ); 


int l,n, m ; 
int number_of_n = (2*lmax+1);//(speed up)
int l_centro;
int s_index;//(speed up)

//////////////////////////

//clock_t begin = clock();

s_index = pindex(smax,s); 


double Tolerance = pow(10,-13);

for (l=abs(s);l<=lmax;l++) {  

    l_centro = l*(l+1);
  
        // for n=0. This case can be simplified by half!!!
      for (m = l%2 ; m<=l ; m+=2 ){//THIS IS AN EXTRA SIMPLIFICATIOON BECAUSE D[l][m][0] = 0  if l+m = odd. 
        
            if ( Tolerance <  cabs( DD[0][l][ s_index ][m][0]  ) ){ 

           aln[ l_centro ] += DD[0][l][ s_index ][m][0] * Jmn[ m* number_of_n ] ;  
           //printf("s=%i,l=%i,m=%i, n=%i, \n ",s,l, m, 0);
           //printf("!DD[0][%i][%i][%i][0]=%f+%fi \n ",s_index,l,m,0, creal( DD[0][l][ s_index ][m][0])+ cimag( DD[0][l][ s_index ][m][0])  );  
           //printf("Jmn[ m* number_of_n ]=%f+%fi\n", creal(Jmn[ m* number_of_n ])+ cimag(Jmn[ m* number_of_n ]) ) ;
         }


          }
     // for the rest n>1
      for ( n =1 ; n<= l; n++ ){

           for (m =0 ; m<=l; m++ ){   

              if ( Tolerance <  cabs( DD[0][l][ s_index ][m][n]  ) ){
            
               // n+
               aln[ l_centro + n ] += DD[0][l][ s_index ][m][n] * Jmn[ m* number_of_n  +  n] ; 
               // n-
               aln[ l_centro - n ] += DD[0][l][ s_index ][m][ 2*l+1 - n ] * Jmn[ m* number_of_n + ( number_of_n - n ) ] ; 


               }
               //printf("s=%i,l=%i,m=%i, n=%i, \n ",s,l, m, n );
               //printf("DD[0][%i][%i][%i][%i]=%f+%fi \n ",s_index,l,m,  n,  creal( DD[0][l][ s_index ][m][  n] ) + cimag( DD[0][l][ s_index ][m][  n] ) );  
               //printf("Jmn[ m*(2*lmax+1)  +  n]=%f+%fi\n", creal(Jmn[ m* number_of_n + (    n ) ])+ cimag(Jmn[ m* number_of_n + (   n ) ]) ) ;
               //printf("fuck aln[%i %i]=%f \n ",l, n,creal(aln[ l*(l+1) + n ] ) ); 
               
               //printf("s=%i,l=%i,m=%i, n=%i, \n ",s,l, m, number_of_n - n );
               //printf("DD[0][%i][%i][%i][%i]=%f+%fi \n ",s_index,l,m,number_of_n - n,  creal( DD[0][l][ s_index ][m][ 2*l+1 - n ] ) + cimag( DD[0][l][ s_index ][m][2*l+1 - n] ) );  
               //printf("Jmn[ m*(2*lmax+1)  +  n]=%f+%fi\n", creal(Jmn[ m* number_of_n + ( number_of_n - n ) ])+ cimag(Jmn[ m* number_of_n + ( number_of_n - n ) ]) ) ;
               //printf("fuck aln[%i %i]=%f \n ",l, number_of_n - n ,creal(aln[ l*(l+1) - n ] ) ); 
               //printf("--------------------------\n" );  
              //printf("s=%i,l=%i,m=%i, n=%i, \n ",s,l, m,number_of_n - n ); 
              //printf("DD[0][%i][%i][%i][%i]=%f \n ",index_s,l,m,n, DD[0][ index_s ][l][m][number_of_n - n] );  
         }
      } 
   }

free(Jmn);

//clock_t end = clock();
//double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
//printf("time_spent_fordward2 =%f \n ", time_spent ); 

}  

