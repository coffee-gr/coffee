//////////////////////////////////////////////////////
// Dr. Leon Escobar
// Departmento de Matematicas,
// Universidad del Valle,
// Colombia
// 2019
//////////////////////////////////////////////////////

#include <eths.h>

///// eth Up /////
void spin_ethU(fftw_complex *f, fftw_complex *alm, int s, int lmax , int smax, int Ntheta , int Nphi,  fftw_complex *****DD, fftw_complex *W ){

spin_forward_transform( f , alm, s,  lmax, smax, Ntheta, Nphi, DD, W); 
int Nlm = lmax*(lmax+2)+1 ; // number of alm per lmax!
int m=0; //contador
int l=0;
int i=0;
//float r;

double square_root_terms_up ;//(for speed up)

while( i < Nlm ){

               square_root_terms_up = -sqrt( ( l-s )*( l+s+1 ) );

               while(i <= m){ 
               //printf("antes aln[%i]=%f \n ",i,creal(alm[i] ) );  
               if(l < abs(s)){ 
                            alm[i] = 0.0;                            
                            }
               else{ 
               alm[i] = square_root_terms_up * alm[i];  
                   } 
               i=i+1;        
               }

               //printf("despues aln[%i]=%f \n ",i,creal(alm[i] ) );   
               l = l+1;
               m = ( 2*l + 1 ) + m ;
   }
  s = s + 1 ; 
  spin_backward_transform( f, alm, s, lmax, smax, Ntheta, Nphi, DD );
}



///// eth Down /////
void spin_ethD(fftw_complex *f, fftw_complex *alm, int s, int lmax, int smax, int Ntheta , int Nphi,  fftw_complex *****DD, fftw_complex *W ){
  
spin_forward_transform(f, alm, s,  lmax, smax, Ntheta, Nphi, DD, W); 
int Nlm = lmax*(lmax+2)+1 ; // number of alm per lmax!
int m=0; //contador
int l=0;
int i=0;

double square_root_terms_down ;//(for speed up)

//float r;
while( i < Nlm ){

               square_root_terms_down = sqrt( ( l+s )*( l-s+1 ) );

               while(i <= m){  
               //printf("antes aln[%i]=%f \n ",i,creal(alm[i] ) ); 
               if( l < abs(s) ){ 
                        alm[i] = 0.0;
                        } 
               else{ 
               alm[i] =  square_root_terms_down * alm[i];
                    } 
               i=i+1; 
               }
               //printf("despues aln[%i]=%f \n ",i,creal(alm[i] ) );   
               l = l+1;
               m = ( 2*l + 1 ) + m ;
   }
  s = s - 1 ; 
  spin_backward_transform( f, alm, s, lmax ,smax, Ntheta, Nphi, DD );

}
/////////////////////////


