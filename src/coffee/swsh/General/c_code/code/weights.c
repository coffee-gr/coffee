//////////////////////////////////////////////////////
// Dr. Leon Escobar
// Department of Mathematics,
// Universidad del Valle,
// Colombia
// 2019
//////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////
// 
// This two functions  are taken from  Kevin M. Huffenberger & Benjamin D. Wandelt code
// This functions is part of spinsfast.
// 
////////////////////////////////////////////////////////////////////////////////////////

#include <weights.h>

//-----------------------------------------------------
int indx2p (int ip, int wsize) {
  return( (ip > wsize/2) ?  ip-wsize : ip );
}


//-----------------------------------------------------
void Compute_quadrature_weights( fftw_complex * restrict W, int wsize ) {  
  // Creating weights in the fourier space
  fftw_complex *w = calloc(wsize, sizeof(fftw_complex)); 
   int ip,p, eo;
  for (ip=0; ip<wsize;ip++) {
    p = indx2p(ip,wsize);    
    eo = abs(p % 2);
    if (p == -1) {
      w[ip] = I * M_PI/2.;
    } else if (p == 1) {
      w[ip] = - I * M_PI/2.;
    } else if ( eo == 0) {
      w[ip] = 2./(1.-p*p);
    } else {
      w[ip] = 0;
    }    
  }  
  fftw_plan wplan = fftw_plan_dft_1d(wsize, w, W, FFTW_BACKWARD, FFTW_ESTIMATE);
 
  fftw_execute(wplan);
  fftw_destroy_plan(wplan);
  free(w);
}



//-------------------s----------------------------------
//Precompute the quadrature weights
fftw_complex *Precompute_quadrature_weights( int Ntheta ){
  //Sampling of the extended function : HW
  int Sampling_on_torus= 2*(Ntheta-1); 
  fftw_complex *W = calloc( Sampling_on_torus , sizeof(fftw_complex)); 
  Compute_quadrature_weights( W , Sampling_on_torus );
  return W;  
}