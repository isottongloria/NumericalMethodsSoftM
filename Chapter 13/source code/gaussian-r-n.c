
#include <sys/time.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "RngStream.h"
#include "array_alloc.h"
#include "vectors2d.h"
#include "gaussian-r-n.h"
#include "sim-gauss-2.h"

#define two_pi 2.0 * M_PI
#define MAXG 10000000

extern RngStream rngs; // setup in the main C

int iG=0,iGnew=1;
//double **vGauss;
double vGauss[2][10000000];

int Nr, Ns, Nturn, icycle=0, iturn=0, irow=0;
double ***rot;
int **perm;
int *check,N_ok=0;

/* Marsaglia (1972) */
void Gaussian_(double *y1, double *y2)
{
  double x1, x2, w;
  do {
    x1 = 2.0 * RngStream_RandU01(rngs) - 1.0;
    x2 = 2.0 * RngStream_RandU01(rngs) - 1.0;
    w = x1 * x1 + x2 * x2;
  } while ( w >= 1.0 );
  
  w = sqrt( (-2.0 * log( w ) ) / w );
  *y1 = x1 * w;
  *y2 = x2 * w;
}


/* Box-Muller alg. from Wikipedia */
void BoxMuller(double *y1, double *y2)
{
  //create two random numbers, make sure u1 is greater than epsilon
  double u1, u2, mag,phi;
  do{ u1 = RngStream_RandU01(rngs);}while(u1 <= 1.e-20);
  u2 = RngStream_RandU01(rngs);
  //compute two Gaussian random numbers
  mag = sqrt(-2.0 * log(u1));
  phi = two_pi * u2;
  *y1 = mag * cos(phi);
  *y2 = mag * sin(phi);
  return;
}


void Fresh_Gaussian(){
  double g1,g2;
  int i;
  /* fill the vector with fresh Gaussian numbers */
  for(i=0;i<Ns;i++){
    BoxMuller(&g1, &g2);
    vGauss[0][i] = g1;
    vGauss[1][i] = g2;
  }
}


void Wallace_Gaussian(){
  double g1,g2;
  int FRESH,i,i1,i2,j1,j2;
  for(i1=0;i1<Ns;i1+=2){
    i2=i1+1;
    j1 = perm[irow][i1];
    j2 = perm[irow][i2];
    vGauss[iGnew][i1] = vGauss[iG][j1]*rot[irow][i1][0]+vGauss[iG][j2]*rot[irow][i1][1];
    vGauss[iGnew][i2] = vGauss[iG][j1]*rot[irow][i2][1]+vGauss[iG][j2]*rot[irow][i2][1];
  }
  FRESH=0;
  irow++;
  if(irow==Nr){
    irow=0;
    iturn++;
    if(iturn==Nturn){
      iturn=0;
      icycle++;
      /* time for a fresh Gaussian set */
      FRESH=1;
      Fresh_Gaussian();
    }
  }
  /* unless a fresh Gaussian set was drawn, copy the new Gaussians to vGauss */
  if(iG){iG=0; iGnew=1;} else{ iG=0; iGnew=0; }
  if(icycle==0 && iturn==0 && irow<6){
    fprintf(stderr,"\nWallace: %lf %lf .. %lf %lf",
	    vGauss[0][0],vGauss[0][1],vGauss[0][Ns-2],vGauss[0][Ns-1]);
    fprintf(stderr,"\nWallace: %lf %lf .. %lf %lf",
	    vGauss[1][0],vGauss[1][1],vGauss[1][Ns-2],vGauss[1][Ns-1]);
    fprintf(stderr,"\n");
  }
  
  return;
}

void Shuffle(int *v, int N){
  int i,j,ii;
  for(i=N-1;i>=0;i--){
    ii = RngStream_RandInt1(rngs,N);
    j=v[i];
    v[i]=v[ii];
    v[ii]=j;
  }
}

/* TO-DO: The Wallace Gaussian Generator, see
https://developer.nvidia.com/gpugems/gpugems3/part-vi-gpu-computing/chapter-37-efficient-random-number-generation-and-application
*/
void Set_Wallace_method(int nslot, int nrow, int nturn){
  int i,j,k,ii;
  double phi;
  Nr=nrow;
  Ns=nslot;
  Ns+=(Ns%2); // even number, Ns
  rot = (double ***)dmatrix3(Nr,Ns,2);
  perm=(int **)imatrix(0,Nr-1,0,Ns-1);
  check=(int *)ivector(0,Ns-1);
  //vGauss=(double **)dmatrix(0,1,0,Ns-1);
  if(Ns>MAXG){
    fprintf(stderr,"\nProblem: Ns=%d is larger than MAXG\n");
    exit(3);
  }

  /* nrow Permutations */
  for(k=0;k<Nr;k++) for(i=0;i<Ns;i++) perm[k][i]=i;
  for(k=0;k<Nr;k++){
    
    /* shuffle 3 times */
    Shuffle(perm[k], Ns);
    Shuffle(perm[k], Ns);
    Shuffle(perm[k], Ns);

    /* Check */
    N_ok=0;
    for(i=0;i<Ns;i++) check[i]=0;
    for(i=0;i<Ns;i++){
      if(check[perm[k][i]] && perm[k][i]){
	fprintf(stderr,"\nERROR 1\n");
	fprintf(stderr,"\nk=%d i=%d perm[k][i]=%d",k,i,perm[k][i]);
	exit(1);}
      check[perm[k][i]]=1;
      N_ok++;
    }
    if(N_ok!=Ns){ fprintf(stderr,"\nERROR 2 , %d\n",N_ok); exit(2);}
  }
  k=0;
  fprintf(stderr,"Example of perm[]: %d %d %d ... %d %d %d (% values)\n",
	  perm[k][0],perm[k][1],perm[k][2],perm[k][Ns-3],perm[k][Ns-2],perm[k][Ns-1],Ns);

  /* nrow 2d-rotations */
  for(k=0;k<Nr;k++){
    for(i=0;i<Ns;i+=2){
      phi=two_pi*RngStream_RandU01(rngs);
      rot[k][i][0] = cos(phi);  rot[k][i+1][0] = sin(phi);
      rot[k][i][1] = -sin(phi); rot[k][i+1][1] = cos(phi);
    }
  }
  k=0;
  fprintf(stderr,"Example of rot[]: %8.4lf %8.4lf ... %8.4lf %8.4lf\n",
	  rot[k][0][0],rot[k][1][0],rot[k][Ns-2][0],rot[k][Ns-1][0]);
  fprintf(stderr,"                : %8.4lf %8.4lf ... %8.4lf %8.4lf\n",
	  rot[k][0][1],rot[k][1][1],rot[k][Ns-2][1],rot[k][Ns-1][1]);

  Fresh_Gaussian();
}
