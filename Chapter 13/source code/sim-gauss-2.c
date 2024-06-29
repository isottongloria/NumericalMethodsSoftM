/* 
   Marco Baiesi
   December 2023
   class of the course
   Numerical Methods in Soft Matter
   for the masters Physics and Physics of Data at the University of Padova
*/
#include <sys/time.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "RngStream.h"
#include "gaussian-r-n.h"
#include "array_alloc.h"
#include "vectors2d.h"
#include "useful-tools.h"
#include "sim-gauss-2.h"

#define D 2
#define MAX 21000000

//#define DEB  // if defined, some debug printing is enabled

/* VARIABLES FOR RANDOM NUMBERS */
unsigned int iseed;
struct timeval tp;
struct timezone zp;
/******* RngStream********/
RngStream rngs;
unsigned long sseed;
unsigned long vseed[6]={1023,12009,999143,71347,96663,3213};
/******* RngStream********/
time_t time_start,time_now; 
double elapsed;


/////
double T = 0.5;
double dt = 5.e-3, Dt = 0.1;
double tt=20000, ta; // max time, total accumulated time
double mu, gam = 1., mu_dt, mu0, gam0 = 20., mu0_dt;

double eps=1., eps0=10.;  // Mathematica: gauss-pot.png
double R=0.125, R0=1.25;
double sigma, sigma0;
double fp, fp0;
double cut, cut0, cutSQ, cut0SQ;
double inv, inv0, stoch, stoch0;
double x_min, x_max, y_min, y_max, y0_min, y0_max;

int Dv0[10000][2], N0=0;

int sim_time_T=0; // global counter of the MD steps
int N = 1000.;  // number of small particles
double *box, half_box[D];

double v_trap, v_trap_ini, v_per_decade, EF;
int Nv;
int nstep_x0;     // =Dt/dt, to perform undersampling of the probe's position
int nstep_save=0; // how many Dt between config. saving, no saving if =0
double k_trap=0.5, *x_trap;
double Lambda;
double f_active;

/* polymer length (Lp), number (Mp) and bond strength (k_pol) */
int Lp,Mp;
double k_pol;
int i_ini, i_fin; // 1, N
// coords., forces, noises
double **x, **F, **dB;
// counting +1 each time the box right boundary is crossed and -1 for the left one 
int **travel;
// cell list
int **ix, **tau_cell, **phi_cell, nu_list[MAX], Dim_a, Dim_b;
// bonds
int *BONDED;
double **vb;

// time series
double *x_rel;

int nrun=0;
char fnO[1000], str0[1000];
int FIRST_SAVE=1, THERMALIZING=1;

char dirDATA[1000],dirCONF[1000];

/* === Gaussian nrs === */
double gauss1, gauss0;
int igau = 0;
/* ==================== */


void Print_elapsed_time(){
  int h,m,s;
  time(&time_now);
  elapsed = difftime(time_now,time_start) / 60.;
  m = (int)floor(elapsed);
  s = (int)floor((elapsed - m)*60.);
  printf("[%d:%d]",m,s);
  return ;
}

void Print_info(FILE *out1)
{
	int i;
	fprintf(out1,"\n# =================================");
	fprintf(out1,"\n#    tt      = %.0lf",tt);
	fprintf(out1,"\n#    dt      = %lf",dt);
	fprintf(out1,"\n#    Dt      = %lf",Dt);
	fprintf(out1,"\n#    T       = %lf",T);
	fprintf(out1,"\n#    v_min   = %lf",v_trap_ini);
	fprintf(out1,"\n#    Nv      = %d",Nv);
	fprintf(out1,"\n#    v_x_dec = %.0lf",v_per_decade);
	fprintf(out1,"\n#    mu0     = %lf",mu0);
	fprintf(out1,"\n#    mu      = %lf",mu);
	fprintf(out1,"\n#  . R0      = %lf",R0);
	fprintf(out1,"\n#  . eps0    = %lf",eps0);
	fprintf(out1,"\n#  . sigma0  = %lf",sigma0);
	fprintf(out1,"\n#  . cut0    = %lf",cut0);
	fprintf(out1,"\n#  . fp0     = %lf",fp0);
	fprintf(out1,"\n#  . inv0    = %lf",inv0);
	fprintf(out1,"\n#    R       = %lf",R);
	fprintf(out1,"\n#    eps     = %lf",eps);
	fprintf(out1,"\n#    sigma   = %lf",sigma);
	fprintf(out1,"\n#    cut     = %lf",cut);
	fprintf(out1,"\n#    fp      = %lf",fp);
	fprintf(out1,"\n#    inv     = %lf",inv);
	fprintf(out1,"\n#    k_trap  = %lf",k_trap);
	fprintf(out1,"\n#    k_pol   = %lf",k_pol);
	fprintf(out1,"\n#  * Dim_a   = %d",Dim_a);
	fprintf(out1,"\n#  * Dim_b   = %d",Dim_b);
	fprintf(out1, "\n#   Lambda  = %lf", Lambda);
	fprintf(out1, "\n#   f_act   = %lf", f_active);
	fprintf(out1, "\n# =================================\n");
  	//getchar();
}

void Make_str0(char *s1){
  sprintf(str0,"N%d-n%d_%.0lfx%.0lf_T%.2lf_k%.1lf-%.1lf_r%.3lf-%.1lf_R%.3lf-%.1lf_dt%.4lf",
	  N,Lp,box[0],box[1],T,k_trap,k_pol, R,eps,R0,eps0,dt);
	  
  strcat(str0,s1);
}



/* https://gcc.gnu.org/onlinedocs/gcc/Inline.html :
   By declaring a function inline, you can direct GCC to make calls to that function faster. One way GCC can achieve this is to integrate that function's code into the code for its callers. This makes execution faster by eliminating the function-call overhead */

inline double Noise()
{
	if (igau)
	{
		igau = 0;
		return gauss1;
	}
	else
	{
		BoxMuller(&gauss0, &gauss1);
		igau = 1;
		return gauss0;
	}
}


void MD_step(){
  double r,r2,v[D],e,force;
  int i,j,ii, a,b,Da,Db,a_p=-1,b_p=-1;
  
	/*
	  -----------------------------------------------------
	  big particle: 0
	  small particles: i_ini=1 --> i_fin=N
	  trap position: N+1 (when saving on file), x_trap here
	  -----------------------------------------------------
	*/

	/*           */
	/* Cell list */
	/*           */

	/*
	  sim_time_T has such a long periodicity that a cell is
	  tagged for sure at least once during the period;
	  hence, there is no danger that tau_cell[a][b] = sim_time_T
	  of 20000000 iterations ago
	*/
  sim_time_T = (sim_time_T+1)%20000000;
  
  for(i=i_ini;i<=i_fin;i++){
    for(j=0;j<D;j++) ix[i][j] = (int)floor(x[i][j]);
    a=ix[i][0]; b=ix[i][1];
    if(tau_cell[a][b] != sim_time_T){
      /* start a list for this cell */
      tau_cell[a][b]=sim_time_T;
      phi_cell[a][b]=i;
      nu_list[i] = -1; // -1: last element of the list of cell [a][b]
    }
    else{
      /* add at the beginning the of the cell list */
      nu_list[i] = phi_cell[a][b];
      phi_cell[a][b]=i;
    }
  }

  
  /* deterministic dynamics of the trap */
  x_trap[0] = x_trap[0] + v_trap * dt;

  /* random noise terms, and initialization of forces */
  for(j=0;j<D;j++){
    dB[0][j] = stoch0 * Noise();
    F[0][j]  = k_trap * (x_trap[j] - x[0][j]);
    for(i=i_ini;i<=i_fin;i++){
      dB[i][j] = stoch * Noise();
      F[i][j] = 0.;
    }
  }

  if(eps0!=0){
    /* big particle (0) interacting with smaller ones sitting in list of (relative) cells "Dv0" */
    a=(int)floor(x[0][0]);
    b=(int)floor(x[0][1]);
    for(ii=0;ii<N0;ii++){
      a_p = Periodic_int(a + Dv0[ii][0], Dim_a);
      b_p = Periodic_int(b + Dv0[ii][1], Dim_b);
      /* check if particles are present in the cell (tau_cell[a_p][b_p]=sim_time_T) */
      if(tau_cell[a_p][b_p]==sim_time_T){
	i=phi_cell[a_p][b_p];
	do{
	  for(j=0;j<D;j++) v[j] = Shortest(x[0][j]-x[i][j], half_box[j], box[j]);
	  r2 = Mod2(v);
	  if(r2<cut0SQ){
	    e = fp0 * exp(-r2 * inv0);
	    for(j=0;j<D;j++){
	      force = e * v[j];
	      F[0][j] += force;
	      F[i][j] -= force;
	    }
	  }
	  i=nu_list[i];
	}while(i>=0);
      }
    }
  }

  
  /* polymer bonds */
  if(Lp>1){ // assuming k_pol>0 if Lp>1
    for(a=0;a<Mp;a++){
      i=i_ini+a*Lp;
      ii=i+1;
      for(b=1;b<Lp;b++){
	/* distance from the periodic image of x[ii] */
	for(j=0;j<D;j++) vb[i][j] = v[j] = Shortest(x[ii][j] - x[i][j], half_box[j], box[j]);
	r = Mod(v);
	// 2023.12.07 CHANGED e = k_pol * (sigma-r)/r; // 1/r to normalize v to a versor
	for(j=0;j<D;j++){
	  // 2023.12.07 CHANGED force = e * v[j];
	  force = -k_pol*(r-Lambda)*v[j]/r; // k_pol * v[j];
	  F[ii][j] += force;
	  F[ii][j] += f_active*v[j]/r;
	  F[i][j]  -= force;
	}
	i=ii;
	ii++;
      }
    }
  }

  /* small-small interactions */
  if(eps!=0){
    for(i=i_ini;i<=i_fin;i++){
      a = ix[i][0]; b = ix[i][1];

      for(Da=-1;Da<=1;Da++){
	a_p = Periodic_int(a+Da, Dim_a);
	for(Db=-1;Db<=1;Db++){
	  b_p = Periodic_int(b+Db, Dim_b);
	  if(tau_cell[a_p][b_p]==sim_time_T){
	    /* particles are present in the cell */
	    ii=phi_cell[a_p][b_p];
	    do{
	      /* ii>i counts pair ii-i only once */
	      if(ii>i){
		/* distance from the periodic image of x[ii] */
		if(ii==i+1 && BONDED[i]) for(j=0;j<D;j++) v[j] = vb[i][j]; // computed above
		else for(j=0;j<D;j++) v[j] = Shortest(x[ii][j] - x[i][j], half_box[j], box[j]);
		/* Gaussian repulsion */
		r2 = Mod2(v);
		if(r2<cutSQ){
		  e = fp * exp(-r2 * inv);
		  for(j=0;j<D;j++){
		    force = e * v[j];
		    F[ii][j] += force;
		    F[i][j]  -= force;
		  }
		}	
#ifdef DEB
		if(!THERMALIZING){
		  fprintf(stderr,
			  "\ni=%d (%.3lf,%.3lf) ii=%d (%.3lf,%.3lf) r=%.3lf  F=%.6lf*(%.4lf,%.4lf)",
			  i,x[i][0],x[i][1],ii,x[ii][0],x[ii][1],r,e,v[0],v[1]);
		  fprintf(stderr,"\n%lf (%lf)  exp=%lf f*exp=%lf  with fp=%lf",
			  inv,-pow(r,2) * inv,exp(-pow(r,2) * inv),e,fp);
		}
#endif
	      }
	      ii=nu_list[ii];
	    }while(ii>=0);
	  }
	}
      }
    }
  }

  /* Euler step */
  for(j=0; j<D; j++){
    x[0][j] += mu0_dt * F[0][j] + dB[0][j];
    for(i=i_ini;i<=i_fin;i++) x[i][j] += mu_dt * F[i][j] + dB[i][j];
  }
  /* taking care of periodic boundary conditions */
  i=0;
  for(j=0;j<D;j++){
    /* trap moved according to the probe's posit., even if the trap goes outside of the box */
    while(x[i][j]<0){      x[i][j] += box[j]; x_trap[j] += box[j]; travel[0][j]--; }
    while(x[i][j]>=box[j]){x[i][j] -= box[j]; x_trap[j] -= box[j]; travel[0][j]++; }
  }
  for(i=i_ini;i<=i_fin;i++){
    for(j=0;j<D;j++){
      while(x[i][j]<0){      x[i][j] += box[j]; travel[i][j]--;}
      while(x[i][j]>=box[j]){x[i][j] -= box[j]; travel[i][j]++;}
    }
  }
}


/* Q=1 is std, Q<1 for thermalization */
void Params(double Q){
  sigma=2*R; // reading R in input -> interaction range = 2R
  sigma0=R0+R;
  cut0= 4.*sigma0;
  fp0 = Q*eps0/(pow(sigma0,2));
  inv0= 1./(2*pow(sigma0,2));
  
  cut = 4.*sigma;
  fp  = Q*eps/(pow(sigma,2));
  inv = 1./(2*pow(sigma,2));
  cut0SQ= cut0*cut0;
  cutSQ= cut*cut;
  
  mu0    = 1./gam0;
  mu     = 1./gam;
  mu0_dt = mu0*dt;
  mu_dt  = mu*dt;
  stoch  = sqrt(2. * mu * T * dt);
  stoch0 = sqrt(2. * mu0 * T * dt);

  /* here not used, they were needed if walls were used */
  x_min = 0.;
  x_max = box[0];
  y_min = 4*R;
  y_max = box[1]-4*R;
  y0_min = 4*R0;
  y0_max = box[1]-4*R0;
}

void Init(){
  int m,n,i,j,a,b;
  double v[D];
  char comm[1000];
  /* RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR */
  gettimeofday(&tp,&zp);
  iseed=(tp.tv_sec+tp.tv_usec)%10000;  
  for(i=0;i<6;i++) vseed[i] = iseed+i*73;
  RngStream_SetPackageSeed(vseed);
  rngs = RngStream_CreateStream ("NMSM");
  for(i=0;i<1000;i++) RngStream_RandU01(rngs);
  /* RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR */

  i_ini=1; i_fin=N;
  x_trap = (double *)dvector(0,D-1);
  x = (double **)dmatrix(0,N,0,D-1);
  F = (double **)dmatrix(0,N,0,D-1);
  dB= (double **)dmatrix(0,N,0,D-1); 
  travel= (int **)imatrix(0,N,0,D-1);
  ix= (int **)imatrix(0,N,0,D-1);
  Dim_a = (int)floor(box[0]+1e-5);
  Dim_b = (int)floor(box[1]+1e-5);
  tau_cell = (int **)imatrix(0,Dim_a,0,Dim_b);
  phi_cell = (int **)imatrix(0,Dim_a,0,Dim_b);

  x_rel = (double *)dvector(0,MAX);
  
  vb = (double **)dmatrix(0,N,0,D-1);
  BONDED = (int *)ivector(0,N+1);
  for(i=0;i<=N;i++) BONDED[i]=0;
  for(a=0;a<Mp;a++){
    i=i_ini+a*Lp;
    for(b=1;b<Lp;b++){
      BONDED[i]=1;
      i++;
    }
  }
  for(a=0;a<50;a++) printf("\nBONDED[%2d] = %d",a,BONDED[a]);

  Params(1.);
  for(j=0;j<D;j++) half_box[j] = box[j]/2;
  printf("\nbox[0]=%.1lf, box[1]=%.1lf",box[0],box[1]);

  /*                 */
  /* place particles */
  /*                 */
  /* placing big particle at the center of the box */
  x[0][0] = box[0]/2.;
  x[0][1] = box[1]/2.;
  x_trap[0] = x[0][0];// + 1.5* v_trap*gam0/k_trap;
  x_trap[1] = box[1]/2.;

  /* list of cells close to the center of the big particle */
  m = (int)(cut0+1.);
  N0=0;
  for(a=-m;a<=m;a++){
    for(b=-m;b<=m;b++){
      v[0]=a; if(v[0]<0) v[0]++; if(v[0]>0) v[0]--; 
      v[1]=b; if(v[1]<0) v[1]++; if(v[1]>0) v[1]--; 
      if(Mod(v)<=cut0){
	Dv0[N0][0]=a;
	Dv0[N0][1]=b;
	N0++;
      }
    }
  }
  printf("\n[%d,%d], [%d,%d], .... [%d,%d]:  N0=%d",
	  Dv0[0][0],Dv0[0][1],Dv0[1][0],Dv0[1][1],Dv0[N0-1][0],Dv0[N0-1][1],N0);
  
  /* randomly place small particles;
     if Lp>1 they are joined in horizontal polymers of length Lp */
  i=i_ini;
  for(m=0;m<Mp;m++){
    for(n=0;n<Lp;n++){
      if(n==0)
	for(j=0;j<D;j++) x[i][j] = box[j] * RngStream_RandU01(rngs);
      else{
	for(j=0;j<D;j++) x[i][j] = x[i-1][j];
	x[i][0] += 2*R;
	while(x[i][0]>box[0]){
	  x[i][0] -= box[0];
	  travel[i][0]++;
	}
      }
      i++;
    }
  }
  
  Make_str0("");
  
  Print_info(stderr);
  printf("\n# -------------------------------------------------------");
  printf("\n# %s",str0);
  printf("\n# -------------------------------------------------------\n");
  sprintf(comm,"mkdir %s; mkdir %s; mkdir OUT; mkdir FIG",dirDATA,dirCONF);
  system(comm);
  sprintf(comm,"mkdir DATA/%s",str0);
  system(comm);
}

/* start with smaller epsilons and gradually increase them */
void Thermalize(int n, int z){
	int t;
	double Q;
	Q = pow(2, z);
	printf("\nx[0]=(%lf,%lf)", x[0][0], x[0][1]);
	while (z <= 0)
	{
		Params(Q);
		printf("\ntherm %2d (%lf)... ", z, Q);
		for (t = 0; t < n; t++)
			MD_step();
		printf(" ok  x[0]=(%lf,%lf)", x[0][0], x[0][1]);
		z++;
		Q = pow(2, z);
		Print_elapsed_time();
	}
	Params(1.);
	THERMALIZING = 0;
}


void Save_time_series(double v, double *x, int imin, int imax){
  FILE *out1;
  char fn[1000];
  int i;
  sprintf(fn,"%s/%s/x-rel_v_%lf.dat",dirDATA,str0,v);
  if(imin==0) out1=fopen(fn,"w");
  else        out1=fopen(fn,"a");
  if(out1==NULL){ fprintf(stderr,"\n OUTPUT FILE ERROR (%s)\n",fn);return;}
  for(i=imin;i<=imax;i++) fprintf(out1,"%.4lf\n",x[i]);
  fclose(out1);
}

void Save_traj(int nrun, double ti, double v){
  FILE *out1;
  int i;
  sprintf(fnO,"%s/traj__%s_v%lf_t%.2lf.dat",dirCONF,str0,v,nstep_save*Dt);
  if(FIRST_SAVE){ out1=fopen(fnO,"w"); FIRST_SAVE=0;}
  else out1=fopen(fnO,"a");  
  if(out1==NULL){ fprintf(stderr,"\n OUTPUT FILE ERROR\n");return;}
  fprintf(out1,"%d %lf 0 0 \n",nrun,ti);
  for(i=0;i<=i_fin;i++) fprintf(out1,"%.4lf %.4lf %d %d\n",
				x[i][0],x[i][1],travel[i][0],travel[i][1]);
  fprintf(out1,"%.6lf %.1lf %d %d\n",x_trap[0],x_trap[1],travel[0][0],travel[0][1]);
  fclose(out1);
}

int int_coor(double x, double SF){ return (int)(SF*x); }

void Save_conf(int nrun){
  FILE *out1;
  int i;
  double SF=1000.;
  sprintf(fnO,"%s/conf__%s_%06d.dat",dirCONF,str0, nrun);
  out1=fopen(fnO,"w"); 
  if(out1==NULL){ fprintf(stderr,"\n OUTPUT FILE ERROR (%s)\n",fnO);return;}
  for(i=0;i<=i_fin;i++) fprintf(out1,"%d %d\n",int_coor(x[i][0],SF),int_coor(x[i][1],SF));
  fprintf(out1,"%d %d\n",int_coor(x_trap[0],SF),int_coor(x_trap[1],SF));
  fclose(out1);
}


/* ----------------------------- Il programma */
int main(int argn, char **argc) 
{
  int s, iv, it, it0, N_therm=100000,z_therm=-2;
  time(&time_start);
  
  printf("\nsyntax:",tt); 
  printf("\nG-1-pp  nstep_save  N L  Bx By  T  v_ini Nv v_x_dec  k kp  R eps  R0 eps0 lambda f_active [dt=%.4lf] [tt=%.0lf]",dt,tt); 
  printf("\n        ++++++++++  **********  -  ^^^^^^^^^^^^^^^^  ****  -----  -------\n"); 
  printf("\ndt=%0.5lf (integration time step)",dt); 
  printf("\nDt=%0.5lf (samplinx-x0 time step)",Dt); 
  printf("\nnstep_save=0 ==> no saving");
  printf("\nnstep_save>0 ==> saving configs every nstep_save steps\n"); 
  
  if(argn<16){ printf("\n%d \n",argn);exit(0);}
  
  s=0;
  s++; nstep_save=atoi(argc[s]); // saving configs if >0
  s++; N=atoi(argc[s]);          // nr of small particles
  i_fin=N; 
  s++; Lp=atoi(argc[s]);         // length of each polymer
  box=(double *)dvector(0,D-1); 
  s++; box[0]=atof(argc[s]);     // box size, x
  s++; box[1]=atof(argc[s]);     // box size, y
  s++; T=atof(argc[s]);          // temperature
  s++; v_trap_ini=atof(argc[s]); // smallest velocity of the trap
  s++; Nv=atoi(argc[s]);         // number of velocities
  s++; v_per_decade=atof(argc[s]); // velocities per decade
  EF = pow(10.,1./v_per_decade);
  s++; k_trap=atof(argc[s]);     // stiffness of the trap 
  s++; k_pol=atof(argc[s]);      // stiffness of the polymer bonds
  s++; R=atof(argc[s]);          // "radius" of each particle
  if(R>0.25){
    fprintf(stderr,"\ninteraction range 8R=%lf > cell side =1: STOP\n",8*R);
    exit(33);
  }
  s++; eps=atof(argc[s]);        // repulsive energy of particles
  s++; R0=atof(argc[s]);         // "radius" of the probe
  s++; eps0=atof(argc[s]);       // repulsive energy of the probe
  
  s++; Lambda=atof(argc[s]);         // lambda constant
  s++; f_active=atof(argc[s]);       // f active dumbells
  
  s++; if(argn>s) dt=atof(argc[s]); // integration time step
  s++; if(argn>s) tt=atof(argc[s]); // total time of the simulation
  
  
  printf("\nSimulating at various v's and saving time series\n\n");
  strcpy(dirDATA,"DATA");
  strcpy(dirCONF,"CONF");
  if(nstep_save) printf("\nSaving configurations every %d steps",nstep_save);

  fprintf(stderr, "N = %d, L = %d, Lambda = %lf, f_active = %lf\nbox_x = %lf, box_y = %lf\ndt=%lf\n", N, Lp, Lambda, f_active, box[0], box[1], dt);

  Mp=(int)(N/Lp);
  if(Mp*Lp != N){
    fprintf(stderr,"\nN=%d is not divisible by Lp=%d: stopping\n",N,Lp);
    exit(55);
  }
  else printf("\nMp = N / Lp = %d",Mp);

  
  v_trap=v_trap_ini;
  Init();
  nstep_x0=(int)(Dt/dt);
  printf("dt=%lf  Dt=%lf  nstep_x0 = %d",dt,Dt,nstep_x0);

  /* NO video = long data acquisition, x-x_trap ----------------------- */
  for(iv=1;iv<=Nv;iv++){
    /* relax to a steady state */
    printf("\n%3d/%d) relaxing %s at v = %lf",iv,Nv,str0,v_trap);
    THERMALIZING=1;
    Thermalize(N_therm,z_therm);
    z_therm=0; // using z_therm<0 only in the first thermalization, to separate overlapped particles
    N_therm=50000;
    
    /* sampling */
    printf("\n starting sampling %s\n(tt=%.0lf)",str0,tt);
    Print_elapsed_time();
    nrun=1;
    it=it0=0;
    ta=0.; // total accumulated time
    while(ta<=tt-Dt+1e-10){
      for(s=0;s<nstep_x0;s++) MD_step();
      x_rel[it] = x[0][0] - x_trap[0];
      ta += Dt;
      it++;
      if(it%10000==0){
	printf("\n%.0lf %.2lf,%.2lf ",ta,x[0][0],x[0][1]);
	Print_elapsed_time();
	Save_time_series(v_trap, x_rel, it0, it-1);
	it0=it;
      }
      if(nstep_save){
	/* saving configuration */
	if(it%nstep_save==0){
	  Save_traj(nrun,ta,v_trap);
	  nrun++;
	}
      }
    }
    if(it0<it) Save_time_series(v_trap, x_rel, it0, it-1);
    printf("\n%3d/%d) done %s at v = %lf",iv,Nv,str0,v_trap);
    v_trap *= EF;
    Print_elapsed_time();
  }

  return 0;
}
