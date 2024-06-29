#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


void Copy(double *u, double *w){
  w[0] = u[0];
  w[1] = u[1];
}

double Dot(double *u, double *v){
  return u[0]*v[0]+u[1]*v[1];
}

void Prod(double *u, double A, double *v){
  v[0] = A*u[0];
  v[1] = A*u[1];
}

void Diff(double *u, double *v, double *w){
  w[0] = u[0]-v[0];
  w[1] = u[1]-v[1];
}

void Sum(double *u, double *v, double *w){
  w[0] = u[0]+v[0];
  w[1] = u[1]+v[1];
}

void Incr(double *u, double *du){
  u[0] += du[0];
  u[1] += du[1];
}

void Decr(double *u, double *du){
  u[0] -= du[0];
  u[1] -= du[1];
}


void Mean(double *u, double *v, double *w){
  w[0] = 0.5*(u[0]+v[0]);
  w[1] = 0.5*(u[1]+v[1]);
}

double Mod2(double *u){
  return u[0]*u[0]+u[1]*u[1];
}

double Mod(double *u){
  return sqrt(u[0]*u[0]+u[1]*u[1]);
}

double Dist(double *u, double *v){
  double w[3];
  w[0] = u[0]-v[0];
  w[1] = u[1]-v[1];
  return sqrt(w[0]*w[0]+w[1]*w[1]);
}

void SetLen(double *u, double L){
  double f; int j;
  f=L/Mod(u);
  for(j=0;j<2;j++) u[j] *= f;
}
