/*************************************************************
*  This file contains routines that allocate by using        *
*  calloc(), matrices of different dimensions
*************************************************************/

/* array_alloc.c */

/* Note: matrix allocation routine assumes whole matrix can be 
   allocated in one call to calloc.  If not true (say for large
   matrices on an IBM PC), dmatrix() and free_matrix() will
   have to be modified */

#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include "array_alloc.h"

/* error routine */
void error(msg)
char *msg;
{
puts(msg);
exit(1);
return;
}

void NRError(error_text)
char *error_text;
{ 
  error(error_text);
}

/*********************************************************************
*
*  function: mycalloc
*
*  purpose: memory allocation with error detection
*/

char *mycalloc(num,size)
int num,size;
{ char *ptr;
  
  if ( (num < 0) || (size <0) )
    error("Trying to allocate negative size memory block.\n");
  if ( num == 0 ) num = 1; /* to prevent problem with zero block size */
  if ( size == 0 ) size = 1;
  ptr = (char *)calloc(num,size);
  if ( ptr == NULL )
    error("Cannot allocate memory.\n");
  return ptr;
}  


/*********************************************************************
*
*  function: kb_realloc
*
*  purpose: memory re-allocation with error detection
*           Needs oldsize as input so can initialize new memory to 0
*/

char *kb_realloc(ptr,size,oldsize)
char *ptr;
int size;
int oldsize;
{
  ptr = (char *)realloc(ptr,size);
  if ( ptr == NULL )
    error("Cannot reallocate memory.\n");
  else if (oldsize < size )  memset(ptr+oldsize,0,size-oldsize);
  return ptr;
}  



/******************************************************************
 *                  character allocation                            *
 ******************************************************************/

char *cvector(lo,hi)
int lo,hi;
/* allocates a character vector with range [lo..hi] */
{
  char *v;

  v = (char *)mycalloc((unsigned)(hi-lo+1),sizeof(char));
  if ( !v ) NRError("allocation failure in cvector().");
  return v-lo;
}

char **cmatrix(rlo,rhi,clo,chi)
int rlo,rhi,clo,chi;
/* allocates a character matrix with range [rlo..rhi][clo..chi] */
{
  int i;
  char **m;

  if ( rhi-rlo+1 == 0 ) return NULL;
  m = (char **)mycalloc((unsigned)(rhi-rlo+1),sizeof(char *));
  if ( !m ) NRError("cmatrix allocation error.");
  m -= rlo;

  m[rlo] = (char *) mycalloc((unsigned)(chi-clo+1),(rhi-rlo+1)*sizeof(char));
  if ( !m[rlo] ) NRError("cmatrix allocation error.");
  for ( i = rlo+1 ; i <= rhi ; i++ )
      m[i] = m[rlo] + (i - rlo)*(chi - clo + 1) - clo;

  return m;
}

char ***cmatrix3(n1,n2,n3)
int n1,n2,n3;
/* allocates a char matrix with range [0..n1-1][0..n2-1][0..n3-1] */
{
  int i,j;
  char ***m;

  /* assumes all pointers same machine size and alignment */
  m = (char ***)mycalloc((n2+1)*n1,sizeof(char *));
  if ( !m ) NRError("cmatrix3 allocation error.");

  for ( i = 0 ; i < n1 ; i++ )
    m[i] = (char **)(m + n1 + i*n2);
  m[0][0] = (char *) mycalloc(n1*n2*n3,sizeof(char));
  if ( !m[0][0] ) NRError("cmatrix3 allocation error.");
  for ( i = 0 ; i < n1 ; i++ )
    for ( j = 0 ; j < n2 ; j++ )
      m[i][j] = m[0][0] + i*n2*n3 + j*n3;

  return m;
}


/******************************************************************
 *                  integer allocation                            *
 ******************************************************************/

int *ivector(lo,hi)
int lo,hi;
/* allocates a int vector with range [lo..hi] */
{
  int *v;

  v = (int *)mycalloc((unsigned)(hi-lo+1),sizeof(int));
  if ( !v ) NRError("allocation failure in ivector().");
  return v-lo;
}

unsigned int *uivector(lo,hi)
int lo,hi;
/* allocates a int vector with range [lo..hi] */
{
  unsigned int *v;

  v = (unsigned int *)mycalloc((unsigned)(hi-lo+1),sizeof(int));
  if ( !v ) NRError("allocation failure in ivector().");
  return v-lo;
}

long int *livector(lo,hi)
int lo,hi;
/* allocates a int vector with range [lo..hi] */
{
  long int *v;

  v = (long int *)mycalloc((unsigned)(hi-lo+1),sizeof(long int));
  if ( !v ) NRError("allocation failure in ivector().");
  return v-(long int)lo;
}


long long int *llivector(lo,hi)
int lo,hi;
/* allocates a int vector with range [lo..hi] */
{
  long long int *v;

  v = (long long int *)mycalloc((unsigned)(hi-lo+1),sizeof(long long int));
  if ( !v ) NRError("allocation failure in llivector().");
  return v-(long long int)lo;
}


int **imatrix(rlo,rhi,clo,chi)
int rlo,rhi,clo,chi;
/* allocates a int matrix with range [rlo..rhi][clo..chi] */
{
  int i;
  int **m;

  if ( rhi-rlo+1 == 0 ) return NULL;
  m = (int **)mycalloc((unsigned)(rhi-rlo+1),sizeof(int *));
  if ( !m ) NRError("imatrix allocation error.");
  m -= rlo;

  m[rlo] = (int *) mycalloc((unsigned)(chi-clo+1),(rhi-rlo+1)*sizeof(int));
  if ( !m[rlo] ) NRError("imatrix allocation error.");
  for ( i = rlo+1 ; i <= rhi ; i++ )
      m[i] = m[rlo] + (i - rlo)*(chi - clo + 1) - clo;

  return m;
}

int ***imatrix3(n1,n2,n3)
int n1,n2,n3;
/* allocates a int matrix with range [0..n1-1][0..n2-1][0..n3-1] */
{
  int i,j;
  int ***m;

  /* assumes all pointers same machine size and alignment */
  m = (int ***)mycalloc((n2+1)*n1,sizeof(int *));
  if ( !m ) NRError("imatrix3 allocation error.");

  for ( i = 0 ; i < n1 ; i++ )
    m[i] = (int **)(m + n1 + i*n2);
  m[0][0] = (int *) mycalloc(n1*n2*n3,sizeof(int));
  if ( !m[0][0] ) NRError("imatrix3 allocation error.");
  for ( i = 0 ; i < n1 ; i++ )
    for ( j = 0 ; j < n2 ; j++ )
      m[i][j] = m[0][0] + i*n2*n3 + j*n3;

  return m;
}

int ****imatrix4(n1,n2,n3,n4)
int n1,n2,n3,n4;
/* allocates a int matrix with range [0..n1-1][0..n2-1][0..n3-1][0..n4-1] */
{
  int i,j,k;
  int ****m;

  /* assumes all pointers same machine size and alignment */
  m = (int ****)mycalloc(n1+n1*n2+n1*n2*n3,sizeof(int *));
  if ( !m ) NRError("imatrix4 allocation error.");

  for ( i = 0 ; i < n1 ; i++ )
    { m[i] = (int ***)(m + n1 + i*n2);
      for ( j = 0 ; j < n2 ; j++ )
        m[i][j] = (int **)(m + n1 + n1*n2 + i*n2*n3 + j*n3);
    }
  m[0][0][0] = (int *) mycalloc(n1*n2*n3*n4,sizeof(int));
  if ( !m[0][0][0] ) NRError("imatrix4 allocation error.");
  for ( i = 0 ; i < n1 ; i++ )
    for ( j = 0 ; j < n2 ; j++ )
      for ( k = 0 ; k < n3 ; k++ )
        m[i][j][k] = m[0][0][0] + i*n2*n3*n4 + j*n3*n4 + k*n4;

  return m;
}

/* routines for initializing matrices declared as local variables */
int imat2d_setup(name,spacename,rows,cols)
int **name;
int *spacename;
int rows,cols;
{ for ( ; rows > 0 ; rows--,spacename += cols,name++ )
    *name = spacename;
  return 0;
}

int imat3d_setup(name,spacename,rows,cols,levels)
int ***name;
int *spacename;
int rows,cols,levels;
{ int i;
  int **row = (int **)(name + rows);
  for ( ; rows > 0 ; rows--,name++ )
  { *name = row;
    for ( i = 0 ; i < cols ; i++,spacename += levels, row++ )
      *row = spacename;
  }
  return 0;
}

int imat4d_setup(name,spacename,rows,cols,levels,tiers)
int ****name;
int *spacename;
int rows,cols,levels,tiers;
{ int i,j;
  int ***row = (int ***)(name + rows);
  int **col = (int **)(name + rows + rows*cols);
  for ( ; rows > 0 ; rows--,name++ )
  { *name = row;
    for ( i = 0 ; i < cols ; i++, row++ )
    { *row = col;
      for ( j = 0 ; j < levels ; j++,spacename += tiers, col++ )
        *col = spacename;
    }
  }
  return 0;
}
/* end local declaration routines */

void free_ivector(v,lo,hi)
int *v,lo,hi;
{
  free((char *)(v+lo));
}

void free_imatrix(m)
int **m;
{
  if ( !m ) return;
  free((char *)m[0]);
  free((char *)m);
}

void free_imatrix3(m)
int ***m;
{
  if ( !m ) return;
  free((char *)m[0][0]);
  free((char *)m);
}

void free_imatrix4(m)
int ****m;
{
  if ( !m ) return;
  free((char *)m[0][0][0]);
  free((char *)m);
}


/******************************************************************
 *                  floating-point allocation                     *
 ******************************************************************/
float *fvector(lo,hi)
int lo,hi;
/* allocates a float vector with range [lo..hi] */
{
  float *v;

  v = (float *)mycalloc((unsigned)(hi-lo+1),sizeof(float));
  if ( !v ) NRError("allocation failure in fvector().");
  return v-lo;
}

float **fmatrix(rlo,rhi,clo,chi)
int rlo,rhi,clo,chi;
/* allocates a float matrix with range [rlo..rhi][clo..chi] */
{
  int i;
  float **m;

  if ( rhi-rlo+1 == 0 ) return NULL;
  m = (float **)mycalloc((unsigned)(rhi-rlo+1),sizeof(float *));
  if ( !m ) NRError("fmatrix allocation error.");
  m -= rlo;

  m[rlo] = (float *) mycalloc((unsigned)(chi-clo+1),(rhi-rlo+1)*sizeof(float));
  if ( !m[rlo] ) NRError("fmatrix allocation error.");
  for ( i = rlo+1 ; i <= rhi ; i++ )
      m[i] = m[rlo] + (i - rlo)*(chi - clo + 1) - clo;

  return m;
}

float ***fmatrix3(n1,n2,n3)
int n1,n2,n3;
/* allocates a float matrix with range [0..n1-1][0..n2-1][0..n3-1] */
{
  int i,j;
  float ***m;

  /* assumes all pointers same machine size and alignment */
  m = (float ***)mycalloc((n2+1)*n1,sizeof(float *));
  if ( !m ) NRError("fmatrix3 allocation error.");

  for ( i = 0 ; i < n1 ; i++ )
    m[i] = (float **)(m + n1 + i*n2);
  m[0][0] = (float *) mycalloc(n1*n2*n3,sizeof(float));
  if ( !m[0][0] ) NRError("fmatrix3 allocation error.");
  for ( i = 0 ; i < n1 ; i++ )
    for ( j = 0 ; j < n2 ; j++ )
      m[i][j] = m[0][0] + i*n2*n3 + j*n3;

  return m;
}

float ****fmatrix4(n1,n2,n3,n4)
int n1,n2,n3,n4;
/* allocates a float matrix with range [0..n1-1][0..n2-1][0..n3-1][0..n4-1] */
{
  int i,j,k;
  float ****m;

  /* assumes all pointers same machine size and alignment */
  m = (float ****)mycalloc(n1+n1*n2+n1*n2*n3,sizeof(float *));
  if ( !m ) NRError("fmatrix4 allocation error.");

  for ( i = 0 ; i < n1 ; i++ )
    { m[i] = (float ***)(m + n1 + i*n2);
      for ( j = 0 ; j < n2 ; j++ )
        m[i][j] = (float **)(m + n1 + n1*n2 + i*n2*n3 + j*n3);
    }
  m[0][0][0] = (float *) mycalloc(n1*n2*n3*n4,sizeof(float));
  if ( !m[0][0][0] ) NRError("fmatrix4 allocation error.");
  for ( i = 0 ; i < n1 ; i++ )
    for ( j = 0 ; j < n2 ; j++ )
      for ( k = 0 ; k < n3 ; k++ )
        m[i][j][k] = m[0][0][0] + i*n2*n3*n4 + j*n3*n4 + k*n4;

  return m;
}

/* routines for initializing matrices declared as local variables */
int fmat2d_setup(name,spacename,rows,cols)
float **name;
float *spacename;
int rows,cols;
{ for ( ; rows > 0 ; rows--,spacename += cols,name++ )
    *name = spacename;
  return 0;
}

int fmat3d_setup(name,spacename,rows,cols,levels)
float ***name;
float *spacename;
int rows,cols,levels;
{ int i;
  float **row = (float **)(name + rows);
  for ( ; rows > 0 ; rows--,name++ )
  { *name = row;
    for ( i = 0 ; i < cols ; i++,spacename += levels, row++ )
      *row = spacename;
  }
  return 0;
}

int fmat4d_setup(name,spacename,rows,cols,levels,tiers)
float ****name;
float *spacename;
int rows,cols,levels,tiers;
{ int i,j;
  float ***row = (float ***)(name + rows);
  float **col = (float **)(name + rows + rows*cols);
  for ( ; rows > 0 ; rows--,name++ )
  { *name = row;
    for ( i = 0 ; i < cols ; i++, row++ )
    { *row = col;
      for ( j = 0 ; j < levels ; j++,spacename += tiers, col++ )
        *col = spacename;
    }
  }
  return 0;
}
/* end local declaration routines */


void free_fvector(v,lo,hi)
float *v;
int lo,hi;
{
  free((char *)(v+lo));
}

void free_fmatrix(m)
float **m;
{
  if ( !m ) return;
  free((char *)m[0]);
  free((char *)m);
}

void free_fmatrix3(m)
float ***m;
{
  if ( !m ) return;
  free((char *)m[0][0]);
  free((char *)m);
}

void free_fmatrix4(m)
float ****m;
{
  if ( !m ) return;
  free((char *)m[0][0][0]);
  free((char *)m);
}



/******************************************************************
 *                  double precision allocation                   *
 ******************************************************************/
double *dvector(lo,hi)
int lo,hi;
/* allocates a double vector with range [lo..hi] */
{
  double *v;

  v = (double *)mycalloc((unsigned)(hi-lo+1),sizeof(double));
  if ( !v ) NRError("allocation failure in dvector().");
  return v-lo;
}

double **dmatrix(rlo,rhi,clo,chi)
int rlo,rhi,clo,chi;
/* allocates a double matrix with range [rlo..rhi][clo..chi] */
{
  int i;
  double **m;

  if ( rhi-rlo+1 == 0 ) return NULL;
  m = (double **)mycalloc((unsigned)(rhi-rlo+1),sizeof(double *));
  if ( !m ) NRError("dmatrix allocation error.");
  m -= rlo;

  m[rlo] = (double *) mycalloc((unsigned)(chi-clo+1),(rhi-rlo+1)*sizeof(double));
  if ( !m[rlo] ) NRError("dmatrix allocation error.");
  for ( i = rlo+1 ; i <= rhi ; i++ )
      m[i] = m[rlo] + (i - rlo)*(chi - clo + 1) - clo;

  return m;
}

double ***dmatrix3(n1,n2,n3)
int n1,n2,n3;
/* allocates a double matrix with range [0..n1-1][0..n2-1][0..n3-1] */
{
  int i,j;
  double ***m;

  /* assumes all pointers same machine size and alignment */
  m = (double ***)mycalloc((n2+1)*n1,sizeof(double *));
  if ( !m ) NRError("dmatrix3 allocation error.");

  for ( i = 0 ; i < n1 ; i++ )
    m[i] = (double **)(m + n1 + i*n2);
  m[0][0] = (double *) mycalloc(n1*n2*n3,sizeof(double));
  if ( !m[0][0] ) NRError("dmatrix3 allocation error.");
  for ( i = 0 ; i < n1 ; i++ )
    for ( j = 0 ; j < n2 ; j++ )
      m[i][j] = m[0][0] + i*n2*n3 + j*n3;

  return m;
}

double ****dmatrix4(n1,n2,n3,n4)
int n1,n2,n3,n4;
/* allocates a double matrix with range [0..n1-1][0..n2-1][0..n3-1][0..n4-1] */
{
  int i,j,k;
  double ****m;

  /* assumes all pointers same machine size and alignment */
  m = (double ****)mycalloc(n1+n1*n2+n1*n2*n3,sizeof(double *));
  if ( !m ) NRError("dmatrix4 allocation error.");

  for ( i = 0 ; i < n1 ; i++ )
    { m[i] = (double ***)(m + n1 + i*n2);
      for ( j = 0 ; j < n2 ; j++ )
        m[i][j] = (double **)(m + n1 + n1*n2 + i*n2*n3 + j*n3);
    }
  m[0][0][0] = (double *) mycalloc(n1*n2*n3*n4,sizeof(double));
  if ( !m[0][0][0] ) NRError("dmatrix4 allocation error.");
  for ( i = 0 ; i < n1 ; i++ )
    for ( j = 0 ; j < n2 ; j++ )
      for ( k = 0 ; k < n3 ; k++ )
        m[i][j][k] = m[0][0][0] + i*n2*n3*n4 + j*n3*n4 + k*n4;

  return m;
}

/* routines for initializing matrices declared as local variables */
int dmat2d_setup(name,spacename,rows,cols)
double **name;
double *spacename;
int rows,cols;
{ for ( ; rows > 0 ; rows--,spacename += cols,name++ )
    *name = spacename;
  return 0;
}

int dmat3d_setup(name,spacename,rows,cols,levels)
double ***name;
double *spacename;
int rows,cols,levels;
{ int i;
  double **row = (double **)(name + rows);
  for ( ; rows > 0 ; rows--,name++ )
  { *name = row;
    for ( i = 0 ; i < cols ; i++,spacename += levels, row++ )
      *row = spacename;
  }
  return 0;
}

int dmat4d_setup(name,spacename,rows,cols,levels,tiers)
double ****name;
double *spacename;
int rows,cols,levels,tiers;
{ int i,j;
  double ***row = (double ***)(name + rows);
  double **col = (double **)(name + rows + rows*cols);
  for ( ; rows > 0 ; rows--,name++ )
  { *name = row;
    for ( i = 0 ; i < cols ; i++, row++ )
    { *row = col;
      for ( j = 0 ; j < levels ; j++,spacename += tiers, col++ )
        *col = spacename;
    }
  }
  return 0;
}
/* end local declaration routines */


void free_dvector(v,lo,hi)
double *v;
int lo,hi;
{
  free((char *)(v+lo));
}

void free_dmatrix(m)
double **m;
{
  if ( !m ) return;
  free((char *)m[0]);
  free((char *)m);
}

void free_dmatrix3(m)
double ***m;
{
  if ( !m ) return;
  free((char *)m[0][0]);
  free((char *)m);
}

void free_dmatrix4(m)
double ****m;
{
  if ( !m ) return;
  free((char *)m[0][0][0]);
  free((char *)m);
}

