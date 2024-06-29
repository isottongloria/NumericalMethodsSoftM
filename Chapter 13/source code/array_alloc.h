/* array_alloc.h */

void error(char *msg);
void nrerror(char *error_text);
char *mycalloc(int num,int size);
char *kb_realloc(char *ptr,int size,int oldsize );

/******************************************************************
 *                  character  allocation                         *
 ******************************************************************/
 char *cvector(int lo,int hi);
 char **cmatrix(int rlo,int rhi,int clo,int chi);
 char ***cmatrix3(int n1,int n2,int n3);
/******************************************************************
 *                  integer  allocation                           *
 ******************************************************************/
unsigned int *uivector(int lo, int hi);
unsigned int **uimatrix(int rlo,int rhi,int clo,int chi);
long long int *llivector(int lo,int hi);
int *ivector(int lo,int hi);
int **imatrix(int rlo,int rhi,int clo,int chi);
long long int **llimatrix(int rlo,int rhi,int clo,int chi);
int ***imatrix3(int n1,int n2,int n3);
int ****imatrix4(int n1,int n2,int n3,int n4);
int imat2d_setup(int **name,int *spacename,int rows,int cols);
int imat3d_setup(int ***name,int *spacename,int rows,int cols,int levels);
int imat4d_setup(int ****name,int *spacename,int rows,int cols,int levels,int tiers);
void free_ivector(int *v,int lo,int hi);
void free_imatrix(int **m);
void free_imatrix3(int ***m);
void free_imatrix4(int ****m);

/******************************************************************
 *                  floating-point allocation                     *
 ******************************************************************/
float *fvector(int lo,int hi);
float **fmatrix(int rlo,int rhi,int clo,int chi);
float ***fmatrix3(int n1,int n2,int n3);
float ****fmatrix4(int n1,int n2,int n3,int n4);
int fmat2d_setup(float **name,float *spacename,int rows,int cols);
int fmat3d_setup(float ***name,float *spacename,int rows,int cols,int levels);
int fmat4d_setup(float ****name,float *spacename,int rows,int cols,int levels,int tiers);
void free_fvector(float *v,int lo,int hi);
void free_fmatrix(float **m);
void free_fmatrix3(float ***m);
void free_fmatrix4(float ****m);
/******************************************************************
 *                  double precision allocation                   *
 ******************************************************************/
double *dvector(int lo,int hi);
double **dmatrix(int rlo,int rhi,int clo,int chi);
double ***dmatrix3(int n1,int n2,int n3);
double ****dmatrix4(int n1,int n2,int n3,int n4);
int dmat2d_setup(double **name,double *spacename,int rows,int cols);
int dmat3d_setup(double ***name,double *spacename,int rows,int cols,int levels);
int dmat4d_setup(double ****name,double *spacename,int rows,int cols,int levels,int tiers);
void free_dvector(double *v,int lo,int hi);
void free_dmatrix(double **m);
void free_dmatrix3(double ***m);
void free_dmatrix4(double ****m);

