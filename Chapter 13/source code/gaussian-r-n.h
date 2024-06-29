void Gaussian_(double *y1, double *y2);
void BoxMuller(double *y1, double *y2);
void Fresh_Gaussian();
void Wallace_Gaussian();
void Shuffle(int *v, int N);
void Set_Wallace_method(int nslot, int nrow, int nturn);

extern int iG,iGnew;
extern double vGauss[2][10000000];
//extern double **vGauss;
