
#define min(a,b)    ({ __typeof__ (a) _a = (a);  __typeof__ (b) _b = (b); _a < _b ? _a : _b; })

double compute_energy_somemodel(int i, int j){

    return 0.;
}


/********************************************************************************************/

double compute_energy_translation(){

    return 0;
}


double compute_energy(double* v){

    int i;
    double en = 0.;
    long int index;
    long int ix, iy, iz;
    long int j;
    double v2=0;

    for(i = 0; i< mySys.NPart; i++){
    	v2 = v[0,i] * v[0,i] + v[1,i] * v[1,i] + v[2,i] * v[2,i];
    	en += v2;   
    }
   
    return en;

}

