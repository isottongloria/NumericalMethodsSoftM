#define min(a,b)    ({ __typeof__ (a) _a = (a);  __typeof__ (b) _b = (b); _a < _b ? _a : _b; })

double compute_energy_somemodel(int i, int j){

    return 0.;
}


/********************************************************************************************/

double compute_energy_translation(){

    return 0;
}


double compute_energy(double **V, double m){
	double en = 0.0;
    	double r = 0.0;  // distance between particle i and j
    	for(int i = 0; i< mySys.NPart; i++){
		//r = sqrt( pow(X[i][0]-X[j][0],2) + pow(X[i][1]-X[j][1],2) + pow(X[i][2]-X[j][2],2));
		// en += 4 * eps * (pow((sigma / r),12) - pow((sigma / r),6));
		en += 0.5 * m * (pow(V[i][0],2) + pow(V[i][1],2) + pow(V[i][2],2));      	
    	} 
   
    return en;

}


