#define min(a,b)    ({ __typeof__ (a) _a = (a);  __typeof__ (b) _b = (b); _a < _b ? _a : _b; })

double compute_energy_somemodel(int i, int j){

    return 0.;
}


/********************************************************************************************/

double compute_energy_translation(){

    return 0;
}


double compute_energy(double **X, double m){
	
	double en = 0.0;
    	// potential energy
    	double r = 0.0;  // distance between particle i and j    	
    	for(int i = 0; i< mySys.NPart; i++){
    	    	for(int j = 0; j< mySys.NPart; j++){
    			r = sqrt( pow(X[i][0]-X[j][0],2) + pow(X[i][1]-X[j][1],2) + pow(X[i][2]-X[j][2],2));
    			if (r <= 2*mySys.radius){
    				en += 1e7;
    			}
    			else if (r > 2*mySys.radius){
    				en += 0;
    			}
		}
	}
	
    return en;

}

double update_energy(double **X, double **X_new, double en_old, int i_star){
	
	double en_new = en_old; // was computed on X
	double en_diff  = 0;
    	// potential energy
    	double r = 0.0;  // distance between particle i and j    	
    	for(int i = 0; i< mySys.NPart; i++){
		r = sqrt( pow(X[i][0]-X[i_star][0],2) + pow(X[i][1]-X[i_star][1],2) + pow(X[i][2]-X[i_star][2],2));
		if (r <= 2*mySys.radius){
			en_new -= 1e7;
			en_diff -= 1e7;
		}
		else if (r > 2*mySys.radius){
			en_new -= 0;
			en_diff -=0;
		}
	}
	
	for(int i = 0; i< mySys.NPart; i++){
		r = sqrt( pow(X_new[i][0]-X_new[i_star][0],2) + pow(X_new[i][1]-X_new[i_star][1],2) + pow(X_new[i][2]-X_new[i_star][2],2));
		if (r <= 2*mySys.radius){
			en_new += 1e7;
			en_diff += 1e7;
		}
		else if (r > 2*mySys.radius){
			en_new += 0;
			en_diff += 0;
		}
	}
	
	
    return en_diff;

}

