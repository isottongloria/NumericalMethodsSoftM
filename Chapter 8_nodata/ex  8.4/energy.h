#define min(a,b)    ({ __typeof__ (a) _a = (a);  __typeof__ (b) _b = (b); _a < _b ? _a : _b; })
#include <iostream>
using namespace std;



double compute_energy(double **X){
	

	double en = 0.0;
	double dx= 0., dy= 0., dz = 0.;
	
    	// lennard jones potential energy
    	double r = 0.0;  // distance between particle i and j    	
    	for(int i = 0; i< mySys.NPart; i++){
    	    	for(int j = 0; j< mySys.NPart; j++){
    	    	
    	    		if (i != j){
	    			//r = sqrt( pow(X[i][0]-X[j][0],2) + pow(X[i][1]-X[j][1],2) + pow(X[i][2]-X[j][2],2));
	    			dx = X[i][0]-X[j][0];
	    	    		dy = X[i][1]-X[j][1];
	    	    		dz = X[i][2]-X[j][2];
	    	    		dx -= round(dx / mySys.box_x) * mySys.box_x;
	    	    		dy -= round(dy / mySys.box_y) * mySys.box_y;
	    	    		dz -= round(dz / mySys.box_z) * mySys.box_z;
	    			r = sqrt( pow(dx,2) + pow(dy,2) + pow(dz,2));

	    			if(r != 0 && r < mySys.sigma_cut){
	    				en += 4 * mySys.eps * (pow((mySys.sigma / r),12) - pow((mySys.sigma / r),6));
	    			}
	    			
	    			else if(r >=  mySys.sigma_cut){
	    				en +=  0;
	    			}
			}
		}
	}

	// Add the tail correction term
        en += mySys.NPart * (8.0 / 3.0) * M_PI * mySys.density * (1. / 3. * pow(mySys.sigma / mySys.sigma_cut, 9) - pow(mySys.sigma / mySys.sigma_cut, 3));

	
    return en;

}

double update_energy(double **X, double **X_new, double en_old, int i_star){
	
	double en_new = en_old; // was computed on X
	double en_diff  = 0;
	double dx= 0., dy= 0., dz = 0.;
	
    	// potential energy
    	double r = 0.0;  // distance between particle i and j    	
    	for(int i = 0; i< mySys.NPart; i++){
    	
    	    	dx = X[i][0]-X[i_star][0];
    		dy = X[i][1]-X[i_star][1];
    		dz = X[i][2]-X[i_star][2];
    		
    		dx -= round(dx / mySys.box_x) * mySys.box_x;
    		dy -= round(dy / mySys.box_y) * mySys.box_y;
    		dz -= round(dz / mySys.box_z) * mySys.box_z;

		r = sqrt( pow(dx,2) + pow(dy,2) + pow(dz,2));
    		//cout << r << endl;	
		//r = sqrt( pow(X[i][0]-X[i_star][0],2) + pow(X[i][1]-X[i_star][1],2) + pow(X[i][2]-X[i_star][2],2));
		if (r !=0 && r <= mySys.sigma_cut && r != 0){
			en_new -= 4 * mySys.eps * (pow((mySys.sigma / r),12) - pow((mySys.sigma / r),6));
			en_diff -= 4 * mySys.eps * (pow((mySys.sigma / r),12) - pow((mySys.sigma / r),6));
		}
		else if (r > mySys.sigma_cut){
			en_new -= 0;
			en_diff -=0;
		}
	}
	
	for(int i = 0; i< mySys.NPart; i++){
	
	    	dx = X_new[i][0]-X_new[i_star][0];
    		dy = X_new[i][1]-X_new[i_star][1];
    		dz = X_new[i][2]-X_new[i_star][2];
    		dx -= round(dx / mySys.box_x) * mySys.box_x;
    		dy -= round(dy / mySys.box_y) * mySys.box_y;
    		dz -= round(dz / mySys.box_z) * mySys.box_z;
		r = sqrt( pow(dx,2) + pow(dy,2) + pow(dz,2));
		
		//r = sqrt( pow(X_new[i][0]-X_new[i_star][0],2) + pow(X_new[i][1]-X_new[i_star][1],2) + pow(X_new[i][2]-X_new[i_star][2],2));
		if (r !=0 && r <= mySys.sigma_cut){
			en_new += 4 * mySys.eps * (pow((mySys.sigma / r),12) - pow((mySys.sigma / r),6));
			en_diff += 4 * mySys.eps * (pow((mySys.sigma / r),12) - pow((mySys.sigma / r),6));
		}
		else if (r > mySys.sigma_cut){
			en_new += 0;
			en_diff += 0;
		}
	}
	
	
    return en_diff;

}
