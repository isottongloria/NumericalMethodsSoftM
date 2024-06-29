#include <iostream>
using namespace std;


double compute_Pressure(double **X){

	double virial = 0.0;          // const + virial term
	double pressure = 0.0 ;  // tail correction 
	double dx= 0., dy= 0., dz = 0.;
    	// lennard jones potential energy
    	double rij = 0.0;  // distance between particle i and j    	
    	for(int i = 0; i< mySys.NPart; i++){
    	    	for(int j = 0; j <i; j++){
    	    		if (j != i) {
	    	    		dx = X[i][0]-X[j][0];
	    	    		dy = X[i][1]-X[j][1];
	    	    		dz = X[i][2]-X[j][2];
	    	    		dx -= round(dx / mySys.box_x) * mySys.box_x;
	    	    		dy -= round(dy / mySys.box_y) * mySys.box_y;
	    	    		dz -= round(dz / mySys.box_z) * mySys.box_z;
	    			rij = sqrt( pow(dx,2) + pow(dy,2) + pow(dz,2));
    
		    		virial +=  48 * mySys.eps * ( 0.5*pow(mySys.sigma / rij,6) -  pow(mySys.sigma / rij, 12));

	    		}
	    	}
	}
	
    	//Pres += mySys.NPart * mySys.T / pow(mySys.box_x,3); //  
	pressure = (mySys.NPart*mySys.T - virial/3.)/(mySys.box_x*mySys.box_y*mySys.box_z);
    
    
	// Add the tail correction term
        pressure += (16.0 / 3.0) * M_PI * pow(mySys.density, 2.) * (2. / 3. * pow(mySys.sigma / mySys.sigma_cut, 9.) - pow(mySys.sigma / mySys.sigma_cut, 3.));

    	// cout << mySys.density << "  " << rij << "   " << mySys.eps << endl;
    return pressure;
}



