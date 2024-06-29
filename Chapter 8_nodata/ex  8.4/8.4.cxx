/*  Off lattice 3D Monte Carlo simulation
* cd /home/gloria/Scaricati/ovito-basic-3.9.4-x86_64/bin
* ./ovito
*
*
*/
#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <cstdlib>
#include <ctime>
#include <numeric>
#include <vector>

#include "definitions.h"
#include "random.h"
#include "energy.h"
#include "pressure2.h"

using namespace std;


void init_Positions(double **X) {
	
	// Initialize the matrix with random values for positions x, y, z
    	for (int i = 0; i < mySys.NPart; i++) {

	    X[i][0] = ran38(&(mySys.seed)) * mySys.box_x; // x position
	    X[i][1] = ran38(&(mySys.seed)) * mySys.box_y; // y position
	    X[i][2] = ran38(&(mySys.seed)) * mySys.box_z; // z position

	}  

}


void copy(double** X, double** X_new) {

	// Function to copy the content of one 2D array to another
    	for (int i = 0; i < mySys.NPart; i++) {
        	for (int j = 0; j < 3; j++) {
            		X_new[i][j] = X[i][j]; 
        	}
    	}
    
}



// Function to initialize a vector to zero
void setVectorToZero(vector<double>& vec) {
    fill(vec.begin(), vec.end(), 0.0);
}

int main(int argc, const char * argv[]) {

	double m = 1 ; 
	mySys.NPart = 125;     // number of particles
	mySys.disp_max = 0.05; // max displacement
	mySys.NSteps = 100;   // monte carlo steps
	mySys.realiz = 100;    // MC realizations
	mySys.T = 0.9;          // temperature
	mySys.energy = 0.;     // initialize energy to zero
	mySys.energy_curr = 0.;
	mySys.seed = time(0); // seed for random number generator 
	mySys.box_x = 5;
	mySys.box_y = 5;
	mySys.box_z = 5;
	mySys.sigma = 1.0;
	mySys.sigma_cut = mySys.box_x / 2.0;
	mySys.eps = 1.0;
	mySys.density = mySys.NPart / pow(mySys.box_x,3);
	mySys.radius = 1.0001; // set hard sphere radius (sigma)
	const double beta = 1.0 / ( mySys.T );
	int i_star = 0;       // index of particle to displace
	double diff_ene = 0.0; // difference of energy between config
	double u = 0.0;       // random uniform number
	double ratio = 0.0;   // metropolis ratio
	double proposed_position = 0.0; 
	bool has_overlap = false; // Do hard spheres have overlap?
	double r = 0.; // distance between 2 spheres
	int n = 100;               // number of radius points in which I evaluate the radial_distribution
	double dr = 0.1; // for radial distribution function computation
	double accept_ratio = 0; // acceptance ratio
	double temp = 0.0; 
	double P = 0. ;  // pressure

	// ********* Allocate memory *********
	
	//  for the position (NPartx3) 
	double **X = (double **)malloc(mySys.NPart * sizeof(double *));    	
   	for (int i = 0; i < mySys.NPart; i++){
    	    X[i] = (double *)malloc(3.0 * sizeof(double));
    	}
    	
    	// copy of the position (for the proposed displacemnet)
	double **X_new = (double **)malloc(mySys.NPart * sizeof(double *));    	
   	for (int i = 0; i < mySys.NPart; i++){
    	    X_new[i] = (double *)malloc(3.0 * sizeof(double));
    	}
    	
    	    	
    	// **** RDF results **** 
    	vector<double> rho_result(n,0.0);   	

    	// **** Acceptance ratio per timestep, avg over realizations **** 
    	vector<double> accept_ratio_t(mySys.NSteps,0.0);   	
    	    	
     	// **** Energy per timestep, avg over realizations **** 
    	vector<double> energy_t(mySys.NSteps,0.0);   	
    	       	   
	    	    	
    	// ********* Call the initialization function *********
    	init_Positions(X); 
    	copy(X,X_new);							
	
	
	
	// ********* Compute the energy of the system before displ *******
	mySys.energy = compute_energy(X);
	
	// ********* Save data to file for each timestep *******
	ofstream filemc;
	filemc.open("LJ_T2_dmax1.txt",ios::out);
	filemc.precision(10);

	ofstream filemc3;
	filemc3.open("LJ_a_T2_dmax1.txt",ios::out);
	filemc3.precision(10);

	ofstream filemc4;
	filemc4.open("LJ_E_T2_dmax1.txt",ios::out);
	filemc4.precision(10);

	ofstream filemc5;
	filemc5.open("LJ_P_T2_dmax1.txt",ios::out);
	filemc5.precision(10);		

	mySys.NPart = 20;
	for(int nP = 0; nP < 10; nP++){
	
		// update density
		mySys.NPart +=  2*nP;
		mySys.density = mySys.NPart / pow(mySys.box_x,3);
		cout << " Running with density = " << mySys.NPart << "  " << mySys.density << endl;
	    	
		//init pressure
		P = 0. ; 
		
		// set energy to zero
		setVectorToZero(energy_t);
		
		// ********* Monte Carlo *********
		time_t startTime = time(0);
		for(int real = 0; real < mySys.realiz; real++){
			
			// ********* Call the initialization function *********
		    	init_Positions(X); 
		    	copy(X,X_new);							
			mySys.energy = compute_energy(X);

			cout << "realiz -> " << real << endl;
			for(int t = 0; t < mySys.NSteps; t++){
				
				for (int i = 0; i < mySys.NPart; i++){
				
					// Pick one rnd particle
					i_star = ran38(&(mySys.seed)) * mySys.NPart; 

					// For each direction propose a displacement rnd in [- disp_max, disp_max]
					for (int j = 0; j < 3; j++){
					    	proposed_position = X[i_star][j] - mySys.disp_max + ran38(&(mySys.seed)) * (2.0 * mySys.disp_max);

					    	// Apply PBC
				    		X_new[i_star][j] = proposed_position  - floor(proposed_position / mySys.box_x) * mySys.box_x;

					}

								
					// Metropolis test
					//mySys.energy_curr = update_energy(X, X_new, mySys.energy, i_star);
					diff_ene = update_energy(X, X_new, mySys.energy, i_star); // mySys.energy_curr - mySys.energy ; 

					ratio = exp((-diff_ene)/beta);
					u = ran38(&(mySys.seed)); // uniform real in [0,1]	
					
					// 1) = check metropolis			
					if (diff_ene <= 0.0 || u < ratio) {  

						for (int j = 0; j < 3; j++){
							temp = X_new[i_star][j];
							X[i_star][j] = temp;
						}
						// update energy
						mySys.energy += diff_ene;
						
						// update acceptance ratio
						accept_ratio_t[t] += 1; 
						
					}
				
					// if 1) = false 
					else {
						for (int j = 0; j < 3; j++){
							// if the test fails, restore the X_new array pre-displacement-values 
							temp = X[i_star][j];
							X_new[i_star][j] = temp;
						}

					}				
					
				}
				
				// keep track of the new energy at every timestep
				energy_t[t] += mySys.energy / mySys.realiz;
				
					
				// write new configuration to file
				if (real==0 && nP==2 ){
					filemc << mySys.NPart << endl; //  ovito format requirement
					filemc << endl;                //  ovito format requirement
					for (int n = 0; n < mySys.NPart; n++){
						filemc << "S" << n  << "\t" <<  X[n][0] << "\t" << X[n][1] << "\t" << X[n][2] << endl;
					}
				}
				

			} // close timesteps		
			

			// compute pressure at the end of every realization	
			P += compute_Pressure(X) / mySys.realiz;

		} // close realiz
		time_t endTime = time(0);

		// total simulation time
		cout << "total simulation time: " << endTime - startTime << endl;
		
		// update pressure
		filemc5 << mySys.density << "\t" << P << endl;
		cout << "pressure -> " << P << endl;
	} // close density loop
	


	// write acceptance ratio to file 
	for (int i = 0; i < mySys.NSteps ; i++){
		filemc3 << accept_ratio_t[i]/ (mySys.realiz * mySys.NPart) << endl;
	}

	// write energy per timestep to file 
	for (int i = 0; i < mySys.NSteps ; i++){
		filemc4 << energy_t[i] << endl;
	}

	
    	// ********* Deallocate memory for the matrix ********* 
    	for (int i = 0; i <  mySys.NPart; i++)
    	{
        	free(X[i]);
    	}
    	free(X);
	
    return 0;
}	


