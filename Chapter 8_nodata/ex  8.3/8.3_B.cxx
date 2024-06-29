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

using namespace std;


void init_Positions(double **X) {
	
	// Initialize the matrix X with ordered values for positions x, y, z
	
	// starting points for placing the particles along x (5 particles), y (4 part.), z (4 part.)
	int start_x = static_cast<int>((mySys.box_x - 5 * mySys.radius) / 2);
	int start_y = static_cast<int>((mySys.box_x - 4 * mySys.radius) / 2);
	int start_z = static_cast<int>((mySys.box_x - 4 * mySys.radius) / 2);

	// number of particles to put on each axes
	int nx = 5;
	int ny = 5;
	int nz = 4;
	
	// place the particles at distance mySys.radius, not overlapped
	int p_count = 0; // particle identifier
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			for (int k = 0; k < nz; k++) {
				// Calculate initial position
				X[p_count][0] = start_x + i* mySys.radius;
				X[p_count][1] = start_y + j* mySys.radius;
				X[p_count][2] = start_z + k* mySys.radius;
				
				p_count += 1;
			}
		}
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





int main(int argc, const char * argv[]) {

	const double K_B = 8.617343e-5; // Boltzmann's constant in natural uni
	double m = 1 ; 
	mySys.NPart = 100;     // number of particles
	mySys.disp_max =1.0; // max displacement
	mySys.NSteps = 6000;   // monte carlo steps
	mySys.realiz = 10;    // MC realizations
	mySys.T = 2;          // temperature
	mySys.energy = 0.;     // initialize energy to zero
	mySys.energy_curr = 0.;
	mySys.seed = time(0); // seed for random number generator 
	mySys.box_x = 5.0;
	mySys.box_y = 5.0;
	mySys.box_z = 5.0;
	mySys.sigma = 1;
	mySys.eps = 1.0;
	mySys.radius = 1.0; // set hard sphere radius (sigma)
	const double beta = 1.0 / ( K_B * mySys.T );
	int i_star = 0;       // index of particle to displace
	double diff_ene = 0.0; // difference of energy between config
	double u = 0.0;       // random uniform number
	double ratio = 0.0;   // metropolis ratio
	double proposed_position = 0.0; 
	bool has_overlap = false; // Do hard spheres have overlap?
	double r = 0.; // distance between 2 spheres
	int n = 50;               // number of radius points in which I evaluate the radial_distribution
	double dr = 0.1; // for radial distribution function computation
	double accept_ratio = 0; // acceptance ratio
	double temp = 0.0; 
	
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
	mySys.energy = compute_energy(X, m);
	
	// ********* Save data to file for each timestep *******
	ofstream filemc;
	filemc.open("C_rho08_dmax1.txt",ios::out);
	filemc.precision(10);

	ofstream filemc3;
	filemc3.open("C_a_rho08_dmax1.txt",ios::out);
	filemc3.precision(10);

	ofstream filemc4;
	filemc4.open("C_E_rho08_dmax1.txt",ios::out);
	filemc4.precision(10);
		
	// ********* Monte Carlo *********
	time_t startTime = time(0);
	for(int real = 0; real < mySys.realiz; real++){
		
		// ********* Call the initialization function *********
	    	init_Positions(X); 
	    	copy(X,X_new);							
		mySys.energy = compute_energy(X, m);
		
		
		// ********* Record the initialization of the particles *********
		if (real==0){
			filemc << mySys.NPart << endl; //  ovito format requirement
			filemc << endl;                //  ovito format requirement
			for (int n = 0; n < mySys.NPart; n++){
				filemc << "S" << n  << "\t" <<  X[n][0] << "\t" << X[n][1] << "\t" << X[n][2] << endl;
			}
		}
		cout << "realiz -> " << real << endl;
		
		// ********* MC sweeps *********
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

				
				// Check for overlaps			
				has_overlap = false;
				for (int a = 0; a < mySys.NPart; a++) {
					r = sqrt(pow(X_new[i_star][0]-X_new[a][0],2) + pow(X_new[i_star][1]-X_new[a][1],2) + pow(X_new[i_star][2]-X_new[a][2],2));
					if (a != i_star && r < mySys.radius) {
						has_overlap = true;
						//cout << "OVERLAP HERE" << endl;
						break;
					}
				}

							
				// Metropolis test
				//mySys.energy_curr = update_energy(X, X_new, mySys.energy, i_star);
				diff_ene = update_energy(X, X_new, mySys.energy, i_star); // mySys.energy_curr - mySys.energy ; 
				ratio = exp((-diff_ene)/beta);
				u = ran38(&(mySys.seed)); // uniform real in [0,1]
				
				// 1) check overlap condition
				if (has_overlap == true) {
					for (int j = 0; j < 3; j++){
						// restore the X_new array pre-displacement-values 
						temp = X[i_star][j];
						X_new[i_star][j] = temp;
					}				
				
				}	
				
				// 2) if 1)=false (no overlap), check metropolis			
				else if (diff_ene <= 0.0 || u < ratio) {  // has_overlap == false && 

					for (int j = 0; j < 3; j++){
						temp = X_new[i_star][j];
						X[i_star][j] = temp;
					}
					// update energy
					mySys.energy += diff_ene;
					
					// update acceptance ratio
					accept_ratio_t[t] += 1; 
					
				}
			
				// if 1)=false and 2)=false, 
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
			if (real==0){
				filemc << mySys.NPart << endl; //  ovito format requirement
				filemc << endl;                //  ovito format requirement
				for (int n = 0; n < mySys.NPart; n++){
					filemc << "S" << n  << "\t" <<  X[n][0] << "\t" << X[n][1] << "\t" << X[n][2] << endl;
				}
			}
			

		} // close timesteps		
		

	} // close realiz
	time_t endTime = time(0);

	// total simulation time
	cout << "total simulation time" << endTime - startTime << endl;
	

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

