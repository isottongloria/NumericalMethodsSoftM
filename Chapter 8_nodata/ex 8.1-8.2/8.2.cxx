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

void init_Positions(double **X)
{

    // Initialize the matrix with random values for positions x, y, z
    for (int i = 0; i < mySys.NPart; i++)
    {
     	X[i][0] = ran38(&(mySys.seed))*mySys.box_x; // x position
        X[i][1] = ran38(&(mySys.seed))*mySys.box_y; // y position
        X[i][2] = ran38(&(mySys.seed))*mySys.box_z; // z position
    }

}


void init_Velocity(double **V, double K_B, double m)
{

    // Initialize the matrix with random values for positions x, y, z
    for (int i = 0; i < mySys.NPart; i++)
    {
     	V[i][0] = gaussrand(&(mySys.seed))* m / (K_B*mySys.T) ; // x position
        V[i][1] = gaussrand(&(mySys.seed))* m / (K_B*mySys.T) ; // y position
        V[i][2] = gaussrand(&(mySys.seed))* m / (K_B*mySys.T); // z position
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
	mySys.NPart = 10;     // number of particles
	mySys.disp_max = 0.1; // max displacement
	mySys.NSteps = 100;   // monte carlo steps
	mySys.T = 2;          // temperature
	mySys.energy = 0.;     // initialize energy to zero
	mySys.energy_curr = 0.;
	mySys.seed = time(0); // seed for random number generator 
	mySys.box_x = 1.;
	mySys.box_y = 1.;
	mySys.box_z = 1.;
	mySys.sigma = 1.0;
	mySys.eps = 1.0;
	mySys.disp_max = 0.2; // set max displacement
	const double beta = 1.0 / ( K_B * mySys.T );
	int i_star = 0;       // index of particle to displace
	double diff_ene = 0.0; // difference of energy between config
	double u = 0.0;       // random uniform number
	double ratio = 0.0;   // metropolis ratio
	double proposed_position = 0.0; 
	int timesteps = 1; // number of timesteps
	
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
    	
	//  for the velocity (NPartx3) 
	double **V = (double **)malloc(mySys.NPart * sizeof(double *));    	
   	for (int i = 0; i < mySys.NPart; i++){
    	    V[i] = (double *)malloc(3.0 * sizeof(double));
    	}
    	    	
    	    	
    	    	
    	// ********* Call the initialization function *********
    	init_Positions(X); 
    	copy(X,X_new);							
	init_Velocity(V, K_B, m);
	
	
	
	// ********* Compute the energy of the system before displ *******
	mySys.energy = compute_energy(V, m);
	
	
	
	// ********* Monte Carlo *********
	for(int t = 0; t < timesteps; t++){
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
			mySys.energy_curr = compute_energy(V, m);
			diff_ene = mySys.energy_curr - mySys.energy ; 
			ratio = exp((-diff_ene)/beta);
			u = ran38(&(mySys.seed)); // uniform real in [0,1]
				
			if (diff_ene <= 0.0 || u < ratio) {
				for (int j = 0; j < 3; j++){
					X[i_star][j] = X_new[i_star][j];
				}
				// update energy
				mySys.energy = mySys.energy_curr;
			}
			
			else {
				for (int j = 0; j < 3; j++){
					// if the test fails, restore the X_new array pre-displacement-values 
					X_new[i_star][j] = X[i_star][j];
				}

			}
		}
		
	}
	




    	// ********* Deallocate memory for the matrix ********* 
    	for (int i = 0; i <  mySys.NPart; i++)
    	{
        	free(X[i]);
    	}
    	free(X);
	
    return 0;
}	


        // ********* to draw a random number *********
        // double random_value = ran38(&(mySys.seed));
  

/*for (int i = 0; i < mySys.NPart; ++i) {
		for (int j = 0; j < mySys.NPart; ++j) {
			for (int k = 0; k < mySys.NPart; ++k) {
				if (icount < mySys.NPart) {
				    L = mySys.box_x;
				} else if (icount < 2 * mySys.NPart) {
				    L = mySys.box_y;
				} else {
				    L = mySys.box_z;
				}

				x[icount] = distribution2(generator) * L;
				cout << x[icount] << " ";
				icount++;
			    }
			cout << endl;
		}
	        cout << endl;
	    }
	    
	    
	    
	    
	    
	        // Print the matrix for verification
    cout << "Matrix X:" << endl;
    for (int i = 0; i < mySys.NPart; i++)
    {
        for (int j = 0; j < 3.0; j++)
        {
            cout << X[i][j] << " ";
        }
        cout << endl;
    }
    */
