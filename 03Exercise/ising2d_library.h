#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <cstdlib>
#include <ctime>
#include <numeric>

void initializeSpins(int nx, int ny, double* spins, int index) {
    std::mt19937 mt_generator(std::time(0)); // Mersenne Twister generator (with period 2^19937 - 1)
    std::uniform_int_distribution<int> init_spins(1, 2); // 1 --> 1, 2 --> -1
    index = 0;

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            int random_number = init_spins(mt_generator);
            if (random_number == 1) {
                spins[index] = 1.0;
            } else {
                spins[index] = -1.0;
            }
            index++;
        }
    }
}


double magnetization(int nx, int ny, double *spins, int index){   
	index=0;
        double magn=0;
        
        for (int j=0;j<ny; j++){
        	for(int i=0; i<nx; i++){
            		magn+=spins[index];
            		index++;
        	}
    	}
    return magn;
}


double energia(double J, int nx, int ny, double* spins, int index)
{
    double ene=0.;
    int icount=0;
    int nn_left,nn_right,nn_up,nn_down; //nearest neighbour indexes
    double nn ; 
    
    for(int j=0;j<ny;j++){
        for( int i=0; i<nx;i++){
 		nn_left = ((i - 1 + nx) % nx) + j * nx; //using periodic boundary conditions
		nn_right = ((i + 1) % nx) + j * nx;
		nn_up = i + ((j + 1) % ny) * nx;
		nn_down = i + ((j - 1 + ny) % ny) * nx;
        	nn = spins[nn_left] + spins[nn_right] + spins[nn_up] + spins[nn_down];
        	ene +=  -0.5*J*spins[i + j*nx] * nn;
        }
    }
    
    return ene/2.0;
}


double delta_E(int nx, int ny, double J, double* spins, int rnd_flip){
	int i = rnd_flip % nx; //compute row and column indexes from rnd_flip in [0,nx*ny)
	int j = floor(rnd_flip / nx);
	int nn_left,nn_right,nn_up,nn_down; //nearest neighbour indexes
        nn_left = ((i - 1 + nx) % nx) + j * nx; //using periodic boundary conditions
	nn_right = ((i + 1) % nx) + j * nx;
	nn_up = i + ((j + 1) % ny) * nx;
	nn_down = i + ((j - 1 + ny) % ny) * nx;
	return 2*J*spins[rnd_flip]*(spins[nn_left]+spins[nn_right]+spins[nn_up]+spins[nn_down]);
}


double* update_state (double* spins, int rnd_flip){
	spins[rnd_flip]*=-1.;
	return spins;
}


double* glauber (double* spins, int nx, int ny, double J, double kbT, int rnd_flip){
	double ratio = 0., r = 0., diff_ene=0;
	mt19937 mt_generator(time(0)); 
    	uniform_real_distribution<double> uniform_dist(0.0, 1.0); // r for metropolis hastings
    	uniform_int_distribution<int> uniform_flip(0, nx*ny); // rnd_flip = proposal spin flip
    	
        for(long int n=0; n < nx*ny; n++){
        	rnd_flip = uniform_flip(mt_generator); //proposed spin flip 
        	diff_ene = delta_E(nx, ny, J, spins, rnd_flip); // delta energy if I accept that flip
        	ratio = exp((-diff_ene)/kbT);
        	
        	r = uniform_dist(mt_generator); // uniform real number in [0,1]
		
		if (diff_ene <= 0.0 || r < ratio) {
			spins = update_state(spins, rnd_flip);
		    }
	}
	return spins;
}

