#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <cstdlib>
#include <ctime>
#include <numeric>

using namespace std;


void initializeSpins(int nx, int ny, double *spins, int index, mt19937& mt_generator, uniform_real_distribution<double>& uniform_dist)
{

    	double r = 0.0;
    	for (int i = 0; i < nx; i++) {
        	for (int j = 0; j < ny; j++) {
            		//r = pseudo_random();
            		r = uniform_dist(mt_generator);
            		if(r<0.5){
                		spins[index]=-1;
            		}else{
                		spins[index]=1;
            		}
            		index++;
        	}
    	}
    
}


vector<double> init_beta_list(double beta_i, double beta_f, int n_chains){
	// initializes a vector of doubles containing evenly spaced values between beta_i and beta_f
	vector<double> beta_list;
	double delta = (beta_f - beta_i) / n_chains; //step size delta between each element
	
	for (int i = 0; i < n_chains; i++) {
		beta_list.push_back(beta_i + i * delta);
	}
	
	return beta_list;
}


double magnetization(int nx, int ny, double *spins, int index){   
	index=0;
        double magn=0;
        
    	for (index = 0; index < nx * ny; ++index) {
        	magn += spins[index];
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
    
    return ene;
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


double* glauber (double* spins, int nx, int ny, double J, double kbT, int rnd_flip, mt19937& mt_generator, uniform_real_distribution<double>& uniform_dist, uniform_int_distribution<int>& uniform_flip){
	double ratio = 0., r = 0., diff_ene=0;

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

    
// Function to initialize a vector to zero
void setVectorToZero(double* vector, int length) {
	fill(vector, vector + length, 0.0);
}

void setVectorToOne(double* vector, int length) {
	fill(vector, vector + length, 1.0);
}

void setVectorToFalse(vector<bool>& vector, int length) {
    fill(vector.begin(), vector.begin() + length, false);
}

int main(int argc, const char * argv[]) {

	int nx = 120;
	int ny = 120;
	long int nsteps = 1e3;
	int burnInSteps = 1e3; 
	double J = 1;
	double const T_critic = 2.0/log(1.0 + sqrt(2.0));
	double kbT = 1.01*T_critic;
	double ratio = 0., r = 0.;
	int index = 0, iter_count = 0, realiz = 1;
	double diff_ene;
	int rnd_flip;
	double energy0, magnetization0, magn = 0;
	double E_ave = 0., M_ave = 0., C_ave = 0., Chi_ave = 0.;
	double E_error = 0., E2_error =0., M_error = 0., M2_error =0., C_error = 0., Chi_error = 0.;
	double b = 100; // batch size
	double k = static_cast<int>(nsteps / b); // number of batches
	int n_chains = 30;  // number of temperatures considered
	double beta_i = 1. / (1.1*T_critic), beta_f = 1. / (0.9*T_critic); // first and last point of the set of temperature array
	double acceptance_prob = 0.; // swap config acceptance probability
	int n_swap = 0; // number of swaps
	
	// utils
	double u = 0.; // uniform rnd number 	
	int k1 = 0, k2 = 0; // temperatures of adjacent chains
	
	cout << "Nx: "<< nx << endl;
        cout << "Ny: "<< ny << endl;
	cout << "beta : [" << beta_i << ", " << beta_f << "]" << endl;
	cout << "N passi: " << nsteps << endl;
	cout << "N burnt in: " << burnInSteps << endl;
	cout << " Number of chains: " << n_chains << endl;
	cout.precision(10);
	
	// Marsenne Twister generator 
        mt19937 mt_generator(time(0)); 
    	uniform_real_distribution<double> uniform_dist(0.0, 1.0); // r for metropolis hastings
    	uniform_int_distribution<int> uniform_flip(0, nx*ny); // rnd_flip = proposal spin flip
    	uniform_int_distribution<int> uniform_swap(0, n_chains-1); // k1 = proposal chain 1 swap with adjacent chain
    	
	// Allocate memory for energies, magnetizations, spins for each chain
	vector<double*> E(n_chains, nullptr);
	vector<double*> M(n_chains, nullptr);
	vector<double*> spins(n_chains, nullptr);

	// Initialize chains
	for (int i = 0; i < n_chains; ++i) {
		E[i] = new double[burnInSteps + nsteps + 1];
		M[i] = new double[burnInSteps + nsteps + 1];
		spins[i] = new double[nx * ny];
		
		// initialize M,E, rnd spins for every chains
		setVectorToZero(M[i], nx * ny);
		setVectorToZero(E[i], nx * ny);
		initializeSpins(nx, ny, spins[i], index, mt_generator, uniform_dist);
	}
    	
        
	// inverse temperature list
	vector<double> beta_list = {0.56523191, 0.53499243, 0.50782422, 0.48328199, 0.46100257, 0.44068679, 0.42208602, 0.40499188, 0.38922844, 0.37464615}; // list containing Tc
	
	
	// Define a vector to track whether each chain has already performed a swap
        vector<bool> chain_swapped(n_chains, false); // false --> not swapped

	// Swapping rate between the chains
	double* swap_rate = new double[burnInSteps + nsteps + 1];
	setVectorToZero(swap_rate, burnInSteps + nsteps + 1);
	
	// Initialize list of pointers for configurations
        vector<double*> configPointers(n_chains, nullptr);
        for (int i = 0; i < n_chains; ++i) {
        	configPointers[i] = spins[i];
    	}
    		    	 	
	// Save data to files
	ofstream filemc1;
	ofstream filemc2;
	filemc1.open("E_noswap_120_Tc.txt",ios::out);
	filemc1.precision(10);
	filemc2.open("M_noswap_120_Tc.txt",ios::out);
	filemc2.precision(10);	
	
	cout << "choosing temperature to save: " << 1/beta_list[6] << endl;
	    	 	
	// Burnt in steps for each chain
        for(long int t=0; t < burnInSteps; t++) {
        	if (t%100 == 0) { cout << " Iteration --> " << t << endl;}
		for (int k = 0; k < n_chains; ++k) {
			// MC step, for each chain k
			kbT = 1. / beta_list[k];
			configPointers[k] = glauber(configPointers[k], nx, ny, J, kbT, rnd_flip,  mt_generator, uniform_dist, uniform_flip);
			M[k][t] = magnetization(nx, ny, configPointers[k], index);
	    		E[k][t] = energia(J, nx, ny, configPointers[k], index);		
		}
		
		// Save data to file
		filemc1 << t << "\t";
                filemc2 << t << "\t";
		for (int i = 0; i < n_chains; ++i) {
            		filemc1 << E[i][t] / (nx * ny);
            		filemc2 << M[i][t] / (nx * ny);
            		if (i != n_chains-1) {
                		filemc1 << "\t";
                		filemc2 << "\t";
            		} else {
                		filemc1 << "\n";
                		filemc2 << "\n";
            		}
        	}
	
	}	
	
	cout << " ***** end burnt in *****" << endl;
        // Monte Carlo steps for each chain 
        for(long int t=burnInSteps; t < nsteps + burnInSteps; t++) {
        	if (t%100 == 0) { cout << " Iteration --> " << t << endl;} 
		for (int k = 0; k < n_chains; ++k) {
			// MC step, for each chain k
			kbT = 1. / beta_list[k];
			configPointers[k] = glauber(configPointers[k], nx, ny, J, kbT, rnd_flip,  mt_generator, uniform_dist, uniform_flip);
			M[k][t] = magnetization(nx, ny, configPointers[k], index);
	    		E[k][t] = energia(J, nx, ny, configPointers[k], index);

			
		}
		
		// Save data to file
		filemc1 << t << "\t";
                filemc2 << t << "\t";
		for (int i = 0; i < n_chains; ++i) {
            		filemc1 << E[i][t] / (nx * ny);
            		filemc2 << M[i][t] / (nx * ny);
            		if (i != n_chains-1) {
                		filemc1 << "\t";
                		filemc2 << "\t";
            		} else {
                		filemc1 << "\n";
                		filemc2 << "\n";
            		}
        	}       	
				
	}
	


	return 0;
}


