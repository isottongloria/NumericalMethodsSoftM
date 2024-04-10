#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <cstdlib>
#include <ctime>
#include <numeric>

using namespace std;

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


void simulate_ising (double* spins, double* M, double* E, int nx, int ny, double J, double kbT, int rnd_flip, int index, int burnInSteps, int nsteps, double realiz){
	// Average over realizations
	for(int a = 0; a < realiz; a++){ 
	
		// Equilibration
		initializeSpins(nx, ny, spins, index);
		for(long int t = 0; t < burnInSteps; t++){
        		spins = glauber(spins, nx, ny, J, kbT, rnd_flip);
    		}
    		
		// Trajectory of monte-carlo steps
		for(long int t=0; t < nsteps; t++){
		
			// Keep trace of energy and magnetisation
            		M[t] = M[t] + magnetization(nx, ny, spins, index)/realiz;
            		E[t] = E[t] + energia(J, nx, ny, spins, index)/realiz;
            		
            		// Evolve state at each mc step
            		spins = glauber (spins, nx, ny, J, kbT, rnd_flip);
		
		}	
	}	
}	


// Function to calculate the mean of a vector
double mean(double* data) {
    double sum = 0.0;
    double nsteps = sizeof(data);
    for (long int t=0; t < nsteps; t++) {
        sum += data[t];
    }
    return sum / sizeof(data);
}

// Function to calculate the standard deviation of a vector
double stdDev(double* data, double mean) {
    double variance = 0.0;
    double nsteps = sizeof(data);
    for (long int t=0; t < nsteps; t++) {
        variance += (data[t] - mean) * (data[t] - mean);
    }
    return sqrt(variance / (sizeof(data) - 1));
}


double* acf (double* data)
{
    double media = mean(data);
    double* autocorrelation = new double[sizeof(data)/2];
    for (int t = 0; t < sizeof(autocorrelation); t ++)
    {
        double n = 0; // Numerator
        double d = 0; // Denominator
        for (int i = 0; i < sizeof(data); i ++)
        {
            double xim = data[i] - media;
            n += xim * (data[(i + t) % sizeof(data)] - media);
            d += xim * xim;
        }
        autocorrelation[t] = n / d;
    }
    return autocorrelation;
}

int main(int argc, const char * argv[]) {

	int nx = 25;
	int ny = 25;
	long int nsteps = 5000;
	int burnInSteps = 1000; 
	double J = 1.;
	double const T_critic = 2.0/log(1.0 + sqrt(2.0));
	double kbT = 0.8*T_critic;
	double ratio = 0., r = 0.;
	int index = 0, iter_count = 0, realiz = 50;
	double diff_ene;
	int rnd_flip;
	double energy0, magnetization0, magn = 0;
	double E_ave = 0., M_ave = 0., C_ave = 0., Chi_ave = 0.;
	
	cout << "Nx: "<< nx << endl;
        cout << "Nx: "<< ny << endl;
	cout << "J: "<< J << endl;
	cout << "kbT: " << kbT << endl;
	cout << "N passi: " << nsteps << endl;
	cout.precision(10);
	
	// Marsenne Twister generator 
        mt19937 mt_generator(time(0)); 
    	uniform_real_distribution<double> uniform_dist(0.0, 1.0); // r for metropolis hastings
    	uniform_int_distribution<int> uniform_flip(0, nx*ny); // rnd_flip = proposal spin flip
    	
    	// allocate spins, energy and magnetization
	double* E = new double[nsteps];
	double* M = new double[nsteps];
	double* spins = new double[nx * ny];   	
    	
    	//initialize energy and magnetization
        initializeSpins(nx, ny, spins, index);
    	energy0 = energia(J, nx, ny, spins, index);
    	magnetization0 = magnetization(nx, ny, spins, index);
        E[0] = energy0;
        M[0] = magnetization0;

	
	// Simulation
	
	//simulate_ising (spins, M, E, nx, ny, J, kbT, rnd_flip, index, burnInSteps, nsteps, realiz);
	// Average over realizations
	for(int a = 0; a < realiz; a++){ 
	
		// Equilibration
		initializeSpins(nx, ny, spins, index);
		for(long int t = 0; t < burnInSteps; t++){
        		spins = glauber(spins, nx, ny, J, kbT, rnd_flip);
    		}
    		
		// Trajectory of monte-carlo steps
		for(long int t=0; t < nsteps; t++){
		
			// Keep trace of energy and magnetisation
            		M[t] = M[t] + magnetization(nx, ny, spins, index)/realiz;
            		E[t] = E[t] + energia(J, nx, ny, spins, index)/realiz;
            		
            		// Evolve state at each mc step
            		spins = glauber (spins, nx, ny, J, kbT, rnd_flip);
		
		}	
	}
		
	// Compute observables
	
	// Calculate the mean of E and M
	/*E_ave = mean(E)/(nx*ny);
  	M_ave = mean(M)/(nx*ny);

   	// Calculate the standard deviation for E and M
   	double stdDevE = stdDev(E, E_ave)/(nx*ny);
   	double stdDevM = stdDev(M, M_ave)/(nx*ny);

   	// Calculate specific heat (C) and susceptibility (Chi)
    	C_ave = (stdDevE * stdDevE) / (kbT * kbT);
    	Chi_ave = (stdDevM * stdDevM) / kbT;

	cout << "Mean of E: " << E_ave << endl;
	cout << "Mean of M: " << M_ave << endl;
	cout << "Specific Heat (C): " << C_ave << endl;
	cout << "Susceptibility (Chi): " << Chi_ave << endl;
   
    	//compute acf function
        double * autocorr = acf (M);  
        */
                   
                   
                     
	// save data to file
	ofstream filemc;
	filemc.open("ising2d_eqT35.txt",ios::out);
	filemc.precision(10);
	
	for(long int t=0; t < nsteps; t++){
		if(t%200==0) filemc << t << "\t" << E[t]/(nx*ny) << "\t" << M[t]/(nx*ny)  << '\n'; //<< "\t" << autocorr[t] 
	}
	
	// free memory
    	delete(E);
        delete(M);
        delete(spins);

    return 0;
}
