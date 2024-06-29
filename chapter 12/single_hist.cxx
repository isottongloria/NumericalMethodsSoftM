#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <cstdlib>
#include <ctime>
#include <numeric>
#include <vector>
#include <iomanip>

using namespace std;


void initializeSpins(int nx, int ny, double *spins, int index, mt19937& mt_generator, uniform_real_distribution<double>& uniform_dist)
{

    	double r = 0.0;
    	for (int i = 0; i < nx * ny; i++) {
    		//r = pseudo_random();
    		r = uniform_dist(mt_generator);
    		if(r<0.5){
        		spins[i]=-1;
    		}else{
        		spins[i]=1;
    		}
      
    	}
    
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



double calculateC(double* energies, long int nsteps, double kbT) {
	// Calculate the average energy
	double avgEnergy = 0.0;
	for (long int i = 0; i < nsteps; ++i) {
		avgEnergy += energies[i];
	}
	avgEnergy /= nsteps;

	// Calculate the variance
	double variance = 0.0;
	for (long int i = 0; i < nsteps; ++i) {
		variance += (energies[i] - avgEnergy) * (energies[i] - avgEnergy);
	}
	variance /= nsteps;

	// Calculate specific heat
	return variance / pow(kbT,2) ;
}

double calculateChi(double* magnetizations, long int nsteps, double kbT) {
	// Calculate the average energy
	double avgMagnet = 0.0;
	for (long int i = 0; i < nsteps; ++i) {
		avgMagnet += magnetizations[i];
	}
	avgMagnet /= nsteps;

	// Calculate the variance
	double variance = 0.0;
	for (long int i = 0; i < nsteps; ++i) {
		variance += (magnetizations[i] - avgMagnet) * (magnetizations[i] - avgMagnet);
	}
	variance /= nsteps;

	// Calculate specific heat
	return variance / kbT ;
}




int main(int argc, const char * argv[]) {

	int nx = 50;
	int ny = 50;
	long int nsteps = 1e4;
	int burnInSteps = 0; 
	double J = 1;
	double const T_critic = 2.0/log(1.0 + sqrt(2.0));
	double kbT = 0.5*T_critic;
	double ratio = 0., r = 0.;
	int index = 0, iter_count = 0, realiz = 1;
	double diff_ene;
	int rnd_flip;
	double energy0, magnetization0, magn = 0;
	double E_ave = 0., M_ave = 0., C_ave = 0., Chi_ave = 0.;
	double E_error = 0., E2_error =0., M_error = 0., M2_error =0., C_error = 0., Chi_error = 0.;
	double b = 100; // batch size
	double k = static_cast<int>(nsteps / b); // number of batches
	
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
	double* E = new double[nsteps+1];
	double* M = new double[nsteps+1];
	double* spins = new double[nx * ny];   	
    	
    	//initialize energy and magnetization
        initializeSpins(nx, ny, spins, index, mt_generator, uniform_dist);
    	energy0 = energia(J, nx, ny, spins, index);
    	magnetization0 = magnetization(nx, ny, spins, index);
        E[0] = energy0;
        M[0] = magnetization0;
	

	// define a vector with all the temperatures of the simulations
	//vector<double> kbT_list = {2.26918531, 2.45828409, 2.64738287, 2.83648164, 3.02558042, 3.2146792 , 3.40377797};
	vector<double> kbT_list = {2.0421, 2.269 , 2.4959, 2.7228, 2.9497, 3.1766, 3.4035};

	// save in separate files the trajectories
	// Save in separate files the trajectories
	vector<string> filenames = {
	"single_hist_1.txt", "single_hist_2.txt", "single_hist_3.txt",
	"single_hist_4.txt", "single_hist_5.txt", "single_hist_6.txt",
	"single_hist_7.txt"
	};

	vector<ofstream*> fileStreams(kbT_list.size());

	for (int i = 0; i < kbT_list.size(); ++i) {
		fileStreams[i] = new ofstream(filenames[i]);
		if (!fileStreams[i]->is_open()) {
		    cerr << "Error: Unable to open file " << filenames[i] << endl;
		    return 1;
		}
		fileStreams[i]->precision(10);
	}
	
	cout << "ok" << endl;

	
	for(int Tt = 0; Tt < kbT_list.size(); Tt++){ 
		
		// choose one temperature
		kbT = kbT_list[Tt];
		 
		// reset M, E
		setVectorToZero(M, nx * ny);
		setVectorToZero(E, nx * ny);
		
		// initialize spins
		setVectorToOne(spins, nx * ny); // random ---> initializeSpins(nx, ny, spins, index);
		
		// Equilibration
		for(long int t = 0; t < burnInSteps; t++){
			spins = glauber(spins, nx, ny, J, kbT, rnd_flip, mt_generator, uniform_dist, uniform_flip);
    		}

		// Monte-carlo steps
		for(long int t=0; t < nsteps; t++){
			if (t%1000 == 0) {cout << t << endl;}
			// Keep trace of energy and magnetisation
	    		M[t] =  magnetization(nx, ny, spins, index); 
	    		E[t] =  energia(J, nx, ny, spins, index);

	    		// Evolve state at each mc step 
	    		spins = glauber(spins, nx, ny, J, kbT, rnd_flip,  mt_generator, uniform_dist, uniform_flip);
	    		
		}	
		
		// write the 3 simulations (at 3 different temperatures) on 3 separate files
		
		
		for (long int t = 0; t < nsteps; t++) {
		    *(fileStreams[Tt]) << t << "\t" << E[t] / (nx * ny) << "\t" << M[t] / (nx * ny) << std::endl;
		}

	
	 }
	 
	// deallocate 
	for (int i = 0; i < kbT_list.size(); ++i) {
		fileStreams[i]->close();
		delete fileStreams[i];
	}
    
    return 0;
}


