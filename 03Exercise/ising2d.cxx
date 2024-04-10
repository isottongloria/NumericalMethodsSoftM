#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <cstdlib>
#include <ctime>
#include <numeric>

using namespace std;

double pseudo_random(void){
    const long int a = 16807;
    const long int c = 0;
    const long m = 2147483647;
    static long int x0 = time(0);
    long int x1;
    double r;
    
    x1 = ( a * x0 + c ) % m;
    r = ((double) x1)/((double) m);
    x0 = x1;
    
    return r;
    
}


void initializeSpins(int nx, int ny, double *spins, int index)
{
    double r;
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            r = pseudo_random();
            if(r<0.5){
                spins[index]=-1;
            }else{
                spins[index]=1;
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

// Function to initialize a vector to zero
void setVectorToZero(double* vector, int length) {
    fill(vector, vector + length, 0.0);
}

void simulate_ising (double* spins, double* M, double* E, int nx, int ny, double J, double kbT, int rnd_flip, int index, int burnInSteps, int nsteps, double realiz){
	// Average over realizations
	for(int a = 0; a < realiz; a++){ 
		setVectorToZero(M, nx * ny);
		setVectorToZero(E, nx * ny);
		
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



int main(int argc, const char * argv[]) {

	int nx = 25;
	int ny = 25;
	long int nsteps = 1e6;
	int burnInSteps = 0; 
	double J = 1;
	double const T_critic = 2.0/log(1.0 + sqrt(2.0));
	double kbT = 1.5*T_critic;
	double ratio = 0., r = 0.;
	int index = 0, iter_count = 0, realiz = 70;
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
        /*initializeSpins(nx, ny, spins, index);
    	energy0 = energia(J, nx, ny, spins, index);
    	magnetization0 = magnetization(nx, ny, spins, index);
        E[0] = energy0;
        M[0] = magnetization0;
        cout << magnetization0/(250*250) << endl;
	for(long int n=0; n < 20; n++){
		cout << pseudo_random() << endl;
	}*/
	
	// Equilibration
	/*ofstream filemc;
	filemc.open("ising2d_eqT35.txt",ios::out);
	filemc.precision(10);
        for(long int n=0; n < nsteps; n++){
        	
        	rnd_flip = uniform_flip(mt_generator); //proposed spin flip 
        	diff_ene = delta_E(nx, ny, J, spins, rnd_flip); // delta energy if I accept that flip
        	ratio = exp((-diff_ene)/kbT);
        	
        	r = uniform_dist(mt_generator); // uniform real number in [0,1]
		
		if (diff_ene <= 0.0 || r < ratio) {
			spins = update_state (spins, rnd_flip);
	        	E[iter_count] = E[iter_count-1] + diff_ene;
                        M[iter_count] = magnetization(nx, ny, spins, index);	
                        if(iter_count%50==0) filemc << iter_count << "\t" << E[iter_count]/(nx*ny) << "\t" << M[iter_count]/(nx*ny) << '\n';
                        iter_count ++;
		    }
	}
	cout << "Acceptance rate: " << static_cast<float>(iter_count) << endl;

        filemc.close();*/
	
	
	// Simulation

	// Average over realizations
	/*for(int a = 0; a < realiz; a++){ 
	
		// Equilibration
		initializeSpins(nx, ny, spins, index);
		for(long int t = 0; t < burnInSteps; t++){
        		spins = glauber(spins, nx, ny, J, kbT, rnd_flip);
    		}
    		
		// Trajectory of monte-carlo steps
		for(long int t=0; t < nsteps; t++){
		
			// Keep trace of energy and magnetisation
			magn = magnetization(nx, ny, spins, index);
            		M[t] = M[t] + magn/realiz;
            		E[t] = E[t] + energia(J, nx, ny, spins, index)/realiz;
            		
            		// Evolve state at each mc step
            		spins = glauber (spins, nx, ny, J, kbT, rnd_flip);
		
		}
	
	}*/
	
	// above Tc
	
	kbT = 1.5*T_critic;
	simulate_ising (spins, M, E, nx, ny, J, kbT, rnd_flip, index, burnInSteps, nsteps, realiz);
	
	ofstream filemc;
	filemc.open("25_a.txt",ios::out);
	filemc.precision(10);
	
	for(long int t=0; t < nsteps; t++){
		if(t%50==0) filemc << t << "\t" << E[t]/(nx*ny) << "\t" << M[t]/(nx*ny) <<  '\n';
	}
	
	// below Tc
	kbT = 0.5*T_critic;
	simulate_ising (spins, M, E, nx, ny, J, kbT, rnd_flip, index, burnInSteps, nsteps, realiz);
	
	ofstream filemc2;
	filemc2.open("25_b.txt",ios::out);
	filemc2.precision(10);
	
	for(long int t=0; t < nsteps; t++){
		if(t%50==0) filemc2 << t << "\t" << E[t]/(nx*ny) << "\t" << M[t]/(nx*ny) <<  '\n';
	}
	
	// at Tc
	kbT = 1.0*T_critic;
	simulate_ising (spins, M, E, nx, ny, J, kbT, rnd_flip, index, burnInSteps, nsteps, realiz);
	
	ofstream filemc3;
	filemc3.open("25_c.txt",ios::out);
	filemc3.precision(10);
	
	for(long int t=0; t < nsteps; t++){
		if(t%50==0) filemc3 << t << "\t" << E[t]/(nx*ny) << "\t" << M[t]/(nx*ny) <<  '\n';
	}
	
	
	
	// free memory
    	delete(E);
        delete(M);
        delete(spins);

    return 0;
}
