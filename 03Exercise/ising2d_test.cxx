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

void setVectorToOne(double* vector, int length) {
	fill(vector, vector + length, 1.0);
}


double calculateC(double* energies, long int nsteps, double kbT) {
	double sum = 0.0;
	double mean = 0.0;

	// Calculate the mean of energies
	for (int i=0; i < nsteps+1; i++) {
	mean += energies[i];
	}
	mean /= (nsteps+1);

	// Calculate the sum of squared deviations
	for (int i=0; i < nsteps+1; i++) {
	sum += pow(energies[i] - mean, 2);
	}

	// Calculate C
	double C = (sum / (nsteps+1)) / pow(kbT, 2);
	return C;
}

double calculateChi(double* magnetizations, long int nsteps, double kbT) {
	double sum = 0.0;
	double mean = 0.0;

	// Calculate the mean of magnetizations
	for (int i=0; i < nsteps+1; i++) {
	mean += magnetizations[i];
	}
	mean /= (nsteps+1);

	// Calculate the sum of squared deviations
	for(int i=0; i < nsteps+1; i++) {
	sum += pow(magnetizations[i] - mean, 2);
	}

	// Calculate Chi
	double Chi = (sum / (nsteps+1)) / kbT;
	return Chi;
}

double error_bar_batch(double* O, int k, double b, int n){
	// this function computes the error on O using batch means
	
	vector<double> Y_k; // vector for the k-batch means
	double batch_mean = 0.0;
	int start = 0;      // auxiliar variable -> marks the point of the vector O to start to do the k-batch mean
	
	for(int i=0; i < k; i++){
	
		for(int j = start; j < start + b; j++){
			batch_mean += O[j] / b;
		}

		start += b; // b is the ith-batch mean
		Y_k.push_back(batch_mean);
		batch_mean = 0.0;
	}
	
	// Calculate the mean of Y_k
	double mean = 0.0;
	for (int i=0; i < k; i++) {
	mean += Y_k[i]/k;
	}

	// Calculate the variance of batck means
	double S2 = 0.0;
	for(int i=0; i < k; i++) {
	S2 += pow(Y_k[i] - mean, 2)/(k-1);
	}

	return sqrt(S2/k);
}


double error_bar_batch2(double* O, int k, double b, int n){
	// this function computes the error on O^2 using batch means
	
	vector<double> O2;  //array of O^2
	vector<double> Y_k; // vector for the k-batch means
	double batch_mean = 0.0;
	int start = 0;      // auxiliar variable -> marks the point of the vector O to start to do the k-batch mean
	
	for(int i=0; i < n + 1; i++){
		O2.push_back(pow(O[i],2));      // computes array of O^2
	}
	
	for(int i=0; i < k; i++){
	
		for(int j = start; j < start + b; j++){
			batch_mean += O2[j] / b;
		}

		start += b; // b is the ith-batch mean
		Y_k.push_back(batch_mean);
		batch_mean = 0.0;
	}
	
	// Calculate the mean of Y_k
	double mean = 0.0;
	for (int i=0; i < k; i++) {
	mean += Y_k[i]/k;
	}

	// Calculate the variance of batck means
	double S2 = 0.0;
	for(int i=0; i < k; i++) {
	S2 += pow(Y_k[i] - mean, 2)/(k-1);
	}

	return sqrt(S2/k);
}

double error_bar_C(double E_ave, double error_E, double error_E2, double kbT){
	// E_ave = <E>    
	// error_E = delta<E>      --> error_bar_batch2 function
	// error_E2 = delta<E^2>   --> error_bar_batch function
	
	double delta = 0.0;
	delta = sqrt( pow(2*E_ave*error_E, 2) + pow(error_E2, 2) ) / ( kbT * kbT );

	return delta;
}

double error_bar_Chi(double M_ave, double error_M, double error_M2, double kbT){
	// M_ave = <M>    
	// error_M = delta<M>      --> error_bar_batch2 function
	// error_M2 = delta<M^2>   --> error_bar_batch function
	
	double delta = 0.0;
	delta = sqrt( pow(2*M_ave*error_M, 2) + pow(error_M2, 2) ) / kbT;

	return delta;
}



int main(int argc, const char * argv[]) {

	int nx = 20;
	int ny = 20;
	long int nsteps = 1e5;
	int burnInSteps = 3000; 
	double J = 1;
	double const T_critic = 2.0/log(1.0 + sqrt(2.0));
	double kbT = 1.5*T_critic;
	double ratio = 0., r = 0.;
	int index = 0, iter_count = 0, realiz = 250;
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
        initializeSpins(nx, ny, spins, index);
    	energy0 = energia(J, nx, ny, spins, index);
    	magnetization0 = magnetization(nx, ny, spins, index);
        E[0] = energy0;
        M[0] = magnetization0;
	
	
	// QUESTA PARTE DI CODICE SERVE PER PRODURRE LA CHAIN
	// Average over realizations
	/*kbT = T_critic;
	setVectorToZero(M, nx * ny);
	setVectorToZero(E, nx * ny);
	
	for(int a = 0; a < realiz; a++){ 
		setVectorToOne(spins, nx * ny); // random ---> initializeSpins(nx, ny, spins, index);
        	
		// Equilibration
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
		cout << "done realiz" << a << endl;
	}	
	

	ofstream filemc;
	filemc.open("100_c.txt",ios::out);
	filemc.precision(10);
	
	for(long int t=0; t < nsteps; t++){
		if(t%500==0) filemc << t << "\t" << E[t]/(nx*ny) << "\t" << M[t]/(nx*ny) <<  '\n';
	}*/
	
	

	// QUESTA PARTE DI CODICE serve per fare FINITE SIZE scaling analysis
	ofstream filemc;
	filemc.open("finite_size_analysis_20_partial2.txt",ios::out);
	filemc.precision(10);
	filemc << "# KbT \t E \t M \t C \t Chi" << endl;
	
	vector<double> kbT_list = {2.51578947,2.55263158, 2.58947368, 2.62631579, 2.66315789, 2.7, 2.75, 2.8, 3}; // Possible temperatures;
	//2.1, 2.03684211, 2.07368421, 2.11052632, 2.14736842, 2.18421053, 2.22105263, 2.25789474, 2.29473684, 2.33157895, 2.36842105, 2.40526316, 2.44210526, 2.47894737,

	for(int a = 0; a < kbT_list.size(); a++){
		kbT = kbT_list[a];
		setVectorToZero(M, nx * ny);
		setVectorToZero(E, nx * ny);
		cout << "Running with T = " << kbT << " ..." << endl;
		for(int b = 0; b < realiz; b++){ 
			setVectorToOne(spins, nx * ny); // random ---> initializeSpins(nx, ny, spins, index);
			
			// Equilibration
			for(long int t = 0; t < burnInSteps; t++){
				spins = glauber(spins, nx, ny, J, kbT, rnd_flip);
	    		}
	    		cout << "Completed Burnt in " << endl;
	    		
			// Trajectory of monte-carlo steps
			for(long int t=0; t < nsteps; t++){
			
				// Keep trace of energy and magnetisation
		    		M[t] = M[t] + magnetization(nx, ny, spins, index)/realiz;
		    		E[t] = E[t] + energia(J, nx, ny, spins, index)/realiz;

		    		// Evolve state at each mc step
		    		spins = glauber (spins, nx, ny, J, kbT, rnd_flip);
				}	
			//cout << "Done realization  " << b << endl;
		}
	
	// compute the mean M, E and Chi,C
	M_ave = 0.0;
	E_ave = 0.0;
	for(long int t=0; t < nsteps; t++){
		M_ave += M[t]/(nx*ny*nsteps);
		E_ave += E[t]/(nx*ny*nsteps);
	}

	// compute susceptibility and specific heat
	C_ave = calculateC(E,nsteps, kbT);
	Chi_ave = calculateChi(M,nsteps, kbT);
	
	// compute obserables per site
	for(long int t=0; t < nsteps; t++){
		E[t] /= (nx*ny);
		M[t] /= (nx*ny);
	}
	C_ave /= (nx*ny);
	Chi_ave/= (nx*ny);
		
	// compute error bars using batch mean method and error propagation
	E_error = error_bar_batch(E, k, b, nsteps); 			//d<E>
	E2_error = error_bar_batch2(E, k, b, nsteps); 			//d<E^2>
	M_error = error_bar_batch(M, k, b, nsteps);  			//d<M>
	M2_error = error_bar_batch2(M, k, b, nsteps); 			//d<M^2>
	C_error = error_bar_C(E_ave, E_error, E2_error, kbT);		//d<C>
	Chi_error = error_bar_Chi(M_ave, M_error, M2_error, kbT);	//d<Chi>
		
	// Store all data in a file
	cout << kbT << "\t" << E_ave << "\t" << E_error << "\t" << M_ave << "\t"  << M_error << "\t"  << C_ave << "\t" << C_error <<  "\t"<<  Chi_ave << "\t" << Chi_error <<  endl;
	
	filemc << kbT << "\t" << E_ave << "\t" << E_error << "\t" << M_ave << "\t"  << M_error << "\t"  << C_ave << "\t" << C_error <<  "\t"<<  Chi_ave << "\t" << Chi_error <<  endl;
	}	



    return 0;
}

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
	
	
	
	
	/*// above Tc
	
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
	*/
	
	
	// free memory
    	/*delete(E);
        delete(M);
        delete(spins);*/
        
