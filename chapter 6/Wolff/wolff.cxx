#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <cstdlib>
#include <ctime>
#include <numeric>
#include <vector>
#include <tuple>


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

bool notInsideVector(int* vec, int size, int num) {
    for (int i = 0; i < size; ++i) {
        if (vec[i] == num) {
            return false; // Number found in the vector
        }
    }
    return true; // Number not found in the vector
}


int* build_cluster(double* spins, int* cluster, int nx, int ny, double beta, double J, int MAX_CLUSTER_SIZE, int starting){

	// Probability to add a spin to the cluster
	double p_add = 1.0 - exp(-2.0 * beta * J);
	
	// Get starting spin coordinates in the form (x,y)
	int x = starting % nx;
	int y = floor(starting / nx);
	
	// Define neighbors of the starting, using periodic boundary conditions		
	vector<int> neighbors = {  
        ((x - 1 + nx) % nx) + y * nx,  // Left
        ((x + 1) % nx) + y * nx,       // Right
        x + ((y + 1) % ny) * nx,       // Up
        x + ((y - 1 + ny) % ny) * nx,  // Down
    	};		
	

	// Get starting spin value
	int starting_spin = spins[starting];
	
	// Loop through neighbors & check if spin can be added
	double r = 0.;
	vector<int> new_elements;
	for (int neighbor : neighbors) {
		if (spins[neighbor] == starting_spin && notInsideVector(cluster, MAX_CLUSTER_SIZE, neighbor)) {
			r = pseudo_random();
		    	if (r < p_add) {
		        	new_elements.push_back(neighbor);
		        	//cout << "new element added --> " << neighbor << endl;
		    	}
		}
	}
	
	
	// find the index of the last non zero element of the cluster
	int last_non_zero_index = 0;
	for (int i = 0; i < MAX_CLUSTER_SIZE; ++i) {
		if (cluster[i] != 0) {
			last_non_zero_index = i;
	    	}
	}
	

	// Add elements to the cluster after the last non-zero element
	for (int new_element : new_elements) {
		cluster[last_non_zero_index + 1] = new_element; // Assign to cluster array after the last non-zero element
		last_non_zero_index++; // Increment index
	}
	
	return cluster;
	
}



tuple<double*, int> wolff(double* spins, int nx, int ny, double J, double beta, int MAX_CLUSTER_SIZE){
	
	// get random spin --> rnd_flip in [0,nx*ny) 
	mt19937 mt_generator(time(0)); 
    	uniform_int_distribution<int> uniform_flip(0, nx*ny); // rnd_flip = proposal spin fl
    	int rnd_flip = uniform_flip(mt_generator);
    	
	//compute row and column indexes from rnd_flip in [0,nx*ny)
	int i = rnd_flip % nx;       
	int j = floor(rnd_flip / nx);
	
	// The random spin is the first element of the cluster
	// Declare a 1D vector for the cluster indexes
   	vector<int> cluster(MAX_CLUSTER_SIZE);
   	cluster[0] = rnd_flip;
	
	// track cluster size
   	int current_cluster_size = 0;
   	
   	//build the cluster
   	for (int element : cluster){
   		if (element != 0) {

        		int* cluster_ptr = build_cluster(spins, cluster.data(), nx, ny, beta, J, MAX_CLUSTER_SIZE, element);
        		/*cout << "current clusters elements: " << endl;
        		for (int element : cluster){
        			if (element != 0){
        				cout << element << endl;
        			}
		   		
		   	}*/
        		
    		}
    
   	}
   	
   	//flip the cluster
   	/*for (int i=0; i<nx; i++){
   		for (int j=0; j<nx; j++){
   			cout << spins[j+i*nx] << "  " ;  	
   		}
   		//cout << endl;
   	}
   	
   	cout << "Flip the cluster ... " << endl;*/
   	
   	//flip the cluster and compute final cluster size
   	for (int element : cluster){
   		if (element != 0) {
   			current_cluster_size +=1 ;
   			spins[element] *=-1.;
   		}
   	}
   	
   	/*for (int i=0; i<nx; i++){
   		for (int j=0; j<nx; j++){
   			cout << spins[j+i*nx] << "  " ;  	
   		}
   		//cout << endl;
   	}*/
   	
	return make_tuple(spins, current_cluster_size);
}


// Function to initialize a vector to zero
void setVectorToZero(double* vector, int length) {
	fill(vector, vector + length, 0.0);
}

void setVectorToOne(double* vector, int length) {
	fill(vector, vector + length, 1.0);
}




int main(int argc, const char * argv[]) {

	int nx = 50;
	int ny = 50;
	long int nsteps = 1e4;
	int burnInSteps = 1e3; 
	double J = 1;
	double const T_critic = 2.0/log(1.0 + sqrt(2.0));
	double kbT = T_critic;
	double ratio = 0., r = 0.;
	int index = 0, iter_count = 0, realiz = 2;
	double diff_ene;
	int rnd_flip;
	double energy0, magnetization0, magn = 0;
	double E_ave = 0., M_ave = 0.;
	double beta = 1.0/kbT;
	int MAX_CLUSTER_SIZE = nx*ny;
	int current_cluster_size = 0;
	
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
    	double* cluster_sizes= new double[nsteps+1];
    	
    	//initialize energy and magnetization
        initializeSpins(nx, ny, spins, index);
    	energy0 = energia(J, nx, ny, spins, index);
    	magnetization0 = magnetization(nx, ny, spins, index);
        E[0] = energy0;
        M[0] = magnetization0;
        
	// Tests
	//vector<int> cluster(MAX_CLUSTER_SIZE);  // Pre-allocate space
	
	//int starting = uniform_flip(mt_generator);
	// int* cluster_ptr = build_cluster(spins, cluster.data(), nx, ny, beta, J, MAX_CLUSTER_SIZE, starting);


	// Run Wolff algorithm
	ofstream filemc;
	filemc.open("wolff_50_Tc2.txt",ios::out);
	filemc.precision(10);

	setVectorToZero(M, nx * ny);
	setVectorToZero(E, nx * ny);

 
	//setVectorToOne(spins, nx * ny); // random ---> initializeSpins(nx, ny, spins, index);
	initializeSpins(nx, ny, spins, index);
	
	// Equilibration
	for(long int t = 0; t < burnInSteps; t++){
		auto result =  wolff(spins, nx, ny, J, beta, MAX_CLUSTER_SIZE); //return tuple
    		double* spins = get<0>(result);
	}

	
	// Trajectory of monte-carlo steps
	for(long int t=0; t < nsteps; t++){
		cout << t << endl;
		// Keep trace of energy and magnetisation
    		M[t] = magnetization(nx, ny, spins, index);
    		E[t] = energia(J, nx, ny, spins, index);

    		// Evolve state at each mc step with 1 Wolf move
    		auto result =  wolff(spins, nx, ny, J, beta, MAX_CLUSTER_SIZE); //return tuple
    		double* spins = get<0>(result);
		int current_cluster_size = get<1>(result);
		cluster_sizes[t] = current_cluster_size;

	}	

	
	// compute the mean M, E 
	M_ave = 0.0;
	E_ave = 0.0;
	for(long int t=0; t < nsteps; t++){
		M_ave += M[t]/(nx*ny*nsteps);
		E_ave += E[t]/(nx*ny*nsteps);
		filemc << t << "\t" << E[t]/(nx*ny) << "\t" << M[t]/(nx*ny)  << "\t" << cluster_sizes[t] << '\n';
	}

    return 0;
}


