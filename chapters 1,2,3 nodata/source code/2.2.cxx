#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <random>

using namespace std;


int main() {

	// Sampling from a gaussian 2d
	
	mt19937 generator(time(0)); // Mersenne Twister random number generator
        uniform_real_distribution<double> distribution(0.0, 1.0);
	int numSamples = 1e4;  // set number of samplings
	double u1, u2 = 0;
	
	//define a vector for storing X, Y 
	vector<double> x;      
	vector<double> y;
	
	ofstream outputFileX("sampl_gauss_X.txt"); // output file
	ofstream outputFileY("sampl_gauss_Y.txt"); // output file
	for (int i = 0; i < numSamples; i++) {
	    u1 = distribution(generator);
	    u2 = distribution(generator);
	    
	    if (u1 == 0.0 || u2 == 0.0) {
		// Handle the case where u1 or u2 is zero to avoid NaN
		i--; // Decrement i and try again
	    } else {
		x.push_back(sqrt(-2 * log(u1)) * cos(2 * M_PI * u2));
		y.push_back(sqrt(-2 * log(u1)) * sin(2 * M_PI * u2));
		outputFileX << x[i] << "\n";
		outputFileY << y[i] << "\n";
	    }
	}
	
 	return 0;
}
