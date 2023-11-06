/*
* 1 ) Sampling uniformly points within a unit radius disk
* r = u1 and theta = 2 Pi u2 with u1, u2 uniformly distributed in [0,1]
*
* 2 ) Sampling from a gaussian 2d (mean 0, Var 1)
*/

#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <random>

using namespace std;


int main() {

	/* 
	*  Sampling unit radius disk naive
	*/
	default_random_engine generator;
        uniform_real_distribution<double> distribution(0.0,1.0);
  
	int numSamples = 1e4;  //set number of samplings
	vector<double> r;      //define a vector for storing radius 
	vector<double> theta;
		
	ofstream outputFile("sampl_disk.txt"); // output file
	for (int i = 0; i < numSamples; i++) {
	
		//sample two rnd variables uniformly in [0,1]
		double u1 = static_cast<double>(distribution(generator)); 
		double u2 = static_cast<double>(distribution(generator));
		
		//sample within the unit radius disk
		r.push_back(u1);
		theta.push_back(2*M_PI*u2);
		outputFile << r.back()*cos(theta.back()) << "\t" << r.back()*sin(theta.back()) << "\n";
	}

	/*
	*  Sampling unit radius disk correct, Box muller
	*/
	double u1, u2 = 0;
	vector<double> X;      //define a vector for storing X, Y 
	vector<double> Y;
	
	ofstream outputFile_a("sampl_disk_corr.txt"); // output file

	for (int i = 0; i < numSamples; i++) {
	    u1 = static_cast<double>(distribution(generator));
	    u2 = static_cast<double>(distribution(generator));
	    X.push_back(sqrt(u1) * cos(2 * M_PI * u2));
	    Y.push_back(sqrt(u1) * sin(2 * M_PI * u2));
	    outputFile_a << X[i] << "\t" << Y[i] << "\n";
	}
	
	
	/* 
	* Sampling from a gaussian 2d
	*/
	vector<double> x;      //define a vector for storing X, Y 
	vector<double> y;
	
	ofstream outputFileX("sampl_gauss_X.txt"); // output file
	ofstream outputFileY("sampl_gauss_Y.txt"); // output file
	for (int i = 0; i < numSamples; i++) {
	    u1 = static_cast<double>(distribution(generator));
	    u2 = static_cast<double>(distribution(generator));
	    
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
