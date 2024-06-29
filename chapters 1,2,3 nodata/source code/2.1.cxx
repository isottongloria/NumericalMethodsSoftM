/*
*  1 ) Sampling uniformly points within a unit radius disk
*  r = u1 and theta = 2 Pi u2 with u1, u2 uniformly distributed in [0,1]
*
*  2 ) Sampling from a gaussian 2d (mean 0, Var 1)
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

	
	// Sampling unit radius disk naive
	
	mt19937 generator(time(0)); // Mersenne Twister random number generator
        uniform_real_distribution<double> distribution(0.0, 1.0);
  
	int numSamples = 1e4;  //set number of samplings
	vector<double> r;      //define a vector for storing radius 
	vector<double> theta;
	double u1, u2 = 0;
	
	ofstream outputFile("disk_sampl_2.1.txt"); // output file
	for (int i = 0; i < numSamples; i++) {
	
		//sample two rnd variables uniformly in [0,1]
		u1 = distribution(generator); 
		u2 = distribution(generator);
		
		//sample within the unit radius disk
		r.push_back(u1);
		theta.push_back(2*M_PI*u2);
		outputFile << r.back()*cos(theta.back()) << "\t" << r.back()*sin(theta.back()) << "\n";
	}

	
	
	//  Sampling unit radius disk correct, Box muller
	
	vector<double> X;      //define a vector for storing X, Y 
	vector<double> Y;
	
	ofstream outputFile_a("disk_correct_2.1.txt"); // output file

	for (int i = 0; i < numSamples; i++) {
	        //sample two rnd variables uniformly in [0,1]
		u1 = distribution(generator); 
		u2 = distribution(generator);
 		
 		// store samples
		X.push_back(sqrt(u1) * cos(2 * M_PI * u2));
		Y.push_back(sqrt(u1) * sin(2 * M_PI * u2));
		outputFile_a << X[i] << "\t" << Y[i] << "\n";
	}
	
 	return 0;
}
