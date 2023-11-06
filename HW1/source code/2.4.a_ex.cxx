/*
* compute an integral in 2 ways :
* 1 - crude montecarlo
* 2 - importance sampling 
*/

#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <random>


using namespace std;

double f(double x){
	return exp(-pow(x,2))*log(x);
}

double g(double x){
	return log(x);
}


int main() {

// 1 - crude montecarlo

	//Random number generator (uniform)
	default_random_engine generator;
        uniform_real_distribution<double> uniform(0.0,1.0);
        
        //parameters for the simulation, [a,b] = extremes of integration
        double a = 0.0;
        double b = 5000;
  	int numSamples = 1e7;
  	double integral_MC = 0.0;
  	double integral_IS = 0.0;
  	double x = 0.0;
  	
  	//crude monte-carlo
	for (int i = 0; i < numSamples; i++) {
		x = a + (b - a) * uniform(generator);
		integral_MC += (b - a) * f(x)/numSamples;
	}
	cout << "The analytical integral estimation is:" <<-0.8700577267283155 << endl;
	cout << "The integral estimation with crude Monte Carlo is: " << integral_MC << endl;

// 2 - importance sampling 
	normal_distribution<double> normal(0.0,1.0/sqrt(2.0)); 
	
	for (int i = 0; i < numSamples; i++) {
		x = sqrt(2.0) * abs(normal(generator)); //sqrt(2.0) *
		integral_IS += sqrt(M_PI/2.0) * g(x)/numSamples;
	}
	cout << "The integral estimation with Importance Sampling is: " << integral_IS << endl;
 	return 0;
}
