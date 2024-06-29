/*
* compute an integral with importance sampling 
*/

#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <random>
#include <utility>

using namespace std;

double f(double x){  //target
	return cos(x);
}

double g(double x, double A, double B, double Z){  //candidate
	return (A + B * x * x)/Z ;
}

double sample_g (double u, double A, double B){ // extremely complicated formula from inversion method
	return -(24 * pow(2.0, 1.0 / 3) * A) / pow(20736 * A * B * B * M_PI * u + 1728 * B * B * B * pow(M_PI, 3) * u + sqrt(764411904 * pow(A, 3) * pow(B, 3) + pow(20736 * A * B * B * M_PI * u + 1728 * B * B * B * pow(M_PI, 3) * u, 2)),(1.0 / 3)) + (1 / (24 * pow(2.0, 1.0 / 3) * B)) * (pow(20736 * A * B * B * M_PI * u + 1728 * B * B * B * pow(M_PI, 3) * u + sqrt(764411904 * pow(A, 3) * pow(B, 3) + pow(20736 * A * B * B * M_PI * u + 1728 * B * B * B * pow(M_PI, 3) * u, 2)), 1.0 / 3));

}


// Calculate the standard deviation and the mean of sampled points
pair<double, double> Mean_Std(const vector<double>& estimates) {
    double sampleMean = 0.0;
    for (int i = 0; i < estimates.size(); i++) {
        sampleMean += estimates[i];
    }
    	
    sampleMean /= estimates.size();

    double sampleVariance = 0.0;
    for (int i = 0; i < estimates.size(); i++) {
        sampleVariance += (estimates[i] - sampleMean) * (estimates[i] - sampleMean);
    }
    sampleVariance /= (estimates.size() - 1);

    return make_pair(sampleMean, sqrt(sampleVariance));
}


int main() {

// 1 - crude montecarlo

	//Random number generator (uniform)
	default_random_engine generator;
        uniform_real_distribution<double> uniform(0.0,1.0);
        
        //parameters for the simulation, [a,b] = extremes of integration
        double a = 0.0;
        double b = M_PI / 2.0;
  	int numSamples = 1e6;
  	double integral_MC = 0.0;
  	double integral_IS = 0.0;
  	double x = 0.0; 
  	
  	// parameters for function g(x) = A + B x^2 and its normalization Z
  	double A = 0.5;
  	double B = 1.0;
  	double Z = (A * M_PI)/2 + B * pow(M_PI,3)/24;
  	
  	//vector for storing sampled values from importance sampling
  	vector<double> estimates;
  	
  	//crude monte-carlo
	for (int i = 0; i < numSamples; i++) {
		x = a + (b - a) * uniform(generator);
		integral_MC += (b - a) * f(x)/numSamples;
	}
	cout << "The integral estimation with crude Monte Carlo is: " << integral_MC << endl;


// 2 - importance sampling 

	ofstream outputFile("mean_std_N.txt");
	numSamples = 1e1;
	for(int n = 0; n <= 6; n++){
	numSamples = numSamples*10;
	integral_IS = 0;
		for (int i = 0; i < numSamples; i++) {
			x = sample_g (uniform(generator), A, B); 
			double estimate = (f(x) / g(x, A, B, Z));
			estimates.push_back(estimate);
			integral_IS += (f(x)/g(x,A,B,Z))/numSamples;
		}
		
	auto result = Mean_Std(estimates);
	double standardDeviation = result.second;
	double sampleMean = result.first;

	cout << "done with " <<numSamples << " samples" << endl;
	
	//save to file
	outputFile <<  static_cast<double>(numSamples) << "\t" << integral_IS << "\t" << standardDeviation << "\n";
	}

 	return 0;
}
