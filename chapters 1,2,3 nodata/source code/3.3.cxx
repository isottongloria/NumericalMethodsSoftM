#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <random>

using namespace std;

double f(double x, double T){  //target
	if (x <= T) {
		return 0;
	    } else {
		return exp(- x);
	    }
}

double sample_g(double u, double a){  //candidate
	return log(1 / (1- u)) / a ;
}

double g(double x, double a){
	return a * exp(- a * x);
}

double rho(double x){
	return exp(-x);
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

	//Random number generator (uniform)
	default_random_engine generator;
        uniform_real_distribution<double> uniform(0.0,1.0);
        
        // Parameters for the simulation
        double a = 0.5; // parameter of g(x)
        double estimate = 0.0;
        double T = 9.0; // parameter of f(x)
        double numSamples = 1e1;
        double x = 0.0;
        double integral_IS = 0.0;
        double integral_MC = 0.0;
  	double A = 0.0; //extremes of integration with crude mc
        double B = 100.0; //extremes of integration with crude mc
        
// 1 - crude montecarlo integration

  	//vector for storing sampled values from monte carlo
  	/*vector<double> estimates_mc;
  	
        // Output file
	ofstream outputFile_mc("A_4_MC_mean_std_N.txt");
	
	// Simulation
	for(int n = 0; n <= 6; n++){
	numSamples = numSamples*10;
	integral_MC = 0;
	
		for (int i = 0; i < numSamples; i++) {
			x = A + (B - A) * uniform(generator);
			integral_MC += (B - A) * f(x, T) * rho(x)/numSamples;
		}
	outputFile_mc <<  static_cast<double>(numSamples) << "\t" << integral_MC << "\n";
}*/

// Search best value of a that minimize the variance

	// Output file
	
	ofstream outputFile_a("A_4_std_a.txt");
	
	// Simulation
	//vector for storing sampled values from importance sampling
  	vector<double> estimates;
	numSamples = 1e7;
	a = 0.0;
	for(int n = 0; n <= 9; n++){
	a = a + 0.1;
	integral_IS = 0;
		for (int i = 0; i < numSamples; i++) {
			x = sample_g (uniform(generator), a);
			estimate = (f(x, T) * rho(x) / g(x, a));  
			estimates.push_back(estimate);
			integral_IS += (f(x, T)* rho(x)/g(x, a))/numSamples;
		}
		
	auto result = Mean_Std(estimates);
	double standardDeviation = result.second;
	double sampleMean = result.first;

	cout << "done with " <<numSamples << " samples" << endl;
	
	//save to file
	outputFile_a << a << "\t" << integral_IS << "\t" << standardDeviation << "\n";
	}
	
// Search best value of T  that minimize the variance	

	// Output file
	ofstream outputFile_b("A_4_std_T.txt");
	
	// Simulation
	numSamples = 1e7;
	a = 0.1;
	double T_list[4] = {3,5,10,20};
	
	for(int n = 0; n <= 3; n++){
	T = T_list[n];
	integral_IS = 0;
		for (int i = 0; i < numSamples; i++) {
			x = sample_g (uniform(generator), a);
			estimate = (f(x, T) * rho(x)/ g(x, a));
			estimates.push_back(estimate);
			integral_IS += (f(x, T) * rho(x)/g(x, a))/numSamples;
		}
		
	auto result = Mean_Std(estimates);
	double standardDeviation = result.second;
	double sampleMean = result.first;

	cout << "done with " <<numSamples << " samples" << endl;
	
	//save to file
	outputFile_b << T << "\t" << integral_IS << "\t" << standardDeviation << "\n";
	}	
	
	
/*
//  Importance sampling 
        
        // Best a and T
        a = 0.18;
        T = 20;
        //vector for storing sampled values from importance sampling
  	vector<double> estimates;
  	
        // Output file
	ofstream outputFile("A_4_mean_std_N.txt");
	
	// Simulation
	numSamples = 1e1;
	for(int n = 0; n <= 6; n++){
	numSamples = numSamples*10;
	integral_IS = 0;
		for (int i = 0; i < numSamples; i++) {
			x = sample_g (uniform(generator), a);
			estimate = (f(x, T) / g(x, a));
			estimates.push_back(estimate);
			integral_IS += (f(x, T)/g(x, a))/numSamples;
		}
		
	auto result = Mean_Std(estimates);
	double standardDeviation = result.second;
	double sampleMean = result.first;

	cout << "done with " <<numSamples << " samples" << endl;
	
	//save to file
	outputFile <<  static_cast<double>(numSamples) << "\t" << sampleMean << "\t" << standardDeviation << "\n";
	}*/

 	return 0;
}
