/*
This C++ program employs the hit-and-miss method to sample random points within D-dimensional domains. It examines two shapes: a rectangle and a unit radius disk. For the rectangle, random points are generated within given bounds [a, b] × [c, d], comparing Monte Carlo area estimates with analytical values. The same comparison is conducted for the unit radius disk.
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

	// RECTANGLE

	// Random number generator
    	mt19937 mt_generator(time(0));
    	uniform_real_distribution<double> uniform(0.0, 1.0);

    	// Define the limits of the 2D bigger rectangle where to sample unif.
	double A = 0;
	double B = 2;
	double C = 0;
	double D = 2;

	// Define the limits of the 2D smaller rectangle
	double a = 0;
	double b = 1;
	double c = 0;
	double d = 1;
	 
	// Simulation data
	vector<int> numSamplesList = {1000, 5000, 10000, 100000, 500000, 1000000}; // Possible throws
	int numSamples = 0.;
	int numSimulations = 50; // Number of simulations to run for each numSamplesList
	int count_rect = 0;
	double X = 0., Y = 0.;
	vector<double> results; // Contain the 50 estimates of area of rectangle via hit miss method 

	ofstream outputFile1("1.1_rect.txt"); // Create an output file
	outputFile1 << "#Area of the rectangle" <<"\n";
	outputFile1 << "#nsamples"  << "\t" <<  "mean"  << "\t" << "std"  << "\n";
	
	for (int n = 0; n < 6; n++) { // Run over numSamplesList
		numSamples = numSamplesList[n];
		results.clear(); // Clear previous results

		for (int s = 0; s < numSimulations; s++) { // Run over realizations
		    count_rect = 0; // Reset the count

		    for (int i = 0; i < numSamples; i++) {
			// Generate random coordinates within the larger rectangle
			double X = A + uniform(mt_generator) * (B - A);
			double Y = C + uniform(mt_generator) * (D - C);

			// Check if (X,Y) falls inside the smaller rectangle
			if (X >= a && X <= b && Y >= c && Y <= d) {
			    count_rect++;
			}
		    }
		    results.push_back((B - A) * (D - C) * count_rect / numSamples);
		}

		// Calculate mean and standard deviation
		double mean = accumulate(results.begin(), results.end(), 0.0) / results.size();
		double var = 0.0;
		for (double value : results) {
		    var += (value - mean) * (value - mean);
		}
		var /= results.size();
		double stdev = sqrt(var);
		outputFile1 << numSamples << "\t" << mean << "\t" << stdev << "\n";
	}
	
	outputFile1.close();


	// DISK
	double R = 1., x = 0., y = 0.;
	int count_disk = 0;
	
	ofstream outputFile2("1.1_disk.txt"); // Create an output file
	outputFile2 << "#Area of the disk of radius R = " << R <<"\n";
	outputFile2 << "#nsamples"  << "\t" <<  "mean"  << "\t" << "std"  << "\n";
	
	for (int n = 0; n < 6; n++) {
		numSamples = numSamplesList[n];
		results.clear(); // Clear previous results

		for (int s = 0; s < numSimulations; s++) {
		    count_disk = 0; // Reset the count

		    for (int i = 0; i < numSamples; i++) {
			// generate random coordinates within a square in the first quadrant (0,R)*(0,R) 
			x = uniform(mt_generator)*R;
			y = uniform(mt_generator)*R;
			
			// Check if (x,y) falls inside the ¼disk
			if (pow(x,2) + pow(y,2) <= pow(R,2)) {
			    count_disk++;
			}
		    }
		    results.push_back(4.*count_disk/numSamples);
		}

		// Calculate mean and standard deviation
		double mean = accumulate(results.begin(), results.end(), 0.0) / results.size();
		double var = 0.0;
		for (double value : results) {
		    var += (value - mean) * (value - mean);
		}
		var /= results.size();
		double stdev = sqrt(var);
		outputFile2 << numSamples << "\t" << mean << "\t" << stdev << "\n";
	}
	
	outputFile2.close();
	
	
    return 0;
}
