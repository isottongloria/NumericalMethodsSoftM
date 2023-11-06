/*
Sampling random points within D- dimensional domains by hit and miss
	1) Rectangle: Generate random points uniformly distributed within a rectangle [a, b] ⇥ [c, d] and compare the
	analytic value of the area A = LabLcd with the Monte Carlo estimate based on the hit-miss method as
	a function of the number of “throws”.
	
Disk Do the same for a unit radius disk.
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
	
        mt19937 mt_generator(time(0)); 
    	uniform_real_distribution<double> uniform(0.0, 1.0);
        
	//define the limits of the 2d bigger rectangle
	double A = 0;
	double B = 2;
	double C = 0;
	double D = 2;
	
	//define the limits of the 2d smaller rectangle
	double a = 0;
	double b = 1;
	double c = 0;
	double d = 1;
	
	// simulation data
	vector<int> numSamplesList =  { 1000, 5000, 10000, 100000, 500000, 1000000};
	int numSamples = 0.;
	int count_rect = 0;
	double X = 0., Y = 0.;
	
	ofstream outputFile1("2.1_rect_c.txt"); // Create an output file

	
	for (int n = 0; n < 6; n++) {
		numSamples = numSamplesList[n];
		count_rect = 0; // Reset the count
	   	for (int i = 0; i < numSamples; i++) {
			// Generate random coordinates within the larger rectangle
			double X = A + uniform(mt_generator) * (B - A);
			double Y = C + uniform(mt_generator) * (D - C);
			
			// Check if (X,Y) falls inside the smaller rectangle
			if (X >= a && X <= b && Y >= c && Y <= d) {
		        	count_rect ++;
			}
		}
		outputFile1 << numSamples << "\t" << (B - A)*(D - C)*count_rect/numSamples << "\n";
	}
	outputFile1.close();
	
	
	// DISK
	double R = 1., x = 0., y = 0.;
	int count_disk = 0;
	
	ofstream outputFile2("2.1_disk_c.txt"); // Create an output file
		
	for (int n = 0; n < 6; n++) {
		numSamples = numSamplesList[n];
		count_disk = 0; // Reset the count
		for (int i = 0; i < numSamples; i++) {
			// Generate random coordinates within a square in the first quadrant (0,1)*(0,1) 
			x = uniform(mt_generator);
			y = uniform(mt_generator);
			if (pow(x,2) + pow(y,2) <= pow(R,2) ) {
				count_disk ++;
			}
		}
		outputFile2<< numSamples << "\t" << 4.*count_disk/numSamples << endl;
	}
	outputFile2.close();
	
	
    return 0;
}
