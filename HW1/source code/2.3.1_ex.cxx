#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <random>

using namespace std;

double f(double x, const double N ){
	return N * exp(-pow(x,2));
}

double g(double x, double p, double A, double sigma) {
    if (x <= p) {
        return A;
    } else {
        return A * exp(-(pow(x - p,2)) / pow(2 * sigma, 2));
    }
}

int main() {

	double p = 0.85;
	double A = 1.0/p;
	double c = A / 0.8;
	double sigma = 1.0 ;
	const double N = sqrt(2.0 / M_PI);
	
	default_random_engine generator; //random number generator
        uniform_real_distribution<double> uniform_1(0.0,1.0);
        uniform_real_distribution<double> uniform_p(0.0,p); 
        normal_distribution<double> distribution(p,sigma);
        
        //rejection method on 0 <= x <= p
        
  	int numSamples = 1e5;  // number of samplings
	double ratio = 0;
	ofstream outputFile1("sampl_rejection_to_p.txt"); // output file
	for (int i = 0; i < numSamples; i++) {
	
		//Generate X according to g(x)=A in [0,p] and a uniform random number
		double u = static_cast<double>(uniform_1(generator)); 
		double X = static_cast<double>(uniform_p(generator)); 
		
		ratio = f(X, N)/(c*g(X, p, A, sigma));
		if (u < ratio){
		outputFile1 << X << "\n";
		}
	}
	
	//rejection method on x > p
	ofstream outputFile2("sampl_rejection_from_p.txt"); // output file
	int count = 0;
	while (count < numSamples) {
	
		//Generate X according to g(x) in [p,inf] and a uniform random number
		double u = static_cast<double>(uniform_1(generator)); 
		double X = abs(static_cast<double>(distribution(generator)));
		
		if (X > p){
			count ++;
			ratio = f(X, N)/(c*g(X, p, A, sigma));
			if (u < ratio){
				outputFile2 << X << "\n";
			}
		}
	}
 	return 0;
}
