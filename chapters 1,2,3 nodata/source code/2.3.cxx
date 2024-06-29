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

double g(double x, double p, double A) {
    if (x <= p) {
        return A;
    } else {
        return A/p *x* exp(p*p-x*x);
    }
}

int main() {

	double p = 0.9;
	double A = (2.*p)/(2*p*p + 1);
	double c = sqrt(2/M_PI)/A;
	double t = A*p;
	const double N = sqrt(2.0 / M_PI);
	
	mt19937 generator(time(0)); // Mersenne Twister rnd number generator (seed current time)
	uniform_real_distribution<double> uniform_1(0.0, 1.0);

        
  	int numSamples = 1e6;  // number of samplings
	double ratio, u, v, X, sample = 0;
	ofstream outputFile1("sampl_rejection_p09.txt"); // output file
	for (int i = 0; i < numSamples; i++) {
	
		//1) Generate X according to g(x) 
		v = uniform_1(generator);
		if (v <= t){ 
			X = v/A;//uniform_p(generator); 
			//cout << X << endl;
		}
		
		else
		{
			X = sqrt(p*p - log(2*p*(1-u)/A));
			//cout << X << endl;
		}
		
		//2) Generate a uniform random number and compute the ratio
		u = uniform_1(generator);
		ratio = f(X, N)/(c*g(X, p, A));
		if (u < ratio){
			// accept the sampling and store it
			outputFile1 << X << "\t" << c*g(X, p, A)*u << "\n";
		}
	}
	
 	return 0;
}
