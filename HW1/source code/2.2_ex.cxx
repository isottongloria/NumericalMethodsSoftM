#include <iostream>
#include <fstream>
#include <ctime>
#include<cmath>
#include <cstdlib>

using namespace std;

double Sampling_power_law (int n, double c){ 				// f(x) = c x^n,  x in [0,1]
        double u = static_cast<double>(rand()) / RAND_MAX;
        double x = static_cast<double>(pow(u*(n+1.0)/c,1.0/(n+1.0)));
	return x;
}

double Sampling_power_law_ext (double c){				// f(x) = c x^n,  x in [0,2]
        double u = static_cast<double>(rand()) / RAND_MAX;
        double x = static_cast<double>(pow(u*3.0/c,1.0/3.0));
	return x;
}

double Sampling_exp_1 (double mu){					// f(x) = mu* e^(-mu*x), x in [0,Inf]
        double u = static_cast<double>(rand()) / RAND_MAX;
        double x = static_cast<double>(log(1.0-u)/(-mu));
	return x;
}

double Sampling_exp_2 (){						// f(x) = 2*x e^(x**2), x in [0,Inf]
        double u = static_cast<double>(rand()) / RAND_MAX;
        double x = static_cast<double>(sqrt(log(1.0/(1.0-u))));
	return x;
}

int main() {
	srand(time(0));
	int n = 4; //power
	double c = 4.0; //normalization
	double X = 0; //sampled values
	int numSamples = 1000;
	
	ofstream outputFile1("sampl_power_3.txt"); // output file
	for (int i = 0; i < numSamples; i++) {
		X = static_cast<double>(Sampling_power_law(n,c));
		outputFile1 << X << "\n";
	}
	
	int m = 5;  //power
	ofstream outputFile2("sampl_power_4.txt"); // output file
	for (int i = 0; i < numSamples; i++) {
		X = static_cast<double>(Sampling_power_law(m,c));
		outputFile2 << X << "\n";
	}
	
	double c_ext = 3./8.; //norm constant
	ofstream outputFile3("sampl_power_ext.txt"); // output file
	for (int i = 0; i < numSamples; i++) {
		X = static_cast<double>(Sampling_power_law_ext(c_ext));
		outputFile3 << X << "\n";
	}
	
	double mu = 1.;
	ofstream outputFile4("sampl_exp.txt"); // output file
	for (int i = 0; i < numSamples; i++) {
		X = static_cast<double>(Sampling_exp_1(mu));
		outputFile4 << X << "\n";
	}
	

	ofstream outputFile5("sampl_exp2.txt"); // output file
	for (int i = 0; i < numSamples; i++) {
		X = static_cast<double>(Sampling_exp_2());
		outputFile5 << X << "\n";
	}
        return 0;
}
