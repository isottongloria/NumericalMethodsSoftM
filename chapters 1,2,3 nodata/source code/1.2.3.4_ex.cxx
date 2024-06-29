#include <iostream>
#include <fstream>
#include <cmath>
#include <random>

using namespace std;

// Exercise 1.2

double Sampling_power_law(int n, double c, mt19937& gen) { 		// f(x) = c x^n,  x in [0,1]
    uniform_real_distribution<double> dist(0.0, 1.0);
    double u = dist(gen);
    double x = pow(u * (n + 1.0) / c, 1.0 / (n + 1.0));
    return x;
}

// Exercise 1.3

double Sampling_power_law_ext(double c, mt19937& gen) {			// f(x) = c x^n,  x in [0,2]
    uniform_real_distribution<double> dist(0.0, 1.0);
    double u = dist(gen);
    double x = pow(u * 3.0 / c, 1.0 / 3.0);
    return x;
}

// Additional exercises 1.4.1

double Sampling_exp_1(double mu, mt19937& gen) {			// f(x) = mu* e^(-mu*x), x in [0,Inf]
    uniform_real_distribution<double> dist(0.0, 1.0);
    double u = dist(gen);
    double x = log(1.0 - u) / (-mu);
    return x;
}

// Additional exercises 1.4.2

double Sampling_exp_2(mt19937& gen) {					// f(x) = 2*x e^(x**2), x in [0,Inf]
    uniform_real_distribution<double> dist(0.0, 1.0);
    double u = dist(gen);
    double x = sqrt(log(1.0 / (1.0 - u)));
    return x;
}

// Additional exercises 1.4.3

double Sampling_exp_3(mt19937& gen) {	               			// f(x) = ab/(a+bx)^n, x in [0,Inf]
    uniform_real_distribution<double> dist(0.0, 1.0);
    double u = dist(gen);
    double x = (2.0) * ( pow((1.0 - u), - 0.5) - 1.0 );
    return x;
}


int main() {
    mt19937 gen(time(0)); // Mersenne Twister PRNG
    int n = 4; // power
    double c = 4.0; // normalization
    double X = 0; // sampled values
    int numSamples = 1e6;

    ofstream outputFile1("sampl_1.2.txt"); // output file
    for (int i = 0; i < numSamples; i++) {
        X = Sampling_power_law(n, c, gen);
        outputFile1 << X << "\n";
    }

    int m = 5;  // power
    ofstream outputFile2("sampl_1.2b.txt"); // output file
    for (int i = 0; i < numSamples; i++) {
        X = Sampling_power_law(m, c, gen);
        outputFile2 << X << "\n";
    }

    double c_ext = 3. / 8.; // norm constant
    ofstream outputFile3("sampl_1.3.txt"); // output file
    for (int i = 0; i < numSamples; i++) {
        X = Sampling_power_law_ext(c_ext, gen);
        outputFile3 << X << "\n";
    }

    double mu = 1.;
    ofstream outputFile4("sampl_1.4.1.txt"); // output file
    for (int i = 0; i < numSamples; i++) {
        X = Sampling_exp_1(mu, gen);
        outputFile4 << X << "\n";
    }

    ofstream outputFile5("sampl_1.4.2.txt"); // output file
    for (int i = 0; i < numSamples; i++) {
        X = Sampling_exp_2(gen);
        outputFile5 << X << "\n";
    }
    

    n = 3;
    numSamples = 1e7;
    ofstream outputFile6("sampl_1.4.3.txt"); // output file
    for (int i = 0; i < numSamples; i++) {
        X = Sampling_exp_3(gen);
        outputFile6 << X << "\n";
    }
    
    return 0;
}

