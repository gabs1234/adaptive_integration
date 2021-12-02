#include "functions.hpp"
#include "utilities.hpp"
#include <cstdlib>
#include <chrono>
#include <limits>

using namespace std;

template <class T>
T f(T x){
	return 2.0/(1+x*x);
}

int main(int argc, char *argv[]) {

	if( argc == 1){
		cout << "usage: ./main <a> <b> <prec> <n>" << endl;
		cout << "[a, b]: integration interval" << endl;
		cout << "prec: precision of floating point" << endl;
		cout << "n: order of gauss-legendre method" << endl;
		return 1;
	}

	int max_prec = numeric_limits<long double>::digits10;

	float a = strtof(argv[1], NULL);
  float b = strtof(argv[2], NULL);
	int nb_digits = strtold(argv[3], NULL);
	if( nb_digits > max_prec ){
		cout << "maximal machine precision of: " << max_prec;
		nb_digits =  max_prec;
	}
	long double eps = pow(10, -nb_digits);
	int n = atoi(argv[4]);

	// Find roots and weights
	long double x[n], w[n];
	long double xp = .5;

	cout.precision(nb_digits);

	// Calculate polynomial in x
	auto tic1 = chrono::steady_clock::now();
  P(n, xp);
	auto tac1 = chrono::steady_clock::now();

	// Calculate derivative of P in x
	auto tic2 = chrono::steady_clock::now();
  dP(n, xp);
	auto tac2 = chrono::steady_clock::now();

	// Calculate the roots
	auto tic3 = chrono::steady_clock::now();
  legendre_roots(n, x, eps);
	auto tac3 = chrono::steady_clock::now();

	// Calculate the weights
	auto tic4 = chrono::steady_clock::now();
  gauss_legendre_weights(n, x, w);
	auto tac4 = chrono::steady_clock::now();

	// Solve

	cout << n << "\t";
	cout << chrono::duration<double, milli>(tac1 - tic1).count() << "\t";
	cout << chrono::duration<double, milli>(tac2 - tic2).count() << "\t";
	cout << chrono::duration<double, milli>(tac3 - tic3).count() << "\t";
	cout << chrono::duration<double, milli>(tac4 - tic4).count() << endl;

	return 0;
}
