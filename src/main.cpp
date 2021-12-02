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

	cout.precision(nb_digits);

	// Calculate the roots
	auto tic1 = chrono::steady_clock::now();
  legendre_roots(n, x, eps);
	auto tac1 = chrono::steady_clock::now();
	// Calculate the weights
	auto tic2 = chrono::steady_clock::now();
  gauss_legendre_weights(n, x, w);
	auto tac2 = chrono::steady_clock::now();

	cout << chrono::duration<double, milli>(tac1 - tic1).count() << "\t";
	cout << chrono::duration<double, milli>(tac2 - tic2).count() << endl;


	// // use adaptive integration
	// auto tic1 = chrono::steady_clock::now();
  // long double sol_adapt = adaptive_integration<long double>(a, b, f, n, eps);
	// auto tac1 = chrono::steady_clock::now();
	//
	//
	// long double exact_solution = 2*( atan(b) - atan(a) );
	//
	// cout << sol_adapt << "\t";
	// cout << chrono::duration<double, milli>(tac1 - tic1).count() << "\t";
	// cout << exact_solution << endl;


	return 0;
}
