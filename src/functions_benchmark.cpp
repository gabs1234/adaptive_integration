#include "functions.hpp"
#include "utilities.hpp"
#include <cstdlib>
#include <chrono>
#include <limits>

using namespace std;

template <class T>
T f(T x){
	return pow(sin(x), 2);
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

	// float a = strtof(argv[1], NULL);
	// float b = strtof(argv[2], NULL);
	float a = -M_PI;
	float b = M_PI;

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

	legendre_roots<long double>(n, x, eps);
	gauss_legendre_weights<long double>(n, x, w);

	// Solve
	auto tic1 = chrono::steady_clock::now();
	long double sol = gauss_legendre<long double>(a, b, f, n, x, w, eps);
	auto tac1 = chrono::steady_clock::now();

	cout << n << "\t";
	cout << chrono::duration<double, milli>(tac1 - tic1).count() << "\t";
	cout << sol << "\t";
	cout << (M_PI - sol)/sol << endl;

	return 0;
}
