#include "headers.h"
#include <chrono>
#include <cstdlib>

using namespace std;

long double f(long double x);

int main(int argc, char *argv[]) {
	if( argc == 1){
		cout << "usage: ./main <a> <b> <prec> <n>" << endl;
		cout << "[a, b]: integration interval" << endl;
		cout << "prec: precision of floating point" << endl;
		cout << "n: order of gauss-legendre method" << endl;
		return 1;
	}

	float a = strtof(argv[1], NULL);
  float b = strtof(argv[2], NULL);
	long double nb_digits = strtold(argv[3], NULL);
  long double eps = pow(10, -nb_digits);
  int n = atoi(argv[4]);

	auto tic1 = chrono::steady_clock::now();
  long double sol = gauss_legendre<long double>(a, b, f, n, eps);
	auto tac1 = chrono::steady_clock::now();

	auto tic2 = chrono::steady_clock::now();
  long double sol_adapt = adaptive_integration<long double>(a, b, f, n, eps, 0);
	auto tac2 = chrono::steady_clock::now();

	cout.precision(12);

	// cout << "sol:\t" << "time1:\t"  << "sol_adapt:\t" << "time2 :"  << endl;
	cout << sol << "\t";
	cout << chrono::duration<double, milli>(tac1 - tic1).count() << "\t";
	cout << sol_adapt << "\t";
	cout << chrono::duration<double, milli>(tac2 - tic2).count() << endl;
  // printf("%lf\n", sol);
  // // printf("%lf\n", sol_adapt);

	return 0;
}

long double f(long double x){
	return pow(x, 4) * exp(-x);
}
