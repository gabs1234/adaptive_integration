#include "functions.hpp"
#include "utilities.hpp"
#include <cstdlib>
#include <chrono>

using namespace std;

template <class T>
T f(T x){
	return std::exp(-std::pow(x, 2.0));
}

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

	// use adaptive integration
	auto tic1 = chrono::steady_clock::now();
  long double sol_adapt = adaptive_integration<float>(a, b, f, n, eps);
	auto tac1 = chrono::steady_clock::now();

	cout.precision(12);

	long double exact_solution = sqrt(M_PI);

	cout << sol_adapt << "\t";
	cout << chrono::duration<double, milli>(tac1 - tic1).count() << "\t";
	cout << exact_solution << endl;


	return 0;
}
