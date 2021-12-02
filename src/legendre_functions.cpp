#include "functions.hpp"
#include <iostream>

using namespace std;


int main(int argc, char const *argv[]) {
	if( argc != 3){
		cout << "usage: ./main <n> <nb_points>" << endl;
		cout << "n: order of gauss-legendre method" << endl;
		cout << "nb_points: number of points in which the legendre function is evaluated" << endl;
		return 1;
	}

  int n = atoi(argv[1]);
  int nb_points = atoi(argv[2]);

	long double x, eps = 10e-6;
	long double roots[n];

	for (int i = 0; i < nb_points; i++) {
		x = -1 + 2.0*i/(nb_points-1);
		cout << x << "\t" << P<long double>(n, x) << endl;
	}
	

	return 0;
}
