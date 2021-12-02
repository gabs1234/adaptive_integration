#include "para_functions.hpp"
#include "utilities.hpp"
#include "mpi.h"
#include <chrono>
#include <limits>

using namespace std;

template <class T>
T f(T x){
	return sin(x)/x;
}

template <class T>
T f_sol(T x){
	return (1/4)*x*x*x*x+(1/3)*x*x*x+.5*x*x+x;
}

int main(int argc, char *argv[]) {
	// Test valid input
	if( argc == 1){
		cout << "usage: ./main <a> <b> <prec> <n>" << endl;
		cout << "[a, b]: integration interval" << endl;
		cout << "prec: precision of floating point" << endl;
		cout << "n: order of gauss-legendre method" << endl;
		return 1;
	}

	// init integration variables
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

	cout.precision(nb_digits);

	// init MPI
	int rank, nb_procs;
  MPI_Status status;
  MPI::Init();
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&nb_procs);

	float start=0, end=0;

	double t1 = MPI_Wtime();

	// Initialise starting intervals for integration
	if( rank == 0 ){
		cout << "The first rank" << endl;
		// Load balancing: find subintervals for each process
		float interval[nb_procs+1];
		load_balance<float>(a, b, nb_procs, f, interval);

		print_list<float>(nb_procs+1, interval);

		for( int i = 1; i < nb_procs; i++ ){
			start = interval[i];
			end = interval[i+1];
			MPI_Send(&start, 1, MPI_FLOAT, i, i, MPI_COMM_WORLD);
			MPI_Send(&end, 1, MPI_FLOAT, i, i, MPI_COMM_WORLD);
		}
		start = interval[0];
		end = interval[1];

	}
	else{
		MPI_Recv(&start, 1, MPI_FLOAT, 0, rank, MPI_COMM_WORLD, NULL);
		MPI_Recv(&end, 1, MPI_FLOAT, 0, rank, MPI_COMM_WORLD, NULL);
	}

	// use adaptive integration per process
  long double sol_adapt = adaptive_integration<long double>(start, end, f, n, eps);

	double t2 = MPI_Wtime();
	cout << rank << ": " << sol_adapt << "\t" << t2 - t1 << endl;

	if( rank != 0 ){
		MPI_Send(&sol_adapt, 1, MPI_LONG_DOUBLE, 0, rank, MPI_COMM_WORLD);
	}
	else{
		long double sol;
		long double final_sol = sol_adapt;
		for( int i = 1; i < nb_procs; i++ ){
			MPI_Recv(&sol, 1, MPI_LONG_DOUBLE, i, i, MPI_COMM_WORLD, NULL);
			final_sol = final_sol + sol;
		}
		cout << "final sol: " << final_sol << endl;
		cout << "exact sol: " << M_PI << endl;
	}


	MPI_Finalize();
	return 0;
}
