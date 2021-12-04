#pragma once
#include "utilities.hpp"
#include <cmath>

#define MAX_ITER 50 //
#define DF_STEP .1 //

/* Gauss Legendre integration functions */
template <class T>
T PO(T x){
	return 1;
}
template <class T>
T P1(T x){
	return x;
}

template <class T>
T P(int n, T x){
	if( n == 0 ){
		return PO<T>(x);
	}
	else if( n == 1 ){
		return P1<T>(x);
	}

	return (P<T>(n-1, x)*x*(2*n - 1) - P<T>(n-2, x)*(n-1))/n;
}

template <class T>
T dP(int n, T x){
	return (P<T>(n-1, x) - x*P<T>(n, x))*n/(1-x*x);
}

template <class T>
T newton_raphson(T x0, int n, T eps){
	T xa = x0;
	T xb = xa - P<T>(n, xa)/dP<T>(n, xa);
	while( std::abs(xb - xa) > eps ){
		xa = xb;
		xb = xa - P<T>(n, xa)/dP<T>(n, xa);
	}
	return xb;
}

template <class T>
T legendre_root_est(int k, int n){
	T res = std::cos( (4*k-1)*M_PI/(4*n+2) );
	return res;
}

template <class T>
void legendre_roots(int n, T* roots, T eps){
	T x0;
	int stop = n >> 1;

	// Take into account the odd n symmetry
	if( n%2 != 0 ){
		roots[stop] = 0.;
	}

	// Thanks to the symmetry of the roots, we only calculate the half of them
	for( int k = 0; k < stop; k++ ){
		x0 = legendre_root_est<T>(k+1, n);
		roots[k] = newton_raphson<T>(x0, n, eps);
		roots[n-k-1] = -roots[k];
	}
}

template <class T>
void gauss_legendre_weights(int n, T* roots, T* weight){
	int stop = n/2;

	for( int i = 0; i <= stop; i++ ){
		weight[i] = 2/( (1-std::pow(roots[i], 2.0))*std::pow(dP<T>(n,roots[i]), 2.0) );
		weight[n-i-1] = weight[i];
	}
}

template <class T>
T gauss_legendre(float a, float b, T (*f)(T), int n, T* x, T* w, T eps){
	float hm = b - a;
	float hp = b + a;

	// Calculate integral over interval a, b
	T sol = 0;
	for( int i = 0; i < n; i++ ){
		sol = sol + w[i]*(*f)( .5*hm*x[i] + .5*hp );
	}

	return .5*hm*sol;
}

template <class T>
T adaptive_integration_rec(float a, float b, T (*f)(T), int n, T* x, T* w, T eps, int iteration){
	// TODO: find solution to not have to re declare variables each time...
	float half = (a+b)/2.;
	T Q_curr = gauss_legendre<T>(a, b, f, n, x, w, eps);
	T Q_next = gauss_legendre<T>(a, half, f, n, x, w, eps) + gauss_legendre<T>(half, b, f, n, x, w, eps);
	T err = std::abs(Q_curr - Q_next);

	if( iteration > MAX_ITER || err < eps ){
		return Q_next;
	}

	return adaptive_integration_rec<T>(a, half, f, n, x, w, eps, iteration + 1) + adaptive_integration_rec<T>(half, b, f, n, x, w, eps, iteration + 1);
}


template <class T>
T adaptive_integration(float a, float b, T (*f)(T), int n, T eps){
	T x[n], w[n];


	// Calculate the roots
	legendre_roots<T>(n, x, eps);;
	// Calculate the weights
	gauss_legendre_weights<T>(n, x, w);

	return adaptive_integration_rec<T>(a, b, f, n, x, w, eps, (int)0);
}

/* Load balancing functions */

template <class T>
void derivative(T (*f)(T), int nb_points, T a, T b, T* df){
	T h = DF_STEP;
	for( int i = 0; i < nb_points; i++ ){
		// or store points in advance ?
		df[i] = ( f((i+1)*h+a)-f(i*h+a) )/h;
	}
}

template <class T>
void sub_intervals(T a, T b, int nb_procs, int nb_points, T* df, T* interval){
	// Divide the interval into nb_procs sub-intervals
	int sep = std::round(nb_points/nb_procs);

	// get sums per sub interval
	T sums[nb_procs];
	for( int k = 0; k < nb_procs; k++ ){
		sums[k] = 0;
		for( int i = k*sep; i < (k+1)*sep; i++ ){
			if( i >= nb_points ){
				break;
			}
			sums[k] = sums[k] + std::abs(df[i]);
		}
	}

	// Get total sum
	T sum = 0;
	for( int i = 0; i < nb_points; i++ ){
		sum = sum + std::abs(df[i]);
	}

	// Get weights for each sub interval
	interval[0] = a;
	T len = (b-a);
	T k = len/((nb_procs-1)*sum);
	// for( int i = 0; i < nb_procs; i++ ){
	// 	interval[i+1] = sums[i]/sum;
	// }
	for( int i = 0; i < nb_procs; i++ ){
		interval[i+1] = k*(sum-sums[i]);
	}

	// Use weights to deduce sub intervals
	for( int i = 1; i < nb_procs+1; i++ ){
		interval[i] = interval[i] + interval[i-1];
	}
	interval[nb_procs] = b;
	// std::cout << "interval list: ";
	// print_list<T>(nb_procs+1, interval);
}

template <class T>
void load_balance(T a, T b, int nb_procs, T (*f)(T), T* interval){
	// approximate derivative of f
	int nb_points = (int)((b-a)/(DF_STEP));
	T df[nb_points];
	derivative<T>(f, nb_points, a, b, df);

	// Deduce balanced partions in interval [a, b] according to f
	interval[0] = a;
	interval[nb_procs] = b;
	sub_intervals<T>(a, b, nb_procs, nb_points, df, interval);
}
