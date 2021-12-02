#include "utilities.hpp"
#include <cmath>

#define MAX_ITER 50

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
void legendre_roots(int n, T roots[], T eps){
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
void gauss_legendre_weights(int n, T roots[], T weight[]){
	int stop = n/2;

  for( int i = 0; i <= stop; i++ ){
		weight[i] = 2/( (1-std::pow(roots[i], 2.0))*std::pow(dP<T>(n,roots[i]), 2.0) );
		weight[n-i-1] = weight[i];
	}
}

template <class T>
T gauss_legendre(float a, float b, T (*f)(T), int n, T (&x)[], T (&w)[], T eps){
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
T adaptive_integration_rec(float a, float b, T (*f)(T), int n, T (&x)[], T (&w)[], T eps, int iteration){
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

	std::cout << "hello!" << std::endl;

	// Calculate the roots
  legendre_roots<T>(n, x, eps);
	std::cout << "roots calculated!" << std::endl;
	// Calculate the weights
  gauss_legendre_weights<T>(n, x, w);
	std::cout << "weights calculated!" << std::endl;

	print_list(n, x);
	print_list(n, w);

	return adaptive_integration_rec<T>(a, b, f, n, x, w, eps, (int)0);
}
