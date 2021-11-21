#include "headers.h"

using namespace std;

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

	return P<T>(n-1, x)*x*(2*n - 1)/n - P<T>(n-2, x)*(n-1)/n;
}

template <class T>
T dP(int n, T x){
	return (P<T>(n-1, x) - x*P<T>(n, x))*n/(1-x*x);
}

template <class T>
T newton_raphson(T x0, int n, T eps){
	T xa = x0;
  T xb = xa - P(n, xa)/dP(n, xa);

  while( abs(xb - xa) > eps ){
		xa = xb;
    xb = xa - P(n, xa)/dP(n, xa);
	}

  return xb;
}

template <class T>
T legendre_root_est(int k, int n){
	T res = cos(4*k-1)*M_PI/(4*n+2);
	return res;
}

template <class T>
void legendre_roots(int n, T roots[], T eps){
  T x0;
	int stop = n/2;

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
		weight[i] = 2/( (1-pow(roots[i], 2))*pow(dP<T>(n,roots[i]),2) );
		weight[n-i-1] = weight[i];
	}
}

template <class T>
T gauss_legendre(float a, float b, T (*f)(T), int n, float eps){
	T x[n], w[n];

	// Calculate the roots
  legendre_roots<T>(n, x, eps);
	// Calculate the weights
  gauss_legendre_weights<T>(n, x, w);

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
T adaptive_integration(float a, float b, T (*f)(T), int n, T eps, int iteration){
	// TODO: find solution to not have to re declare variables each time...
	float half = (a+b)/2.;
  T Q_curr = gauss_legendre<T>(a, b, f, n, eps);
  T Q_next = gauss_legendre<T>(a, half, f, n, eps) + gauss_legendre<T>(half, b, f, n, eps);
  T err = abs(Q_curr - Q_next);

  if( iteration > MAX_ITER || err < eps ){
    return Q_next;
	}

  return adaptive_integration<T>(a, half, f, n, eps, iteration + 1) + adaptive_integration<T>(half, b, f, n, eps, iteration + 1);
}
