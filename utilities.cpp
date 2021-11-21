#include "headers.h"

using namespace std;

void print_list(int n, long double list[]){
	cout << "[ ";
	for (size_t i = 0; i < n; i++) {
		cout << list[i] << " ";
	}
	cout << "]" << endl;
}
