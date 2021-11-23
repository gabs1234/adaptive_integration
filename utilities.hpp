#pragma once
#include <iostream>

/* utilities */
template <class T>
void print_list(int n, T* list){
	std::cout << "[ ";
	for (size_t i = 0; i < n; i++) {
		std::cout << list[i] << " ";
	}
	std::cout << "]" << std::endl;
}
