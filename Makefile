# TODO: make a good makefile

main: src/main.cpp src/functions.hpp src/utilities.hpp
	g++ src/main.cpp -o bin/main -lm

mpi_main: src/para_main.cpp src/para_functions.hpp src/utilities.hpp
	mpic++ src/para_main.cpp -o bin/para_main -lm

plot: src/legendre_functions.cpp src/functions.hpp src/utilities.hpp
	g++ src/legendre_functions.cpp -o bin/plot_legendre -lm
	for i in {2..12}; do \
		./plot_legendre $$i 500 >> data/legendre_function_$$i.dat; \
	done

.PHONY: clean_plot

clean_plot:
	rm data/legendre_function_*
