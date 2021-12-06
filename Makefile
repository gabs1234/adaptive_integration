# TODO: use variables

main: src/main.cpp src/functions.hpp src/utilities.hpp
	@g++ src/main.cpp -o bin/main -lm

mpi_main: src/mpi_main.cpp src/mpi_functions.hpp src/utilities.hpp
	@mpic++ src/mpi_main.cpp -o bin/mpi_main -lm

mpi_bench: src/mpi_benchmark.cpp src/mpi_functions.hpp src/utilities.hpp
	@mpic++ src/mpi_benchmark.cpp -o bin/mpi_bench -lm

plot: src/legendre_functions.cpp src/functions.hpp src/utilities.hpp
	g++ src/legendre_functions.cpp -o bin/plot_legendre -lm
	for i in {2..12}; do \
		./plot_legendre $$i 500 >> data/legendre_function_$$i.dat; \
	done

.PHONY: clean_plot

clean_plot:
	rm data/legendre_function_*
