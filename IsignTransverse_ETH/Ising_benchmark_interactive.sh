#!/bin/bash
#SBATCH --job-name=benchmark
#SBATCH --cpus-per-task=64
#SBATCH --mem=200G
#SBATCH --time=300:00:00

module purge
#module load OpenBLAS/0.3.18-GCC-11.2.0
module load intel/2022.00

icc main.cpp IsingModel.cpp IsingModel_disorder.cpp IsingModel_sym.cpp tools.cpp user_interface.cpp -o Ising_bench.o\
 -I../LIBRARIES_CPP/armadillo-10.8.2/include/ -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -pthread -lstdc++fs -qopenmp -qmkl=parallel -std=c++17 -O3 >& compile_genchmark.log

lscpu > benchmark.txt
echo "#cores\t\tchain length\t\tdim\t\twith vec 'dc' [s]\t\twith vec 'std' [s]\t\twithout vec [s]" >> benchmark.txt
for L in {8..16}; do
	for nth in 1 2 4 8 12 16 20 24 32 40 48 64; do 
		export OMP_NUM_THREADS=$nth
		export MKL_NUM_THREADS=$nth
		./Ising_bench.o -L $L -Ln 1 -Ls 1 -g 0.9 -h 0.8 -w 0.01 -th $nth -m 0 -w 0.01 -r 1 -op 0 -fun 5 -s 1 -b 0 >> benchmark.txt
	done
	"\n" >> benchmark.txt
done
