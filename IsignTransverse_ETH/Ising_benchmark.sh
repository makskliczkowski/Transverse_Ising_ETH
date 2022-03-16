#!/bin/bash
#SBATCH --job-name=benchmark
#SBATCH --cpus-per-task=64
#SBATCH --mem=300G
#SBATCH --time=300:00:00
#SBATCH --constraint=rack-6

#export MKL_NUM_THREADS=64
		
module purge
module load OpenBLAS/0.3.18-GCC-11.2.0
module load intel/2022.00
#module load foss/2021b
#module load gomkl/2021a
icc main.cpp IsingModel.cpp IsingModel_disorder.cpp IsingModel_sym.cpp tools.cpp user_interface.cpp -o Ising_bench.o\
 -I../LIBRARIES_CPP/armadillo-10.8.2/include/ -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -lblas -lopenblas\
 -pthread -lstdc++fs -qopenmp -qmkl=parallel -std=c++17 -O3 >& compile_benchmark.log

#gcc main.cpp IsingModel.cpp IsingModel_disorder.cpp IsingModel_sym.cpp tools.cpp user_interface.cpp -o Ising_bench.o\
# -I../LIBRARIES_CPP/armadillo-10.8.2/include/ -llapack -lopenblas\
# -lrt -pthread -lstdc++fs -fopenmp -fomit-frame-pointer -std=c++17 -O3 >& compile_benchmark.log

now="$(date +'date=%d_%m_time=%H_%M')"
./Ising_bench.o -L 15 -Ln 8 -Ls 1 -g 0.6 -h 0.8 -w 0.01 -th 64 -m 0 -w 0.01 -r 1 -op 0 -fun 5 -s 1 -b 0 >& benchmark_${now}.txt
