#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --mem=1000
#SBATCH --cpus-per-task=1
#SBATCH --job-name=compile

module purge
#unset MKL_SERIAL
#export CPATH="/usr/include/hdf5/serial"
#export LD_LIBRARY_PATH="${MKLROOT}/lib/intel64"
export OMP_NUM_THREADS=16

#module load Armadillo/9.900.1-foss-2020a
module load OpenBLAS/0.3.18-GCC-11.2.0
module load foss/2021b
#module load intel/2021a

#g++ main.cpp IsingModel.cpp IsingModel_disorder.cpp IsingModel_sym.cpp tools.cpp user_interface.cpp -o Ising.o\
#  -I ../LIBRARIES_CPP/armadillo-9.850.1/include/ -larmadillo -pthread -fopenmp -std=c++20 -std=c++11 -O3 >& log_compile.txt


g++ -std=c++2a main.cpp IsingModel.cpp IsingModel_disorder.cpp IsingModel_sym.cpp tools.cpp user_interface.cpp -o Ising.o\
 -DARMA_DONT_USE_WRAPPER -I../LIBRARIES_CPP/armadillo-10.8.0/include/ -llapack -lopenblas -fopenmp -lpthread -lm -lstdc++fs -fomit-frame-pointer -O2 >& log_compile.txt

#clang++ main.cpp IsingModel.cpp IsingModel_disorder.cpp IsingModel_sym.cpp tools.cpp user_interface.cpp -o Ising.o -pthread\
#  -I ../LIBRARIES_CPP/armadillo-9.850.1/include/ -L${MKLROOT}/lib/intel64  -fopenmp -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -lstdc++fs -llapack -fomit-frame-pointer -lblas -std=c++17 -Ofast

#icpc main.cpp IsingModel.cpp IsingModel_disorder.cpp IsingModel_sym.cpp tools.cpp user_interface.cpp -o Ising.o\
# -I../LIBRARIES_CPP/armadillo-10.8.2/include/ -L${MKLROOT}/lib/intel64 -DARMA_DONT_USE_WRAPPER -lblas -llapack\
# -lmkl_rt -lpthread -lm -ldl -lstdc++fs -qopenmp -mkl=parallel -std=c++20 -O3
