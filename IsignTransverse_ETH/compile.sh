#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --mem=1000
#SBATCH --cpus-per-task=1
#SBATCH --job-name=compile

export MKL_SERIAL=yes
export CPATH="/usr/include/hdf5/serial"
export LD_LIBRARY_PATH="${MKLROOT}/lib/intel64"

module purge
module load Armadillo/9.900.1-foss-2020a
module load intel/2021b

#module load Clang/11.0.1-GCCcore-10.3.0
#module load foss/2021b

#g++ main.cpp IsingModel.cpp IsingModel_disorder.cpp IsingModel_sym.cpp tools.cpp user_interface.cpp -o Ising.o\
#  -I ../LIBRARIES_CPP/armadillo-9.850.1/include/ -larmadillo -pthread -fopenmp -std=c++20 -std=c++11 -O3 >& log_compile.txt

g++ main.cpp IsingModel.cpp IsingModel_disorder.cpp IsingModel_sym.cpp tools.cpp user_interface.cpp -o Ising.o -pthread\
  -I../LIBRARIES_CPP/armadillo-9.850.1/include/ -L${MKLROOT}/lib/intel64 -fopenmp -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -lstdc++fs -llapack -fcx-fortran-rules -fomit-frame-pointer -lblas -std=c++17 -O3

#clang++ main.cpp IsingModel.cpp IsingModel_disorder.cpp IsingModel_sym.cpp tools.cpp user_interface.cpp -o Ising.o -pthread\
#  -I ../LIBRARIES_CPP/armadillo-9.850.1/include/ -L${MKLROOT}/lib/intel64 -fopenmp -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -lstdc++fs -llapack -fomit-frame-pointer -lblas -std=c++17 -Ofast

#icc main.cpp IsingModel.cpp IsingModel_disorder.cpp IsingModel_sym.cpp tools.cpp user_interface.cpp -o Ising.o\
# -I ../LIBRARIES_CPP/armadillo-9.850.1/include/ -L${MKLROOT}/lib/intel64 -lmkl_rt -lpthread -lm -ldl -lstdc++fs -qopenmp -larmadillo -qmkl=parallel -std=c++17 -O3
