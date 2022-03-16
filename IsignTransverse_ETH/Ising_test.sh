#!/bin/bash
#SBATCH --job-name=test
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=3:00:00

module purge
#unset MKL_SERIAL
#export CPATH="/usr/include/hdf5/serial"
#export LD_LIBRARY_PATH="/share/apps/eb2/software/imkl/2022.0.1/mkl/2022.0.1/lib/intel64"
#export LIBRARY_PATH="/home/rswietek/LIBRARIES_CPP/armadillo-10.8.2/lib64"
#export OMP_NUM_THREADS=64

#module load Armadillo/9.900.1-foss-2020a
#module load foss/2021a

#module load OpenBLAS/0.3.18-GCC-11.2.0
module load intel/2022

#icc -g main.cpp IsingModel.cpp IsingModel_disorder.cpp IsingModel_sym.cpp tools.cpp user_interface.cpp -o Ising_test.o\
# -I../LIBRARIES_CPP/armadillo-10.8.2/include/ -llapack -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -pthread -lstdc++fs -qopenmp -qmkl=parallel -std=c++17 -O0 >& compile_test.log

icc -g ARMA_TEST.cpp -o arma_test.o\
 -I../LIBRARIES_CPP/armadillo-10.8.2/include/ -lopenblas -llapack -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread\
 -pthread -lstdc++fs -qopenmp -qmkl=parallel -std=c++17 -O0 >& arma_compile_test.log

#g++ -g ARMA_TEST.cpp -o arma_test.o -m64 -fcx-fortran-rules -fomit-frame-pointer -std=c++20 \
# -pthread -fopenmp -I/home/rswietek/LIBRARIES_CPP/armadillo-10.8.2/include\
# -lpthread -lm -ldl -lstdc++fs -O0 >& arma_compile_test.log

# main.cpp IsingModel.cpp IsingModel_disorder.cpp IsingModel_sym.cpp tools.cpp user_interface.cpp -o Ising.o -O3 -march=native >& compile_test.log

#-lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -DARMA_DONT_USE_WRAPPER

#g++ -std=c++2a main.cpp IsingModel.cpp IsingModel_disorder.cpp IsingModel_sym.cpp tools.cpp user_interface.cpp -o Ising.o\
# -DARMA_DONT_USE_WRAPPER -I/home/rswietek/LIBRARIES_CPP/armadillo-10.8.2/include -llapack -lopenblas -fopenmp -lpthread -lm -lstdc++fs -fomit-frame-pointer -Ofast >& compile_test.log
 
now="$(date +'date=%d_%m_time=%H_%M')"
#valgrind --leak-check=yes --track-origins=yes ./Ising_test.o -L 6 -Ln 1 -Ls 1 -g 0.9 -h 0.8 -w 0.01 -th 64 -m 0 -w 0.01 -r 1 -op 0 -fun 0 -s 1 -b 0 >& run_test_${now}.log
valgrind --leak-check=yes --track-origins=yes ./arma_test.o >& arma_run_test_${now}.log
