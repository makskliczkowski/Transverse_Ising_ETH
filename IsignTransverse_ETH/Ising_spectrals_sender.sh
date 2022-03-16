#!/bin/bash
h=$1	

#SIGMA^Z_j=..
#sbatch --cpus-per-task=8 --mem=4G --time=3:00:00 --array=0-239%20 Ising_spectrals.sh 10 $h 0 0.05 8
#sbatch --cpus-per-task=16 --mem=8G --time=3:00:00 --array=0-239%20 Ising_spectrals.sh 11 $h 0 0.05 16
#sbatch --cpus-per-task=16 --mem=8G --time=3:00:00 --array=0-279%20 Ising_spectrals.sh 12 $h 0 0.05 16
#sbatch --cpus-per-task=20 --mem=16G --time=48:00:00 --array=0-279%20 Ising_spectrals.sh 13 $h 0 0.05 20
#sbatch --cpus-per-task=20 --mem=20G --time=48:00:00 --array=0-319%20 Ising_spectrals.sh 14 $h 0 0.05 20
#sbatch --cpus-per-task=40 --mem=40G --time=96:00:00 --array=0-319%20 Ising_spectrals.sh 15 $h 0 0.05 40
#sbatch --cpus-per-task=40 --mem=120G --time=168:00:00 --array=0-359%20 Ising_spectrals.sh 16 $h 0 0.05 40

#SIGMA^Z_q=..
#sbatch --cpus-per-task=8 --mem=1G --time=3:00:00 --array=0-119%20 --constraint=rack-6 Ising_spectrals.sh 10 $h 3 0.05 8
#sbatch --cpus-per-task=8 --mem=2G --time=3:00:00 --array=0-119%20 --constraint=rack-6 Ising_spectrals.sh 11 $h 3 0.05 8
#sbatch --cpus-per-task=8 --mem=4G --time=3:00:00 --array=0-139%20 --constraint=rack-6 Ising_spectrals.sh 12 $h 3 0.05 8
#sbatch --cpus-per-task=20 --mem=8G --time=48:00:00 --array=0-139%20 --constraint=rack-6 Ising_spectrals.sh 13 $h 3 0.05 20
#sbatch --cpus-per-task=32 --mem=16G --time=48:00:00 --array=0-159%20 --constraint=rack-6 Ising_spectrals.sh 14 $h 3 0.05 32
sbatch --cpus-per-task=40 --mem=60G --time=96:00:00 --array=0-159%20 --constraint=rack-6 Ising_spectrals.sh 15 $h 3 0.05 40
sbatch --cpus-per-task=40 --mem=150G --time=168:00:00 --array=0-179%20 --constraint=rack-6 Ising_spectrals.sh 16 $h 3 0.05 40


#H_j=..
#sbatch --mem=1G --time=3:00:00 --array=0-119%5 Ising_spectrals.sh 10 $h 2 0.1
#sbatch --mem=2G --time=3:00:00 --array=0-119%5 Ising_spectrals.sh 11 $h 2 0.1
#sbatch --mem=4G --time=3:00:00 --array=0-139%5 Ising_spectrals.sh 12 $h 2 0.1
#sbatch --mem=8G --time=48:00:00 --array=0-139%10 --constraint=rack-6 Ising_spectrals.sh 13 $h 2 0.1
#sbatch --mem=16G --time=48:00:00 --array=0-159%10 --constraint=rack-6 Ising_spectrals.sh 14 $h 2 0.1
#sbatch --mem=36G --time=48:00:00 --array=0-159%5 --constraint=rack-6 Ising_spectrals.sh 15 $h 2 0.1
#sbatch --mem=90G --time=168:00:00 --array=0-179%5 --constraint=rack-6 Ising_spectrals.sh 16 $h 2 0.1

#H_q=..
#sbatch --cpus-per-task=8 --mem=1G --time=3:00:00 --array=0-119%20 --constraint=rack-6 Ising_spectrals.sh 10 $h 5 0.05 8
#sbatch --cpus-per-task=8 --mem=4G --time=3:00:00 --array=0-119%20 --constraint=rack-6 Ising_spectrals.sh 11 $h 5 0.05 8
#sbatch --cpus-per-task=8 --mem=4G --time=3:00:00 --array=0-139%20 --constraint=rack-6 Ising_spectrals.sh 12 $h 5 0.05 8
#sbatch --cpus-per-task=20 --mem=8G --time=48:00:00 --array=0-139%20 --constraint=rack-6 Ising_spectrals.sh 13 $h 5 0.05 20
#sbatch --cpus-per-task=32 --mem=16G --time=48:00:00 --array=0-159%20 --constraint=rack-6 Ising_spectrals.sh 14 $h 5 0.05 32
sbatch --cpus-per-task=40 --mem=60G --time=96:00:00 --array=0-159%20 --constraint=rack-6 Ising_spectrals.sh 15 $h 5 0.05 40
sbatch --cpus-per-task=40 --mem=150G --time=168:00:00 --array=0-179%20 --constraint=rack-6 Ising_spectrals.sh 16 $h 5 0.05 40