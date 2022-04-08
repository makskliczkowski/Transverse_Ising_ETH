#!/bin/bash
sbatch --cpus-per-task=24 --mem=12G --time=24:00:00 --array=0-319%20 Ising_spectrals.sh 14 $1 6 0.05 24
sbatch --cpus-per-task=24 --mem=12G --time=24:00:00 --array=0-319%20 Ising_spectrals.sh 14 $1 6 0.05 24
exit;

#SIGMA^Z_j=..
#sbatch --cpus-per-task=8 --mem=4G --time=3:00:00 --array=0-239%20 Ising_spectrals.sh 10 $1 0 0.05 8
#sbatch --cpus-per-task=16 --mem=8G --time=3:00:00 --array=0-239%20 Ising_spectrals.sh 11 $1 0 0.05 16
#sbatch --cpus-per-task=16 --mem=8G --time=3:00:00 --array=0-279%20 Ising_spectrals.sh 12 $1 0 0.05 16
#sbatch --cpus-per-task=20 --mem=16G --time=48:00:00 --array=0-279%20 Ising_spectrals.sh 13 $1 0 0.05 20
#sbatch --cpus-per-task=20 --mem=20G --time=48:00:00 --array=0-319%20 Ising_spectrals.sh 14 $1 0 0.05 20
#sbatch --cpus-per-task=40 --mem=40G --time=96:00:00 --array=0-319%20 Ising_spectrals.sh 15 $1 0 0.05 40
#sbatch --cpus-per-task=40 --mem=120G --time=168:00:00 --array=0-359%20 Ising_spectrals.sh 16 $1 0 0.05 40

#SIGMA^Z_q=..
#sbatch --cpus-per-task=8 --mem=1G --time=3:00:00 --array=0-119%20 --constraint=rack-6 Ising_spectrals.sh 10 $1 3 0.05 8
#sbatch --cpus-per-task=8 --mem=2G --time=3:00:00 --array=0-119%20 --constraint=rack-6 Ising_spectrals.sh 11 $1 3 0.05 8
#sbatch --cpus-per-task=8 --mem=4G --time=3:00:00 --array=0-139%20 --constraint=rack-6 Ising_spectrals.sh 12 $1 3 0.05 8
#sbatch --cpus-per-task=20 --mem=8G --time=48:00:00 --array=0-139%20 --constraint=rack-6 Ising_spectrals.sh 13 $1 3 0.05 20
#sbatch --cpus-per-task=32 --mem=16G --time=48:00:00 --array=0-159%20 --constraint=rack-6 Ising_spectrals.sh 14 $1 3 0.05 32
sbatch --cpus-per-task=40 --mem=60G --time=96:00:00 --array=0-159%20 --constraint=rack-6 Ising_spectrals.sh 15 $1 3 0.05 40
sbatch --cpus-per-task=40 --mem=150G --time=168:00:00 --array=0-179%20 --constraint=rack-6 Ising_spectrals.sh 16 $1 3 0.05 40


#H_j=..
#sbatch --mem=1G --time=3:00:00 --array=0-119%5 Ising_spectrals.sh 10 $1 2 0.1
#sbatch --mem=2G --time=3:00:00 --array=0-119%5 Ising_spectrals.sh 11 $1 2 0.1
#sbatch --mem=4G --time=3:00:00 --array=0-139%5 Ising_spectrals.sh 12 $1 2 0.1
#sbatch --mem=8G --time=48:00:00 --array=0-139%10 --constraint=rack-6 Ising_spectrals.sh 13 $1 2 0.1
#sbatch --mem=16G --time=48:00:00 --array=0-159%10 --constraint=rack-6 Ising_spectrals.sh 14 $1 2 0.1
#sbatch --mem=36G --time=48:00:00 --array=0-159%5 --constraint=rack-6 Ising_spectrals.sh 15 $1 2 0.1
#sbatch --mem=90G --time=168:00:00 --array=0-179%5 --constraint=rack-6 Ising_spectrals.sh 16 $1 2 0.1

#H_q=..
#sbatch --cpus-per-task=8 --mem=1G --time=3:00:00 --array=0-119%20 --constraint=rack-6 Ising_spectrals.sh 10 $1 5 0.05 8
#sbatch --cpus-per-task=8 --mem=4G --time=3:00:00 --array=0-119%20 --constraint=rack-6 Ising_spectrals.sh 11 $1 5 0.05 8
#sbatch --cpus-per-task=8 --mem=4G --time=3:00:00 --array=0-139%20 --constraint=rack-6 Ising_spectrals.sh 12 $1 5 0.05 8
#sbatch --cpus-per-task=20 --mem=8G --time=48:00:00 --array=0-139%20 --constraint=rack-6 Ising_spectrals.sh 13 $1 5 0.05 20
#sbatch --cpus-per-task=32 --mem=16G --time=48:00:00 --array=0-159%20 --constraint=rack-6 Ising_spectrals.sh 14 $1 5 0.05 32
sbatch --cpus-per-task=40 --mem=60G --time=96:00:00 --array=0-159%20 --constraint=rack-6 Ising_spectrals.sh 15 $1 5 0.05 40
sbatch --cpus-per-task=40 --mem=150G --time=168:00:00 --array=0-179%20 --constraint=rack-6 Ising_spectrals.sh 16 $1 5 0.05 40