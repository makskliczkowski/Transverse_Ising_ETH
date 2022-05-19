#!/bin/bash

x=$1
y=$2	
fun=$3
num=3	 #minimum number of required input
echo $1, $2, $3, $4, $5, $6, $7
operator=0; site=0; jobid=0;
name=parameter_scaling;
if [[ $# -lt 3 ]]; then   
  echo "Too few input parameters! Required ${num}"
  echo "INPUTS TO SCRIPT ARE: h,  fun, operator, site"
  #exit 1
elif [[ $# == 4 ]]; then
  operator=$4;
elif [[ $# == 5 ]]; then
  operator=$4; site=$5;  
elif [[ $# == 6 ]]; then
  operator=$4; site=$5;  name=$6
else
  operator=$4; site=$5;  name=$6; jobid=$7;
fi
echo $x, $y, $fun, $operator, $site, $name, $jobid, "${@:8}"

#sbatch --job-name=$name --output=logs/${name}_log-%A-%a.out --cpus-per-task=20 --mem=2G --time=168:00:00 --array=0-23%4 Ising_parameter_scaling.sh 16 $x $y 0.05 $fun $operator $site 20 $jobid "${@:8}"
#sbatch --job-name=$name --output=logs/${name}_log-%A-%a.out --cpus-per-task=20 --mem=2G --time=168:00:00 --array=0-23%4 Ising_parameter_scaling.sh 17 $x $y 0.05 $fun $operator $site 20 $jobid "${@:8}"
#sbatch --job-name=$name --output=logs/${name}_log-%A-%a.out --cpus-per-task=20 --mem=4G --time=168:00:00 --array=0-23%4 Ising_parameter_scaling.sh 18 $x $y 0.05 $fun $operator $site 20 $jobid "${@:8}"
#sbatch --job-name=$name --output=logs/${name}_log-%A-%a.out --cpus-per-task=20 --mem=6G --time=168:00:00 --array=0-23%4 Ising_parameter_scaling.sh 19 $x $y 0.05 $fun $operator $site 20 $jobid "${@:8}"
#sbatch --job-name=$name --output=logs/${name}_log-%A-%a.out --cpus-per-task=20 --mem=10G --time=168:00:00 --array=0-23%4 Ising_parameter_scaling.sh 20 $x $y 0.05 $fun $operator $site 20 $jobid "${@:8}"
#sbatch --job-name=$name --output=logs/${name}_log-%A-%a.out --cpus-per-task=32 --mem=20G --time=168:00:00 --array=0-23%4 Ising_parameter_scaling.sh 21 $x $y 0.05 $fun $operator $site 32 $jobid "${@:8}"
#sbatch --job-name=$name --output=logs/${name}_log-%A-%a.out --cpus-per-task=32 --mem=40G --time=168:00:00 --array=0-23%4 Ising_parameter_scaling.sh 22 $x $y 0.05 $fun $operator $site 32 $jobid "${@:8}"
#sbatch --job-name=$name --output=logs/${name}_log-%A-%a.out --cpus-per-task=32 --mem=60G --time=168:00:00 --array=0-23%6 Ising_parameter_scaling.sh 23 $x $y 0.05 $fun $operator $site 32 $jobid "${@:8}"
#sbatch --job-name=$name --output=logs/${name}_log-%A-%a.out --cpus-per-task=32 --mem=120G --time=168:00:00 --array=0-23%6 Ising_parameter_scaling.sh 24 $x $y 0.05 $fun $operator $site 32 $jobid "${@:8}"
#exit;

# #SBATCH --output=logs/g_scale_log-%j-%a.out
#sbatch --job-name="${name}_L=8" --output=logs/${name}_log-%A-%a.out --cpus-per-task=1 --mem=500 --time=1:00:00 --array=0-19%1 Ising_parameter_scaling.sh 8 $x $y 0.05 $fun $operator $site 1 $jobid "${@:8}"
#sbatch --job-name="${name}_L=9" --output=logs/${name}_log-%A-%a.out --cpus-per-task=1 --mem=500 --time=1:00:00 --array=0-19%1 Ising_parameter_scaling.sh 9 $x $y 0.05 $fun $operator $site 1 $jobid "${@:8}"
#sbatch --job-name="${name}_L=10" --output=logs/${name}_log-%A-%a.out --cpus-per-task=1 --mem=1G --time=1:00:00 --array=0-19%4 Ising_parameter_scaling.sh 10 $x $y 0.05 $fun $operator $site 1 $jobid "${@:8}"
#sbatch --job-name="${name}_L=11" --output=logs/${name}_log-%A-%a.out --cpus-per-task=1 --mem=1G --time=3:00:00 --array=0-19%4 Ising_parameter_scaling.sh 11 $x $y 0.05 $fun $operator $site 1 $jobid "${@:8}"
#sbatch --job-name="${name}_L=12" --output=logs/${name}_log-%A-%a.out --cpus-per-task=2 --mem=2G --time=24:00:00 --array=0-19%4 Ising_parameter_scaling.sh 12 $x $y 0.05 $fun $operator $site 2 $jobid "${@:8}"
#sbatch --job-name="${name}_L=13" --output=logs/${name}_log-%A-%a.out --cpus-per-task=4 --mem=8G --time=72:00:00 --array=0-19%4 Ising_parameter_scaling.sh 13 $x $y 0.05 $fun $operator $site 4 $jobid "${@:8}"
#sbatch --job-name="${name}_L=14" --output=logs/${name}_log-%A-%a.out --cpus-per-task=4 --mem=16G --time=72:00:00 --array=0-19%4 Ising_parameter_scaling.sh 14 $x $y 0.05 $fun $operator $site 4 $jobid "${@:8}"
#sbatch --job-name="${name}_L=15" --output=logs/${name}_log-%A-%a.out --cpus-per-task=4 --mem=40G --time=168:00:00 --array=0-19%3 Ising_parameter_scaling.sh 15 $x $y 0.05 $fun $operator $site 8 $jobid "${@:8}"
#sbatch --job-name="${name}_L=16" --output=logs/${name}_log-%A-%a.out --cpus-per-task=4 --mem=70G --time=200:00:00 --array=0-19%4 Ising_parameter_scaling.sh 16 $x $y 0.05 $fun $operator $site 8 $jobid "${@:8}"
#exit;
#symmetries
sbatch --job-name="${name}_L=14_sym" --output=logs/${name}_log-%A-%a.out --cpus-per-task=1 --mem=1G --time=3:00:00 --array=0-19%2 Ising_parameter_scaling.sh 14 $x $y 0.05 $fun $operator $site 1 $jobid "${@:8}"
sbatch --job-name="${name}_L=15_sym" --output=logs/${name}_log-%A-%a.out --cpus-per-task=1 --mem=1G --time=3:00:00 --array=0-19%2 Ising_parameter_scaling.sh 15 $x $y 0.05 $fun $operator $site 1 $jobid "${@:8}"
sbatch --job-name="${name}_L=16_sym" --output=logs/${name}_log-%A-%a.out --cpus-per-task=2 --mem=2G --time=24:00:00 --array=0-19%2 Ising_parameter_scaling.sh 16 $x $y 0.05 $fun $operator $site 2 $jobid "${@:8}"
sbatch --job-name="${name}_L=17_sym" --output=logs/${name}_log-%A-%a.out --cpus-per-task=2 --mem=8G --time=72:00:00 --array=0-19%2 Ising_parameter_scaling.sh 17 $x $y 0.05 $fun $operator $site 2 $jobid "${@:8}"
sbatch --job-name="${name}_L=18_sym" --output=logs/${name}_log-%A-%a.out --cpus-per-task=4 --mem=16G --time=72:00:00 --array=0-19%2 Ising_parameter_scaling.sh 18 $x $y 0.05 $fun $operator $site 4 $jobid "${@:8}"
sbatch --job-name="${name}_L=19_sym" --output=logs/${name}_log-%A-%a.out --cpus-per-task=4 --mem=20G --time=72:00:00 --array=0-19%2 Ising_parameter_scaling.sh 18 $x $y 0.05 $fun $operator $site 4 $jobid "${@:8}"


#sbatch --mem=8G --time=3:00:00 --array=0-29%15 Ising_parameter_scaling.sh 10 $x $y 0.05 $fun
#sbatch --mem=16G --time=24:00:00 --array=0-29%15 Ising_parameter_scaling.sh 11 $x $y 0.05 $fun
#sbatch --mem=16G --time=24:00:00 --array=0-29%15 Ising_parameter_scaling.sh 12 $x $y 0.05 $fun
#sbatch --mem=20G --time=24:00:00 --array=0-29%15 Ising_parameter_scaling.sh 13 $x $y 0.05 $fun
#sbatch --mem=20G --time=96:00:00 --array=0-29%15 Ising_parameter_scaling.sh 14 $x $y 0.05 $fun
#sbatch --mem=30G --time=168:00:00 --array=0-29%15 Ising_parameter_scaling.sh 15 $x $y 0.05 $fun
#sbatch --mem=120G --time=168:00:00 --array=0-29%15 Ising_parameter_scaling.sh 16 $x $y 0.05 $fun
