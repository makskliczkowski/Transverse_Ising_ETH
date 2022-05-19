	#R_ARR=(200 200 200 100 100 40 40 20 10);	r=${R_ARR[`expr ${1}-8`]}
	R_ARR=(200 200 200 100 100 40 40 20 10);	r=${R_ARR[`expr ${1}-12`]}
	
# for L >=14 r_max = 50
for realis in {0..50}; do
	id=$(echo $r $realis | awk '{printf "%d", $1 * $2}')
    sh Ising_parameter_scaling_sender.sh 0.9 0.8 3 0 0 sff $id
	#sh Ising_parameter_scaling_sender.sh 0.9 0.8 0 0 0 diag2 $id
done