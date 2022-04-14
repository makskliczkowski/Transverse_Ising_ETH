for h in $(seq 0.05 0.05 1.5)
do
    sh Ising_g_scaling_sender.sh $h 0 2 0 diag 
done
