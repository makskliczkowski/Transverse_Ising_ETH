for h in $(seq 0.2 0.2 3.0)
do
    sh Ising_g_scaling_sender.sh $h 0 2 0 ${h}_spectrals    
done
