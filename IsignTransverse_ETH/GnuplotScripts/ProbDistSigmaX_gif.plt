set autoscale xy
set logscale y
set terminal gif animate delay 100
set output 'ProbDistSigmaX.gif'
set fit quiet

w=0.1
L=12
set size square

set title sprintf("Probability distribution of O = {/Symbol s}_0^x for\nL=%d J=1, J0=0.2, g=1, g0=0, h=0", L)

set key outside center above
set key spacing 4

while(w <= 5.0){
	unset title


f(x) = a/(sigma*sqrt(2.*pi))*exp(-(x-mu)**2./(2.*sigma**2))
name = sprintf("ProbDistSigmaX_L=%d,J0=0.20,g=1.00,g0=0.00,h=0.00,w=%0.2f.dat", L, w)

fit f(x) name u 1:2 via a, mu, sigma

set xrange [-0.6:0.6]
set yrange [0.1:100]
set ylabel '{/*1.3 P(O_{nn} - <O_{nn}>)}' offset 1.2, 0
set xlabel '{/*1.3 O_{nn} - <O_{nn}>}'


plot f(x) w l lw 4 title sprintf("\ngaussian fit:\n a=%0.2f\n sigma=%0.2f\n mu=%0.2f", a, sigma, mu),\
name u 1:2 with lp pt 6 ps 0.6 title sprintf("distribution\n for w = %0.2f", w)
	w = w + 0.1
}