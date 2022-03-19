set autoscale xy

set terminal gif animate delay 100
set output 'ProbDistGap.gif'
set fit quiet

w=0.1
L=12
set size square

set title sprintf("Probability distribution of r for\nL=%d J=1, J0=0.2, g=1, g0=0, h=0", L)

set key outside center above
set key spacing 4

while(w <= 5.0){
	unset title


f(x) = a * (c + x + x**2) / sqrt( (b + x + x**2)**5 )
name = sprintf("ProbDistGap_L=%d,J0=0.20,g=1.00,g0=0.00,h=0.00,w=%0.2f.dat", L, w)

fit f(x) name u 1:2 via a, b, c

set xrange [0:1]

set ylabel '{/*1.3 P(r)}' offset 1.2, 0
set xlabel '{/*1.3 r}'


plot f(x) w l lw 4 title sprintf("\nGOE fit:\n a=%0.2f\n b=%0.2f\n c=%0.2f", a, b, c),\
name u 1:2 with lp pt 6 ps 0.6 title sprintf("distribution\n for w = %0.2f", w)
	w = w + 0.1
}