f(x) = 27 / 4 * (x+x**2) / sqrt( (1 + x + x**2)**5 )

set ylabel 'P(r)'
set xlabel 'r'
set key outside


plot f(x) title 'GOE - get get', 'ProbDistGap_L=13,J0=0.20,g=1.00,g0=0.00,h=0.00,w=1.00.dat' u 1:(($2/8.3)) w boxes title 'boom boom boom'