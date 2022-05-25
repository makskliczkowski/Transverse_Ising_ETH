set autoscale xy
f(x) = a/(sigma*sqrt(2.*pi))*exp(-(x-mu)**2./(2.*sigma**2))
name = 'ProbDist_sigma_x_L=10,J0=0.20,g=1.00,g0=0.00,h=0.00,w=1.00.dat'

fit f(x) name u 1:2 via a, mu, sigma

set xrange[-1:1]
set ylabel 'P({/Symbol s}_{nn} - <{/Symbol s}_{nn}>)'
set xlabel '{/Symbol s}_{nn} - <{/Symbol s}_{nn}>'
set title 'L=10, J=1.0 ± 0.2, g=1.0 ± 0.0, h=0.0 ± 1.00' 

set key inside top right
set key spacing 4

plot f(x) w l lw 4 title sprintf("\ngaussian fit:\n a=%0.2f\n sigma=%0.2f\n mu=%0.2f", a, sigma, mu),\
name u 1:2 with boxes title sprintf("distribution")