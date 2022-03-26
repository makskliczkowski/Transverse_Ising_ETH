set autoscale xy
unset logscale
set size square
set xrange[0:0.5]
set ylabel 'P(|r|)'
set xlabel '|r|'
set key vertical inside right top

plot 'rSigmaXDist_L=10,J0=0.20,g=1.00,g0=0.00,h=0.00,w=1.00.dat' w lp pt 6 ps 0.75 title 'L=10',\
'rSigmaXDist_L=11,J0=0.20,g=1.00,g0=0.00,h=0.00,w=1.00.dat' w lp pt 6 ps 0.75 title 'L=11',\
'rSigmaXDist_L=12,J0=0.20,g=1.00,g0=0.00,h=0.00,w=1.00.dat' w lp pt 6 ps 0.75 title 'L=12',\
'rSigmaXDist_L=13,J0=0.20,g=1.00,g0=0.00,h=0.00,w=1.00.dat' w lp pt 6 ps 0.75 title 'L=13'
