set size square
set xrange[-0.6:0.6]
set ylabel '<n|{/Symbol s}_0^x|n>'
set xlabel 'E_n/L'
set key outside above

plot 'SigmaX_L=10,J0=0.20,g=1.00,g0=0.00,h=0.00,w=1.00.dat' w lp pt 6 ps 0.5 title 'L=10',\
'SigmaX_L=11,J0=0.20,g=1.00,g0=0.00,h=0.00,w=1.00.dat' w lp pt 6 ps 0.5 title 'L=11',\
'SigmaX_L=12,J0=0.20,g=1.00,g0=0.00,h=0.00,w=1.00.dat' w lp pt 6 ps 0.5 title 'L=12',\
'SigmaX_L=13,J0=0.20,g=1.00,g0=0.00,h=0.00,w=1.00.dat' w lp pt 6 ps 0.5 title 'L=13'
