set autoscale xy
set size square
set logscale y
set ylabel 'z_n'
set xlabel 'L - system size'
set key outside above

plot 'SpectrumRapScalingSigmaX,J0=0.20,g=1.00,g0=0.00,h=0.00,w=1.00.dat' u (log10($1))/log10(2):2 w lp pt 6 ps 1 title 'average',\
'SpectrumRapScalingSigmaX,J0=0.20,g=1.00,g0=0.00,h=0.00,w=1.00.dat' u (log10($1))/log10(2):3 w lp pt 6 ps 1 title '1^{st} outlier',\
'SpectrumRapScalingSigmaX,J0=0.20,g=1.00,g0=0.00,h=0.00,w=1.00.dat' u (log10($1))/log10(2):4 w lp pt 6 ps 1 title '2^{nd} outlier',\
'SpectrumRapScalingSigmaX,J0=0.20,g=1.00,g0=0.00,h=0.00,w=1.00.dat' u (log10($1))/log10(2):6 w lp pt 6 ps 1 title '4^{th} outlier',\
'SpectrumRapScalingSigmaX,J0=0.20,g=1.00,g0=0.00,h=0.00,w=1.00.dat' u (log10($1))/log10(2):8 w lp pt 6 ps 1 title '6^{th} outlier'
