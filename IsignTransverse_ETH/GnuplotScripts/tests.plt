reset 
##--PREAMBLE
set autoscale
set term qt size 1050, 900 font "Helvetica,14"
set mxtics
set mytics
set style line 12 lc rgb '#ddccdd' lt 1 lw 1.5
set style line 13 lc rgb '#ddccdd' lt 1 lw 0.5
set grid xtics mxtics ytics mytics back ls 12, ls 13
set size square
set xtics nomirror
set key inside bottom right font ",16" spacing 2

set style line 1 dt (3,5,10,5) lc rgb "black" lw 1.5
set style line 2 dt (3,3) lc rgb "red" lw 1.5
set style line 3 dt (8,8) lc rgb "blue" lw 1.5
set style line 4 dt (1,1) lc rgb "green" lw 1.5

fileexist(name)=1#int(system("if exist \"".name."\" ( echo 1) else (echo 0)"))

# Margins for each row resp. column
TMARGIN = "set tmargin at screen 0.99; set bmargin at screen 0.56"
BMARGIN = "set tmargin at screen 0.46; set bmargin at screen 0.10"
LMARGIN = "set lmargin at screen 0.10; set rmargin at screen 0.46"
RMARGIN = "set lmargin at screen 0.46; set rmargin at screen 0.82"
RANGE = "set xrange[0:1]; set yrange[0:2.0]"
UNSET = "unset tics; unset xlabel; unset ylabel; unset title; unset key; unset border;"
#-- PARAMETERS
w = 0.0
g = 0.4
L = 16
h = 0.0
h0 = 80
hend = 400
dh = 40
E_or_O = 0  	#energies(0) or operator(1)
var = 1			# variance?
smoothen = 0	# smooth data?
scaling = 0		# use size scaling?
g_scaling = 1	# scaling with g
y_ax = 2		# =0 - ||O_diag||; =1 - ||O_off||; =2 - <r>

#-- GRAPHICS
my_title = "level statistics for"
if(y_ax < 2) my_title = "operator Hilbert-Schmidt norm ||O||^2"
if(g_scaling && !scaling){
	my_title = my_title.sprintf(" L = %d and h = %.2f", L, h)
} else{
	if(scaling){
		my_title = my_title.sprintf(" g = %.2f and h = %.2f", g, h)
	} else{
		my_title = my_title.sprintf(" L = %d and g = %.2f", L, g)
	}
}
#set title my_title

#set y2tics -100, 10
set ytics nomirror
set xlabel "{/*1.25q/{/Symbol p}}"
set ylabel "{/*1.25{/Symbol \341} A_q\267{/Symbol s}^z_{L/2}{/Symbol \361}}"

#set y2tics 0.37, 0.05
#set y2range[0.37:0.54]
set xrange[0.1:5.0]
set ytics nomirror
set autoscale xy
#set autoscale y2

set lmargin at screen 0.13; set rmargin at screen 0.96;
set tmargin at screen 0.96; set bmargin at screen 0.10
#set logscale x
#plot for[gx=10:50:10] '.\Fidelity\'.sprintf("FidelityDecay_L=10,J0=0.00,g=%.2f,g0=0.00,h=0.80,w=0.01.dat", 0.01*gx) u 1: (sqrt($2*$2+$3*$3)) w l title sprintf("g=%.2f", 0.01*gx)

#plot for[hx=20:120:20] '.\Fidelity\'.sprintf("FidelityDecay_L=10,J0=0.00,g=0.20,g0=0.00,h=%.2f,w=0.01.dat", 0.01*hx) u 1: (sqrt($2*$2+$3*$3)) w l title sprintf("h=%.2f", 0.01*hx)

#set xrange[0:2]
#plot "OperatorProductHq_L=16,J0=0.00,g=0.05,g0=0.00,h=0.80,w=0.01.dat" u ($1*2./18.):(($2)) w lp lc rgb "black" pt 7 ps 1.5 title "A_q=H_q",\
#	 "OperatorProductHq_L=16,J0=0.00,g=0.05,g0=0.00,h=0.80,w=0.01.dat" u ($1*2./18.):(($3)) w lp ls 1 pt 6 ps 1.5 notitle,\
#	 "OperatorProductSigmaZq_L=16,J0=0.00,g=0.05,g0=0.00,h=0.80,w=0.01.dat" u ($1*2./18.):(sqrt($2 * $2 + $3 * $3)) w l lc rgb "red" lw 1 notitle,\
#	 "OperatorProductSigmaZq_L=16,J0=0.00,g=0.05,g0=0.00,h=0.80,w=0.01.dat" u ($1*2./18.):(($2)) w lp lc rgb "red" pt 7 ps 1.5 title "A_q={/Symbol s}_q^z",\
#	 "OperatorProductSigmaZq_L=16,J0=0.00,g=0.05,g0=0.00,h=0.80,w=0.01.dat" u ($1*2./18.):(($3)) w lp ls 2 pt 6 ps 1.5  notitle

i0=20; iend=80; di=20
set xlabel 't / t_H'
#set logscale y
#set format x '10^{%L}'
#set format y '10^{%L}'

dir_base='../results/disorder/PBC/'
dir = dir_base.'Magnetization/'
_name(Lx) = dir.sprintf("IsingQPT_L=%d,J0=0.00,g0=0.00,h=%.2f,w=%.2f.dat", Lx, h, w);
set xrange[0:2]
set yrange[0:1]
plot for [Lx=10:18:2] _name(Lx) using 1:2 w lp title sprintf("L=%d", Lx)


exit;
#_name = dir.sprintf("quench_g_init=0.00,h_init-0.50_L=%d,g=%.2f,h=%.2f,k=0,p=1,x=1.dat", L, g, h)
eh=0
nejm = eh? "quench_ferromagnet" : "quench_domain_wall";
_name = dir.nejm.sprintf("_L=%d,J0=0.00,g=%.2f,g0=0.00,h=%.2f,w=%.2f.dat", L, g, h, w);
#plot for[gx=20:60:10] dir.sprintf("meson_coverage_L=%d,J0=0.00,g=%.2f,g0=0.00,h=%.2f,w=%.2f.dat", L, gx*0.01, h, w) u 1:2 w lp title sprintf("g=%.2f", 0.01*gx)
#exit;
 heatmap = 0
if(heatmap){
	set xrange[-0.5:L-0.5]
	set yrange[0:50]
	set view map
	set pm3d interpolate 0,0
	#set cbrange[0:0.1]
	set yrange
	splot _name u 1:2:3 with image
} else {
	set key outside right
	set xrange[0:60];
	set logscale x
	i = 0;
	plot _name u ($1 == 1? $2 : NaN):3 w lp title 'i=1',\
		_name u ($1 == 2? $2 : NaN):3 w lp title 'i=2',\
		_name u ($1 == L / 2? $2 : NaN):3 w lp title 'i=L/2',\
		_name u ($1 == L / 2 + 1? $2 : NaN):3 w lp title 'i=L/2+1'
}

exit;
#dir = dir_base.'SpectralFormFactor/'
#ssf_name(gx,Lx) = dir.sprintf("_L=%d,J0=0.00,g=%.2f,g0=0.00,h=0.20,w=0.01.dat", Lx, gx)
#key_title(i) = sprintf("g=%.2f", 0.01*i)
#	array tH[(iend-i0)/di+1]
#	array LTA[(iend-i0)/di+1]
#	array val[(iend-i0)/di+1]
#do for[i=i0:iend:di]{
#		_name = ssf_name(0.01*i, L)	
#		if(fileexist(_name)){
#			stats _name every ::0::1 using 2 nooutput;	val[(i-i0)/di+1] = STATS_min
#			stats _name every ::0::1 using 3 nooutput; 	tH[(i-i0)/di+1] = STATS_min;
#			stats _name every ::0::1 using 4 nooutput; 	LTA[(i-i0)/di+1] = STATS_min
#			print key_title(i),"  ", val[(i-i0)/di+1], tH[(i-i0)/di+1], LTA[(i-i0)/di+1]
#		}
#		else{
#			print key_title(i)," -----------------------------------------------------------------------------"
#		}
#}
#plot for[i=i0:iend:di] ssf_name(0.01*i, L) u ($1 / tH[(i-i0)/di + 1]):2 w lp lw 1.5 t sprintf("g=%.2f",0.01*i), (x < 1? 2*x - x*log(1+2*x) : 2-x*log( (2*x+1) / (2*x-1))) t 'GOE'
#
#exit;
#set logscale x; set xrange[*:150]
dir = dir_base.'Entropy/'
set ylabel 'S(t)'
set xlabel 't'
ssf_name(gx,Lx,x,M) = dir.sprintf("TimeEvolution_L=%d,J0=0.00,g=%.2f,g0=0.00,h=0.80,w=0.01,x=%.4f,M=%d.dat", Lx, gx, x, M)
M=6
dt=0.16
L=10
M_list = '2 3 4 5 6 7 8 9'
dt_list = '1e-2 2e-2 4e-2 8e-2 16e-2 64e-2 128e-2 256e-2 512e-2'
MARGIN = "set lmargin at screen 0.10; set rmargin at screen 0.99; set bmargin at screen 0.10; set tmargin at screen 0.99;"
set xrange[1e-3:1e2]; set yrange[1e-6:4.0]
set key right bottom
#plot for[Mx in M_list] ssf_name(0.6, L, dt, 1.*Mx) u 1:3 w lp pt 6 ps 0.5 t sprintf("x=%.3f, M=%d", dt, 1.*Mx), ssf_name(0.6, L, dt, M) u 1:2 w lp ls 1 title 'ED'
plot for[dt in dt_list] ssf_name(0.6, L, 1.*dt, M) u 1:3 w lp pt 6 ps 1.5 t sprintf("x*{/Symbol w}_{max}=%.3f, M=%d", 1.*dt, M), ssf_name(0.6, L, 1e-2, M) u 1:2 w lp ls 1 title 'ED'
