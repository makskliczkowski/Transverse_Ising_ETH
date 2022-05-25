reset 
use_png = 0		# 1 if use png output, and 0 for qt output
if(use_png) { set term pngcairo size 1200, 1200 font sprintf("Helvetica,%d",18); }
else {set term qt size 950, 900 font sprintf("Helvetica,%d",18); }
load './gnuplot-colorbrewer-master/diverging/RdYlGn.plt'

FORMAT = "set style line 12 lc rgb '#ddccdd' lt 1 lw 1.5; \
set style line 13 lc rgb '#ddccdd' lt 1 lw 0.5; \
set grid xtics mxtics ytics mytics back ls 12, ls 13;\
set size square;\
set xtics mirror;\
set ytics mirror; set border;"
@FORMAT
set autoscale
set logscale y
set format y "10^{%L}"
set size square
#set xrange[0.05:5]

MARGIN = "set lmargin at screen 0.10; set rmargin at screen 0.99; set bmargin at screen 0.10; set tmargin at screen 0.99;"
UNSET = "unset tics; unset xlabel; unset ylabel; unset title; unset border;"

    glist2 = '0.025 0.05 0.075 0.10 0.125 0.2 0.3 0.35 0.4 0.45 0.50 0.55 0.60 0.65 0.70 0.75'
	glist = '0.5 0.6 0.7 0.8 1.5'# 0.3 0.35 0.4 0.45 0.5'# 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0'
	h0 = 20;	hend = 180; 	dh = 20;
heatmap = 0 		# ==1 plot 2D heatmap, else plot cuts at specific values
h_vs_g = 0;			# ==0 --> as function of h on x-axis
relax_vs_approx = 0	# pick initial relax time = 1 or approx with renormalized peak = 0
plot_thouless = 0	# plot only thouless times
scaling = 3			# = 0 q/j-scaling, =1-size scaling, =2-h/g, =3-compare_operators
user_defined = 1

q_vs_j = 1
site = 1
q = 5
operator = 1	 		# 1-SigmaZ , 0-Hq :local
rescale_times = 0
L = 12
g = 0.8
h=0.8
J=1.0
w=0.01

fileexist(name)=system("[ -f '".name."' ] && echo '1' || echo '0'") + 0
set lmargin at screen 0.15
set rmargin at screen 0.85
set bmargin at screen 0.15
set tmargin at screen 0.9

dir_base='../results/disorder/PBC/'
dir = dir_base.'RelaxationTimes/'

str = (h_vs_g? "_h" : "_g")