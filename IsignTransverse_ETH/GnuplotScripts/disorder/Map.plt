reset 
use_png = 1		# 1 if use png output, and 0 for qt output
if(use_png) { set term pngcairo size 1250, 1200 font sprintf("Helvetica,%d",16); }
else {set term qt size 950, 900 font sprintf("Helvetica,%d",14); }

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

MARGIN = "set lmargin at screen 0.10; set rmargin at screen 0.95; set bmargin at screen 0.10; set tmargin at screen 0.95;"
UNSET = "unset tics; unset xlabel; unset ylabel; unset title; unset border;"

    glist2 = '0.025 0.05 0.075 0.10 0.125 0.2 0.3 0.35 0.4 0.45 0.50 0.55 0.60 0.65 0.70 0.75'
	glist = '0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5'# 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0'
	h0 = 20;	hend = 180; 	dh = 20;
heatmap = 0 		# ==1 plot 2D heatmap, else plot cuts at specific values
h_vs_g = 0;			# ==0 --> as function of g
relax_vs_th = 1		# pick relaxation-time=1 ot thouless-time=0
idx = relax_vs_th? 3 : 5 
scaling = 2			# = 0 q/j-scaling, =1-size scaling, =2-h/g


q_vs_j = 0
site = 0
operator = 1	 		# 1-SigmaZ , 0-Hq :local
rescale_times = 0
L = 14
g = 0.1
h=0.8

set lmargin at screen 0.15
set rmargin at screen 0.85
set bmargin at screen 0.15
set tmargin at screen 0.9

dir_base='../../results/disorder/PBC/'
dir = dir_base.'RelaxationTimes/'
load dir_base.'gnuplot-colorbrewer-master/diverging/RdYlGn.plt'

out_dir = "Relaxation_Times/".(h_vs_g? "vs_h/" : "vs_g/")
command = "mkdir ".out_dir; 
system command

str = (h_vs_g? "_h" : "_g")
op = operator? "SigmaZ" : "H";
_name(Lx, s) = str.op."_".(q_vs_j? "q" : "j").sprintf("=%d_L=%d,J0=0.00,g0=0.00,w=0.01", s, Lx);
	
nu=2
rescale(t,L) = rescale_times? t/L**nu : t

#do for[gx in glist]{
#	g = 1. * gx
#do for[site=-1:7]{
#do for[L=10:15]{
#do for[hx=h0:hend:dh]{
#h = 0.01*hx

set encoding utf8

	my_title = "Relaxation time for"
	if(scaling != 1) { my_title = my_title.sprintf(" L=%d", L); }
	if(scaling != 2) { my_title = my_title.(h_vs_g? sprintf(" g=%.3f", g) : sprintf(" h=%.2f", h)); }
	tmp_title = my_title."\n\n for operator A = ";
	if(operator){ my_title = tmp_title.(q_vs_j? "L^{-1/2}{/Symbol S}_j e^{iqj} {/Symbol s}^z_j" : "{/Symbol s}^z_j");
	} else{ my_title = tmp_title.(q_vs_j? "L^{-1/2}{/Symbol S}_j cos(2{/Symbol p}/L\267qj) H^j" : "H^j");}
	
	str(s) = (q_vs_j? "q" : "j").(s<0? "=L/2" : sprintf("=%d",s))
	if(scaling != 0){ my_title = my_title."; ".(str(site)); }
	#set title my_title
	if(use_png){
		output_name = out_dir.op."_".(q_vs_j? "q" : "j");
		if(scaling != 0) { output_name = output_name.sprintf("=%d", site); }
		if(scaling != 1) { output_name = output_name.sprintf("_L=%d", L); }
		output_name = output_name.",J0=0.00,";
		if(scaling != 2){ output_name = output_name.(h_vs_g? sprintf("g=%.2f,g0=0.00,", g) : sprintf("g0=0.00,h=%.2f,", h))}
		set output output_name."w=0.01.png";
	}
	name = _name(L,site)
	
	name = dir.name.".dat"
	
	set key outside right top
	nejm = relax_vs_th? 'rel' : 'Th'
	label_y=rescale_times? sprintf("{/*1.5t_{%s}/L^{%d}}", nejm, nu) : '{/*1.5t_{'.nejm.'}}'
	set ylabel label_y #rotate by 0 offset 2,0.
	#set ylabel '{/*1.5t_{rel}/t_H}' rotate by 0 offset 2,0.
	set arrow from 0,1 to 2.5,1 nohead
	if(h_vs_g){
		set xrange[0:2.5]
		#set yrange[0.01:1000]
		set xlabel "{/*1.5h/J}"
		if(scaling == 2){
			plot for[gx in glist] name u ($2 == gx + 0.0? $1 : NaN):($3/$4) w lp ls ((1.*gx - 0.05)/0.05) pt 6 ps 1.5 title "g=".gx
		} else{
			if(scaling == 0){
				plot for[i=0:(q_vs_j? L/2 : L-1)] dir._name(L,i).".dat" u ($2 == g? $1 : NaN):($3/$4) w lp ls (i+1) pt (i+3) ps 1.5 title str(i)
			} else{
				plot for[i=10:15] dir._name(i,site<0? i / 2 : site).".dat" u ($2 == g? $1 : NaN):($3/$4) w lp ls ((i-7)) pt (i-5) ps 1.5 title sprintf("L=%d", i)
			}
		}
	} else{
		set xrange[0:1];
		RANGE="set yrange[1e-1:1e5]; set xrange[0.1:1.5];"
		set xlabel "{/*1.5g/J}"
		if(scaling == 2){
			plot for[i=h0:hend:dh] name u (100*$1 == i? $2 : NaN):($3/$4) w lp ls ((i-h0)/dh) pt 6 ps 1.5 title sprintf("h=%.2f", 0.01*i)
		} else{
			set multiplot
			@MARGIN; @FORMAT; @RANGE;
			if(scaling == 0){
				plot for[i=0:(q_vs_j? L/2 : L-1)] dir._name(L,i).".dat" u ($1 == h? $2 : NaN):(column(idx)) w lp ls (i+1) pt (i+4) ps 1.5 title str(i)
				@UNSET; @MARGIN; @RANGE; plot dir._name(L,1).".dat" u ($1 == h? $2 : NaN):4 w l ls 0 lw 3 notitle
			} else{
				L_list = '10 11 12 13 15'; set logscale x
				plot for[i in L_list] dir._name(1.0*i,site<0? i / 2 : site).".dat" u ($1 == h? $2 : NaN):(rescale((column(idx)),1.0*i)) w lp ls ((1.0*i-7)) pt (2*(1.0*i-10)+5) ps 1 title sprintf("L=%d", 1.0*i) #, 1e2*x**(-2) notitle
				if(0&&!rescale_times){
					@UNSET; @MARGIN;
					plot for[i in L_list] dir._name(1.0*i,site<0? i / 2 : site).".dat" u ($1 == h? $2 : NaN):4 w lp ls 0 lw 2 pt (2*(1.0*i-10)+4) ps 1.5 notitle
				}
			}
			unset multiplot
		}	
	}
#}