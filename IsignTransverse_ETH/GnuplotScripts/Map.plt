reset 
use_png = 0		# 1 if use png output, and 0 for qt output
if(use_png) { set term pngcairo size 1250, 1200 font sprintf("Helvetica,%d",18); }
else {set term qt size 1050, 900 font sprintf("Helvetica,%d",18); }

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
	glist = '0.5 0.6 0.7 0.8 1.5'# 0.3 0.35 0.4 0.45 0.5'# 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0'
	h0 = 20;	hend = 180; 	dh = 20;
heatmap = 1 		# ==1 plot 2D heatmap, else plot cuts at specific values
h_vs_g = 1;			# ==0 --> as function of h on x-axis
relax_vs_th = 1		# pick relaxation-time=1 ot thouless-time=0
idx = relax_vs_th? 3 : 5 
scaling = 0		# = 0 q/j-scaling, =1-size scaling, =2-h/g, =3-compare_operators
user_defined = 0

q_vs_j = 1
site = 1
q = 1
operator = 2	 		# 1-SigmaZ , 0-Hq :local
rescale_times = 0
L = 14
g = 0.8
h=0.8

fileexist(name)=system("[ -f '".name."' ] && echo '1' || echo '0'") + 0
set lmargin at screen 0.15
set rmargin at screen 0.85
set bmargin at screen 0.15
set tmargin at screen 0.9

dir_base='../results/disorder/PBC/'
dir = dir_base.'RelaxationTimes/'
load './gnuplot-colorbrewer-master/diverging/RdYlGn.plt'

str = (h_vs_g? "_h" : "_g")
op = ""
if(operator == 0) {op = "H"; }
if(operator == 1) {op = "SigmaZ";}
if(operator == 2) {op = "TFIM_LIOM_plus";}
if(operator == 3) {op = "TFIM_LIOM_minus";}
_str(x) = (q_vs_j? "q" : "j").sprintf("=%d",x);
if(operator > 1){ _str(x) = "n".sprintf("=%d",x); }

_base(Lx) = sprintf("_L=%d,J0=0.00,g0=0.00,w=0.01.dat", Lx)
_name(Lx, s) = str.op."_"._str(s)._base(Lx);

nu=2
rescale(t,L) = rescale_times? t/L**nu : t
use_fit = 0
#do for[gx in glist]{
#	g = 1. * gx
#do for[site=-1:7]{
#do for[L=10:15]{
#do for[hx=h0:hend:dh]{
#h = 0.01*hx

set encoding utf8
if(user_defined == 0)
{
	if(use_png){
		out_dir = "Relaxation_Times/".(h_vs_g? "vs_h/" : "vs_g/")
		command = "mkdir ".out_dir; 
		system command
		output_name = out_dir.op."_".(q_vs_j? "q" : "j");
		if(scaling != 0) { output_name = output_name.sprintf("=%d", q_vs_j? q : site); }
		if(scaling != 1) { output_name = output_name.sprintf("_L=%d", L); }
		output_name = output_name.",J0=0.00,";
		if(scaling != 2){ output_name = output_name.(h_vs_g? sprintf("g=%.2f,g0=0.00,", g) : sprintf("g0=0.00,h=%.2f,", h))}
		set output output_name."w=0.01.png";
	}
	name = _name(L,( q_vs_j? q : site))
	name = dir.name
	
	set key inside right top
	nejm = relax_vs_th? 'rel' : 'Th'
	label_y=rescale_times? sprintf("{/*1.5t_{%s}/L^{%d}}", nejm, nu) : '{/*1.5t_{'.nejm.'}}'
	set ylabel label_y #rotate by 0 offset 2,0.
	#set ylabel '{/*1.5t_{rel}/t_H}' rotate by 0 offset 2,0.
	#set arrow from 0,1 to 2.5,1 nohead
	alfa = -6
	f(x) = 7e5 * (x/0.1)**(alfa)
	if(h_vs_g){
		set xrange[0.01:3.0]
		set logscale x
		#set yrange[0.01:1000]
		set xlabel "{/*1.5h/J}"
		if(scaling == 0){ 	plot for[i=0:(q_vs_j? L/2 : L-1)] dir._name(L,i) u ($2 == g? $1 : NaN):($3) w lp ls (i+1) pt (i+3) ps 1.5 title _str(i)
		} else {
			if(scaling == 1){ plot for[i=L0:Lend:dL] dir._name(i, (q_vs_j? (q<0? i / 2 : q) : site)) u ($2 == g? $1 : NaN):($3/$4) w lp ls ((i-7)) pt (i-5) ps 1.5 title sprintf("L=%d", i)
			} else {
				if(scaling == 2){ plot for[gx in glist] name u ($2 == gx + 0.0? $1 : NaN):($3/$4) w lp ls ((1.*gx - 0.05)/0.05) pt 6 ps 1.5 title "g=".gx
				} else {
					if(scaling == 3){ 
						plot dir."_hSigmaZ_j=".sprintf("%d", site)._base(L) u ($2 == g? $1 : NaN):($3) w lp pt 6 ps 1.5 title sprintf("{/Symbol s}_{j=%d}}", site),\
							dir."_hSigmaZ_q=".sprintf("%d", 1)._base(L) u ($2 == g? $1 : NaN):($3) w lp pt 6 ps 1.5 title sprintf("{/Symbol s}_{q=%d}}", 1),\
							dir."_hSigmaZ_q=".sprintf("%d", L / 2.)._base(L) u ($2 == g? $1 : NaN):($3) w lp pt 6 ps 1.5 title sprintf("{/Symbol s}_{q=%d}}", L / 2.),\
							dir;"_hH_j=".sprintf("%d", site)._base(L) u ($2 == g? $1 : NaN):($3) w lp pt 4 ps 1.5 title sprintf("H_{j=%d}}", site),\
							dir;"_hH_q=".sprintf("%d", 1)._base(L) u ($2 == g? $1 : NaN):($3) w lp pt 4 ps 1.5 title sprintf("H_{q=%d}}", 1),\
							dir;"_hH_q=".sprintf("%d", L / 2.)._base(L) u ($2 == g? $1 : NaN):($3) w lp pt 4 ps 1.5 title sprintf("H_{q=%d}}", L / 2.)
					} else{ 
						plot name u ($2 == g? $1 : NaN):($3/$4) w lp pt 6 ps 1.5 title sprintf("L=%d,g=%.2f,h=%.2f",L,g,h)
		}}}}
	} else{
		set xrange[0:1.5];
		set xlabel "{/*1.5g/J}"
		#set logscale x
		#unset logscale y; set format y '%g'
		#set yrange[-3:10];
		if(scaling == 0){	plot for[i=0:(q_vs_j? L/2 : L-1)] dir._name(L,i) u ($1 == h? $2 : NaN):((log($3)*$2**0.0)) w lp ls (i+1) pt (i+4) ps 1.5 title _str(i),\
								dir._name(L,1) u ($1 == h? $2 : NaN):4 w l ls 0 lw 3 notitle#, f(x) w l ls 0 lw 4 lc rgb 'blue' notitle,
		} else {
			if(scaling == 1){ plot for[i=L0:Lend:dL] dir._name(i, (q_vs_j? (q<0? i / 2. : q) : site)) u ($1 == h? $2 : NaN):(rescale(($3),i)) w lp ls ((i+3-L0)) pt ((i-L0)/dL+1) ps 1 title sprintf("L=%d", i)
			} else {
				if(scaling == 2){ plot for[i=h0:hend:dh] name u (100*$1 == i? $2 : NaN):($3/$4) w lp ls ((i-h0)/dh) pt 6 ps 1.5 title sprintf("h=%.2f", 0.01*i)
				} else {
					if(scaling == 3){ 
						set key spacing 2
						plot dir."_gSigmaZ_j=".sprintf("%d", site)._base(L) u ($1 == h? $2 : NaN):($3) w lp pt 6 ps 1.5 title sprintf("{/Symbol s}_{j=%d}", site),\
							dir."_gSigmaZ_q=".sprintf("%d", 1)._base(L) u ($1 == h? $2 : NaN):($3) w lp pt 6 ps 1.5 title sprintf("{/Symbol s}_{q=%d}", 1),\
							dir."_gSigmaZ_q=".sprintf("%d", L / 2.)._base(L) u ($1 == h? $2 : NaN):($3) w lp pt 6 ps 1.5 title sprintf("{/Symbol s}_{q=%d}", L / 2.),\
							dir."_gH_j=".sprintf("%d", site)._base(L) u ($1 == h? $2 : NaN):($3) w lp pt 4 ps 1.5 title sprintf("H_{j=%d}", site),\
							dir."_gH_q=".sprintf("%d", 1)._base(L) u ($1 == h? $2 : NaN):($3) w lp pt 4 ps 1.5 title sprintf("H_{q=%d}", 1),\
							dir."_gH_q=".sprintf("%d", L / 2.)._base(L) u ($1 == h? $2 : NaN):($3) w lp pt 4 ps 1.5 title sprintf("H_{q=%d}", L / 2.)
					} else{ 
						plot name u ($1 == h? $2 : NaN):($3/$4) w lp pt 6 ps 1.5 title sprintf("L=%d,g=%.2f,h=%.2f",L,g,h)
		}}}}}



} else {

set style line 1 dt (3,5,10,5) lc rgb "violet" lw 1.5
set style line 2 dt (3,3) lc rgb "red" lw 1.5
set style line 3 dt (8,8) lc rgb "blue" lw 1.5
set style line 4 dt (3,3) lc rgb "green" lw 1.5
set style line 5 dt (8,8) lc rgb "orange" lw 1.5
set style line 10 dt (8,8)

	info(s, gx, hx) = "_".(q_vs_j? "q" : "j").sprintf("=%d,J0=0.00,g=%.2f,g0=0.00,h=%.2f,w=0.01", s, gx, hx);
	_name(s, gx, hx) = "_L".op.info(s, gx, hx)
	set key right outside spacing 2
	RANGE="set yrange[5e-1:4e3]; set xrange[1./16.:0.1];"
	set ylabel '{/Symbol t}_{rel}' rotate by 0 font ",22"
	set xlabel '1 / L' font ",22"
	set xtics add("1/10" 1./10.); set xtics add("1/12" 1./12.); set xtics add("1/14" 1./14.); set xtics add("1/16" 1./16.);
	if(scaling == 2){ 
		i0 = 30; di = 10; iend = 120;

		a=1; b=1;
		f(x) = a * x**b
		f_plot(as, bs, x) = as * x**bs
		size = (iend - i0) / di+1
		array a_list[size]
		array b_list[size]

		do for[i=i0:iend:di]{
			idx = (i-i0)/di + 1
			if(use_fit){ fit f(x) dir._name(site, 0.01*i, h).".dat" u 1:2 via a, b; }
			a_list[idx] = a; b_list[idx] = b;
			print a, b
		}
		set multiplot
			@MARGIN; @RANGE;
			plot for[i=i0:iend:di] dir._name(site, 0.01*i, h).".dat" u (1./$1):2 w p pt 7 ps 2 lc ((i-i0)/di+1) title sprintf("g=%.2f", 0.01*i)
			if(use_fit){
				@UNSET; @MARGIN; @RANGE; set key outside bottom
				plot for[i=1:size] f_plot(a_list[i], b_list[i], 1./x) w l ls 10 lc i lw 2 t (abs(b_list[i]) < 0.2? "O(1)" : sprintf("O( L^{%.2f} )", b_list[i]))
			}
		unset multiplot
	} else {
		if(scaling == 0){ 
			plot for[i=-1:1:1] dir._name(i, g, h).".dat" u 1:2 w lp ls (i+1) pt (i+4) ps 1.5 title (i < 0? "q=L/2" : sprintf("q=%d", i))
		} else{
			@MARGIN; @RANGE;
			plot dir."_LSigmaZ".info(-1, g, h).".dat" u 1:2 w p pt 7 ps 2 lc 1 title "{/Symbol s}_{q=L/2}", 2.5 w l ls 1 lw 2 t "O(1)",\
				dir."_LSigmaZ".info(1, g, h).".dat" u 1:2 w p pt 7 ps 2 lc 2 title "{/Symbol s}_{q=1}", 0.7*x w l ls 4 lw 2 t "O(L)",\
				dir."_LH".info(-1, g, h).".dat" u 1:2 w p pt 7 ps 2 lc 7 title "H_{q=L/2}", 0.7 w l ls 2 lw 2 t "O(1)",\
				dir."_LH".info(1, g, h).".dat" u 1:2 w p pt 7 ps 2 lc 4 title "H_{q=1}", 0.63*x**1.5 w l ls 5 lw 2 t "O( L^{3/2} )",\
				dir."_LH".info(2, g, h).".dat" u 1:2 w p pt 7 ps 2 lc 3 title "H_{q=2}", 6*(x/10)**1.5 w l ls 3 lw 2 t "O( L^{3/2} )"
			#0.9
			#plot dir."_LSigmaZ".info(-1, g, h).".dat" u 1:2 w p pt 7 ps 2 lc 1 title "{/Symbol s}_{q=L/2}", 0.98 w l ls 1 lw 2 t "O(1)",\
			#	dir."_LSigmaZ".info(1, g, h).".dat" u 1:2 w p pt 7 ps 2 lc 2 title "{/Symbol s}_{q=1}", 1.1*(x/10)**0 w l ls 4 lw 2 t "O(1)",\
			#	dir."_LH".info(-1, g, h).".dat" u 1:2 w p pt 7 ps 2 lc 7 title "H_{q=L/2}", 0.5 w l ls 2 lw 2 t "O(1)",\
			#	dir."_LH".info(1, g, h).".dat" u 1:2 w p pt 7 ps 2 lc 4 title "H_{q=1}", 7.2*(x/10)**2 w l ls 5 lw 2 t "O(L^2)",\
			#	dir."_LH".info(2, g, h).".dat" u 1:2 w p pt 7 ps 2 lc 3 title "H_{q=2}", 2.1*(x/10)**2 w l ls 3 lw 2 t "O(L^2)"
		}
	}
}	
#}