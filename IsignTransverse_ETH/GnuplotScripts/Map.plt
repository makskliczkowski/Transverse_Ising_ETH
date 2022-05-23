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
	w_num = 5;	array w_list[w_num];
	w_list[1] = 0.01;	w_list[2] = 0.05;	w_list[3] = 0.1;	w_list[4] = 0.3;	w_list[5] = 0.5;

	h0 = 20;	hend = 180; 	dh = 20;
	g0 = 10;	gend = 70; 	dg = 10;
    J0 = 35;    Jend = 95;     dJ = 5
heatmap = 1 		# ==1 plot 2D heatmap, else plot cuts at specific values
h_vs_g = 0;			# ==0 --> as function of h on x-axis
relax_vs_approx = 0	# pick initial relax time = 1 or approx with renormalized peak = 0
plot_thouless = 0	# plot only thouless times
scaling = 0			# = 0 q/j-scaling, =1-size scaling, =2-h/g, =3-compare_operators
user_defined = 0

q_vs_j = 1
site = 1
q = 1
operator = 1	 		# 1-SigmaZ , 0-Hq :local
rescale_times = 0
rescale_x_axis = 0
exponent = 1

L = 13
g = 0.9
h=0.8
J=1.0
w = 0.01
w2 = 0.1

fileexist(name)=system("[ -f '".name."' ] && echo '1' || echo '0'") + 0
set lmargin at screen 0.15
set rmargin at screen 0.85
set bmargin at screen 0.15
set tmargin at screen 0.9

dir_base='../results/disorder/PBC/'
dir = dir_base.'RelaxationTimes/'

str = (h_vs_g? "_h" : "_g")
op = ""
if(operator == 0) {op = "H"; }
if(operator == 1) {op = "SigmaZ";}
if(operator == 2) {op = "TFIM_LIOM_plus";}
if(operator == 3) {op = "TFIM_LIOM_minus";}
_str(x) = (q_vs_j? "q" : "j").sprintf("=%d",x);
if(operator > 1){ _str(x) = "n".sprintf("=%d",x); }

_base(Jx, Lx, dis) = sprintf("_L=%d,J=%.2f,J0=0.00,g0=0.00,w=%.2f.dat", Lx, Jx, dis)
_name(Jx, Lx, s) = str.(!plot_thouless? op."_"._str(s) : "")._base(Jx, Lx, w);

_name_th(Jx, Lx, dis) = dir_base.'ThoulessTime/'.str._base(Jx, Lx, dis);
_name_th_L(Jx, gx, hx, dis) = dir_base.'ThoulessTime/'.sprintf("_L,J=%.2f,J0=0.00,g=%.2f,g0=0.00,h=%.2f,w=%.2f.dat", Jx, gx, hx, dis);

nu=2
rescale(t,L) = rescale_times? t/L**nu : t
rescaleX(x, i) = rescale_x_axis? 1./x**exponent : x;
use_fit = 0

set encoding utf8
if(user_defined == 0){
	if(use_png){
		out_dir = (!plot_thouless? 'Relaxation_Times/' : 'Thouless_Times/').(h_vs_g? "vs_h/" : "vs_g/")
		command = "mkdir ".out_dir; 
		system command
		output_name = out_dir.op."_".(q_vs_j? "q" : "j");
		if(scaling != 0) { output_name = output_name.sprintf("=%d", q_vs_j? q : site); }
		if(scaling != 1) { output_name = output_name.sprintf("_L=%d", L); }
		output_name = output_name.",J0=0.00,";
		if(scaling != 2){ output_name = output_name.(h_vs_g? sprintf("g=%.2f,g0=0.00,", g) : sprintf("g0=0.00,h=%.2f,", h))}
		set output output_name."w=0.01.png";
	}
	name = _name(J, L,( q_vs_j? q : site))
	name = dir.name
	
	set key inside right top
	nejm = !plot_thouless? 'rel' : 'Th'
	label_y=rescale_times? sprintf("{/*1.5t_{%s}/L^{%d}}", nejm, nu) : '{/*1.5t_{'.nejm.'}}'
	set ylabel label_y #rotate by 0 offset 2,0.
	#set ylabel '{/*1.5t_{rel}/t_H}' rotate by 0 offset 2,0.
	#set arrow from 0,1 to 2.5,1 nohead
	alfa = -6
		f(y1, y2) = relax_vs_approx? y1 : y2;
	if(h_vs_g){
		set xrange[0.01:1.5]
		#set logscale x
		#set yrange[0.01:1000]
		set xlabel (rescale_x_axis? sprintf("{/*1.5 1/h^{%d}}", exponent) : "{/*1.5h}")
		if(scaling == 0){ 	plot for[i=1:(q_vs_j? L/2 : L-1)] dir._name(J, L,i) u ($2 == g? $1 : NaN):(f($3,$5)) w lp ls (i+1) pt (i+3) ps 1.5 title _str(i)
		} else {
		if(scaling == 1){ plot for[i=L0:15:dL] dir._name(J, i, (q_vs_j? (q<0? i / 2 : q) : site)) u ($2 == g? $1 : NaN):(f($3,$5)) w lp ls ((i-7)) pt (i-5) ps 1.5 title sprintf("L=%d", i)
		} else {
		if(scaling == 2){ plot for[gx in glist] name u ($2 == gx + 0.0? $1 : NaN):(f($3,$5)) w lp ls ((1.*gx - 0.05)/0.05) pt 6 ps 1.5 title "g=".gx
		} else {
		if(scaling == 3){ 
			plot dir."_hSigmaZ_j=".sprintf("%d", site)._base(J, L, w) u ($2 == g? $1 : NaN):(f($3,$5)/$4) w lp pt 6 ps 1.5 title sprintf("{/Symbol s}_{j=%d}}", site),\
				dir."_hSigmaZ_q=".sprintf("%d", 1)._base(J, L, w) u ($2 == g? $1 : NaN):(f($3,$5)/$4) w lp pt 6 ps 1.5 title sprintf("{/Symbol s}_{q=%d}}", 1),\
				dir."_hSigmaZ_q=".sprintf("%d", L / 2.)._base(J, L, w) u ($2 == g? $1 : NaN):(f($3,$5)/$4) w lp pt 6 ps 1.5 title sprintf("{/Symbol s}_{q=%d}}", L / 2.),\
				dir."_hH_j=".sprintf("%d", site)._base(J, L, w) u ($2 == g? $1 : NaN):(f($3,$5)/$4) w lp pt 4 ps 1.5 title sprintf("H_{j=%d}}", site),\
				dir."_hH_q=".sprintf("%d", 1)._base(J, L, w) u ($2 == g? $1 : NaN):(f($3,$5)/$4) w lp pt 4 ps 1.5 title sprintf("H_{q=%d}}", 1),\
				dir."_hH_q=".sprintf("%d", L / 2.)._base(J, L, w) u ($2 == g? $1 : NaN):(f($3,$5)/$4) w lp pt 4 ps 1.5 title sprintf("H_{q=%d}}", L / 2.),\
				_name_th(J, L, w2) u ($2 == g? $1 : NaN):($3) w lp pt 4 ps 1.5 title "{/Symbol t}_{Th}"
		} else{ 
			plot name u ($2 == g? $1 : NaN):(f($3,$5)/$4) w lp pt 6 ps 1.5 title sprintf("L=%d,g=%.2f,h=%.2f",L,g,h)
		}}}}
	} else{
		set xrange[rescaleX(0.05,0):rescaleX(1.5,0)];
		set xlabel (rescale_x_axis? sprintf("{/*1.5 1/g^{%d}}", exponent) : "{/*1.5g}")
		#set logscale x
		#unset logscale y; set format y '%g'
		set yrange[0.1:*];
		ef(x)=0.9/x**2.
		#ef(x) = 1000*exp(-9*x)
		if(scaling == 0){	plot for[i=1:(q_vs_j? L/2 : L-1)] dir._name(J, L,i) u ($1 == h? rescaleX($2, i) : NaN):(f($3,$6)) w lp ls (i+1) pt (i+4) ps 1.5 title _str(i),\
								dir._name(J, L,1) u ($1 == h? rescaleX($2, 1) : NaN):4 w l ls 0 lw 3 notitle,\
							_name_th(J, L, w2) u ($1 == h? rescaleX($2, 0) : NaN):($3*$4) w lp pt 4 ps 1.5 title "{/Symbol t}_{Th}",\
							ef(x) w l dt (3,5,10,5) lc rgb "black" lw 1.5 notitle,\
							1000*exp(-8.5*x) w l dt (3,5,10,5) lc rgb "blue" lw 1.5 notitle
		} else {
		if(scaling == 1){ plot for[i=L0:Lend:dL] dir._name(J, i, (q_vs_j? (q<0? i / 2. : q) : site)) u ($1 == h? rescaleX($2, i) : NaN):(rescale((f($3,$5)),i)) w lp ls ((i+3-L0)) pt ((i-L0)/dL+1) ps 1 title sprintf("L=%d", i),\
							_name_th(J, L, w2) u ($1 == h? rescaleX($2, L) : NaN):($3*$4) w lp pt 4 ps 1.5 title "{/Symbol t}_{Th}"
		} else {
		if(scaling == 2){ plot for[i=h0:hend:dh] name u (100*$1 == i? rescaleX($2, i) : NaN):(f($3,$5)) w lp ls ((i-h0)/dh) pt 6 ps 1.5 title sprintf("h=%.2f", 0.01*i),\
							_name_th(J, L, w2) u ($1 == 0.01*h? rescaleX($2, hend) : NaN):($3*$4) w lp pt 4 ps 1.5 title "{/Symbol t}_{Th}"
		} else {
		if(scaling == 3){ 
			set key spacing 2
			plot dir."_gSigmaZ_j=".sprintf("%d", site)._base(J, L, w) u ($1 == h? $2 : NaN):(f($3,$5)/$4) w lp pt 6 ps 1.5 title sprintf("{/Symbol s}_{j=%d}", site),\
				dir."_gSigmaZ_q=".sprintf("%d", 1)._base(J, L, w) u ($1 == h? $2 : NaN):(f($3,$5)/$4) w lp pt 6 ps 1.5 title sprintf("{/Symbol s}_{q=%d}", 1),\
				dir."_gSigmaZ_q=".sprintf("%d", L / 2.)._base(J, L, w) u ($1 == h? $2 : NaN):(f($3,$5)/$4) w lp pt 6 ps 1.5 title sprintf("{/Symbol s}_{q=%d}", L / 2.),\
				dir."_gH_j=".sprintf("%d", site)._base(J, L, w) u ($1 == h? $2 : NaN):(f($3,$5)/$4) w lp pt 4 ps 1.5 title sprintf("H_{j=%d}", site),\
				dir."_gH_q=".sprintf("%d", 1)._base(J, L, w) u ($1 == h? $2 : NaN):(f($3,$5)/$4) w lp pt 4 ps 1.5 title sprintf("H_{q=%d}", 1),\
				dir."_gH_q=".sprintf("%d", L / 2.)._base(J, L, w) u ($1 == h? $2 : NaN):(f($3,$5)/$4) w lp pt 4 ps 1.5 title sprintf("H_{q=%d}", L / 2.),\
				_name_th(J, L, w2) u ($1 == h? $2 : NaN):($3) w lp pt 4 ps 1.5 title "{/Symbol t}_{Th}"
			} else{ 
				plot name u ($1 == h? $2 : NaN):(f($3,$5)/$4) w lp pt 6 ps 1.5 title sprintf("L=%d,g=%.2f,h=%.2f",L,g,h)
		}}}}}



} else {
	set key right top
	SCALE="set xrange[0.3:1.0]; set yrange[0.01:0.4];"
	@SCALE;
	#unset logscale y; set format y '%g';
	set logscale x
rescale_thouless = 0
conductance = 0
if(conductance){ set key left top; unset logscale y; set format y '%g'; set xrange[0.2:0.6];}
fanc(tau, tH) = conductance? (log10(1.0 / $3)) : (rescale_thouless? $3 * $4 : $3)

set xlabel 'g'
set ylabel (conductance? "log_{10}(t_{H}/t_{Th})" : (rescale_thouless? "t_{Th}" : "{/Symbol t}_{Th}")) rotate by 0
if(scaling == 1){ plot for[L=9:14] _name_th(J, L, w2) u ($1 == h? $2 : NaN):(fanc($3,$4)) w lp pt 4 ps 0.75 title sprintf("L=%d", L)#, 0.1*x**-4 w l lc rgb 'black'
} else {
set key left top
	if(scaling == 2){
		if(h_vs_g){ 
			plot for[i=h0:hend:dh] _name_th_L(J, g, 0.01*i, w) u 1:($2*$3) w lp pt 4 ps 0.75 title sprintf("h=%.2f", 0.01*i), _name_th_L(J, g, 0.01*hend, w) u 1:($3) w l ls 1 dt (3,5,10,5) lc rgb 'black' lw 2 title "t_H"
		} else {
			plot for[i=g0:gend:dg] _name_th_L(J, 0.01*i, h, w) u 1:($2*$3) w lp pt 4 ps 0.75 title sprintf("g=%.2f", 0.01*i), _name_th_L(J, 0.01*gend, h, w) u 1:($3) w l ls 1 dt (3,5,10,5) lc rgb 'black' lw 2 title "t_H"
		}
	} else {
		Lx = 11
		set multiplot
		ss = "\n";
		while(Lx <= 14){
			@MARGIN; @SCALE
			sizeee = (Jend-J0)/dJ + 1
			array tau[sizeee];	array Jarr[sizeee];
    		do for[i=J0:Jend:dJ]{
				idx = (i - J0)/dJ + 1
        	    stats _name_th(0.01*i, Lx, w2) using ($1 == h && $2 == g? fanc($3,$4) : NaN) every ::0::0 nooutput;   tau[idx] = STATS_min; 
				Jarr[idx] = 0.01*i;
				print Jarr[idx], tau[idx]
			}
			#plot for[i=1:sizeee] '+' u (Jarr[i]):(tau[i]) w lp pt 6 ps 0.75 notitle
			plot Jarr using (Jarr[$1]):(tau[$1]) w lp pt 6 ps 1.75 lc (Lx-9) t ss.sprintf("L=%d", Lx)
			Lx = Lx+1
			unset xlabel; unset ylabel; unset title;
			ss = ss."\n\n\n" 
		}
		unset multiplot
		#set key left bottom
		#plot for[i=1:w_num] _name_th(J, L, w_list[i]) u ($1 == h? $2 : NaN):(fanc($3,$4)) w lp pt 4 ps 0.75 title sprintf("w=%.2f", w_list[i]), _name_th(J, L, w_list[w_num]) u ($1 == h? $2 : NaN):4 w l ls 1 dt (3,5,10,5) lc rgb 'black' lw 2 title "t_H"
	}
}


exit;

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