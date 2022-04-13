dir_base='../results/disorder/PBC/'
dir = dir_base.'Entropy/'
out_dir = 'Time_Evolution/'

reset 

#------------------------------------ PREAMBLE
set autoscale
set mxtics
set mytics

FORMAT = "set style line 12 lc rgb '#ddccdd' lt 1 lw 1.5; \
set style line 13 lc rgb '#ddccdd' lt 1 lw 0.5; \
set grid xtics mxtics ytics mytics back ls 12, ls 13;\
set size square;\
set xtics mirror;\
set ytics mirror;"
@FORMAT
LOG_LOG="set logscale xy; set format x '10^{%L}'; set format y '10^{%L}'";
LIN_LOG="unset logscale xy; set logscale x; set format x '10^{%L}'; set format y '%g';";
LOG_LIN="unset logscale xy; set logscale y; set format x '10^{%L}'; set format y '%g';";
LIN_LIN="unset logscale xy; set format x '%g'; set format y '%g';";
UNSET = "unset tics; unset xlabel; unset ylabel; unset title; unset border;"

set style line 1 dt (3,5,10,5) lc rgb "black" lw 2.5
set style line 2 dt (3,3) lc rgb "red" lw 2.5
set style line 3 dt (8,8) lc rgb "blue" lw 2.5
set style line 4 dt (1,1) lc rgb "green" lw 2.5
set fit quiet

#------------------------------------ PARAMETERS
L = 17
g = 0.40
g2 = 0.9 
h = 0.8;
J0 = 0.; g_knot = 0.; 
w = 0.01;

subsystem_size=3		# subsystem size
scaling = 2		# size scaling=1 or h-scaling=0 or 	g-scaling=2 or subsystem_size=3
what_to_plot = 1   		# subsystem size=0, time evolution=1  or  eigenstates parabola=2

nu=-0.4
fX(x) = exp(nu / x**2.0)
#fX(x) = x**2.3
rescale_x_axis = 0		# rescale x ax0s?
rescale_by_page=0		# rescale byu Page value

compare_scales = 0		# plot 4 panels, 2g's anf log-log and lin-log for each
compare_to_lanczos = 0	# compare all results to lanczos eovlution (if available)
plot_exponent = 1

if(plot_exponent == 0 && scaling == 2) compare_scales = 0;
	h0 = 20;	hend = 300;		dh = 20;
	g0 = 30;	gend = 90;		dg = 5;
	L0 = 12;	Lend = 17; 		dL = 1;

plot_only_lanczos = 0
#if(plot_exponent) plot_only_lanczos = 0;
str_name = (what_to_plot==1? "TimeEvolution" : (what_to_plot==0? "SubsystemSize" : "Eigenstates"));

fit = 0
x_min = 3e0; x_max = 1e1;

fileexist(name)=system("[ -f '".name."' ] && echo '1' || echo '0'") + 0#int(system("if exist \"".name."\" ( echo 1) else (echo 0)"))

f(t) = a*log((t-t0)) + b
f_plot(a,b, t0,t) = (t < x_min || t > x_max)? NaN : a*log((t-t0)) + b

	#------------------------------------ DATA, FIT AND PLOT
	
	page(x) = scaling == 1? (x * log(2.0) - 0.5) / 2. : (L * log(2.0) - 0.5) / 2.
	rescale(x,i) = rescale_by_page? x / page(i) : x;
	#dir = '../results/disorder/PBC/TimeEvolution/j=0/exponent/'
	#str_name="SigmaZ_j=0"
	name(x) = 0; key_title(x) = 0; rescale_X(x,y) = x;
	name2(x) = 0; 		name_lancz(x) = 0;		name_exp(x) = 0;
	i0 = 0; iend = 0; di = 1;
	if(plot_only_lanczos || (L >= 16 && scaling != 1 && !compare_to_lanczos)){ dir=dir.'Lanczos/'; };
	
		if(scaling == 0){
			base(x) = sprintf("_L=%d,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", L, J0, g, g_knot, 0.01*x, w); 		key_title(x) = sprintf("h=%.2f", x/100.)
			name(x) = dir.str_name.base(x);		name_lancz(x) = dir.'Lanczos/'.str_name.base(x); 	name_exp(x) = dir.'exponent/'.base(x);
			i0 = h0; iend = hend; di = dh;
		} else{
			if(scaling == 1){
				base(x, gx) = sprintf("_L=%d,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", x, J0, gx, g_knot, h, w);		key_title(x) = sprintf("L=%d",x);
				suf(x) = (x >= 16 && !plot_only_lanczos? 'Lanczos/' : '')
				name(x) = dir.suf(x).str_name.base(x, g);		name2(x) = dir.suf(x).str_name.base(x, g2)
				name_lancz(x) = dir.'Lanczos/'.str_name.base(x, g);		name_exp(x) = dir.suf(x).'exponent/'.base(x, g);
				i0 = L0; iend = Lend; di = dL;
			} else{
				if(scaling == 2){
					rescale_X(x,y) = rescale_x_axis? x * fX(0.01*y) : x;
					base(x) = sprintf("_L=%d,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", L, J0, 0.01*x, g_knot, h, w);		key_title(x) = sprintf("g=%.2f",0.01*x)
					name(x) = dir.str_name.base(x);		name_lancz(x) = dir.'Lanczos/'.str_name.base(x);	name_exp(x) = dir.'exponent/'.base(x);
					i0 = g0; iend = gend; di = dg;
				} else{
					if(scaling == 3){
						base = sprintf("_L=%d,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", L, J0, g, g_knot, h, w);  		key_title(x) = sprintf("L_A=%d", L / 2. + x);
						name(x) = dir.str_name.base;		name_lancz(x) = dir.'Lanczos/'.str_name.base;		name_exp(x) = dir.'exponent/'.base;
						i0 = -1; iend = 1; di = 1;
					} else {
						base = sprintf("_L=%d,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", L, J0, g, g_knot, h, w);  		key_title(x) = sprintf("L=%d, g=%.2f,h=%.2f", L, g, h);
						name(x) = dir.str_name.base;		name_lancz(x) = dir.'Lanczos/'.str_name.base;		name_exp(x) = dir.'exponent/'.base;
						i0 = 0; iend = 0; di = 1;
					}
				}
			}
		}
		array xVal[(iend-i0)/di+1];
		array S1[(iend-i0)/di+1]; array S1_page[(iend-i0)/di+1];
		array S2[(iend-i0)/di+1]; array S2_page[(iend-i0)/di+1];
		array S3[(iend-i0)/di+1]; array S3_page[(iend-i0)/di+1];
		array a_list[(iend-i0)/di+1]; array b_list[(iend-i0)/di+1]; array t0_list[(iend-i0)/di+1];
		
		if(what_to_plot == 0){ print "g=", g, "\npar\tpage analytical values: L/2-1, L/2, L/2+1" };

		do for[i=i0:iend:di]{
			name = name(i)
			if(fileexist(name)){
				if(fit){
					fit[x_min:x_max][*:*] f(x) name u 1:(rescale($2,i)) via a, b, t0
					a_list[(i-i0)/di+1] = a; b_list[(i-i0)/di+1] = b; t0_list[(i-i0)/di+1] = t0;
				}
				if(what_to_plot == 0){
						stats name every ::0::0 using 2 nooutput; S1[(i-i0)/di+1] = rescale(STATS_min,i);
						stats name every ::1::1 using 2 nooutput; S2[(i-i0)/di+1] = rescale(STATS_min,i);
						stats name every ::2::2 using 2 nooutput; S3[(i-i0)/di+1] = rescale(STATS_min,i);
						stats name every ::0::0 using 3 nooutput;	S1_page[(i-i0)/di+1] = rescale(STATS_min,i);
						stats name every ::1::1 using 3 nooutput; 	S2_page[(i-i0)/di+1] = rescale(STATS_min,i);
						stats name every ::2::2 using 3 nooutput; 	S3_page[(i-i0)/di+1] = rescale(STATS_min,i);
						xVal[(i-i0)/di+1] = scaling == 1? 1.0 / (i**3) : 0.01*i;
					print key_title(i),"  ", S1_page[(i-i0)/di+1], S2_page[(i-i0)/di+1], S3_page[(i-i0)/di+1], "\t\t", xVal[(i-i0)/di+1], S1[(i-i0)/di+1], S2[(i-i0)/di+1], S3[(i-i0)/di+1]
				}
			}
			else{
				print key_title(i)," -----------------------------------------------------------------------------"
			}		
		}
	
	#set logscale y
	#set xrange[0:1]
	#plot dir.sprintf("compare_to_disorder_L=%d,g=%.2f,h=0.80,k=0,p=1,x=1.dat", L, g) u ($1 / (L+0.0)):($2 / page($1)) w p ps 2 t 'w=1e-2',\
	#	dir.sprintf("compare_to_disorder_L=%d,g=%.2f,h=0.80,k=0,p=1,x=1.dat", L, g) u ($1 / (L+0.0)):($3 / page($1)) w p ps 2 t 'w=1e-3',\
	#	dir.sprintf("compare_to_disorder_L=%d,g=%.2f,h=0.80,k=0,p=1,x=1.dat", L, g) u ($1 / (L+0.0)):($4 / page($1)) w p ps 2 t 'w=1e-4',\
	#	dir.sprintf("compare_to_disorder_L=%d,g=%.2f,h=0.80,k=0,p=1,x=1.dat", L, g) u ($1 / (L+0.0)):($5 / page($1)) w p ps 2 t 'k=0,p=1',\
	#	dir.sprintf("compare_to_disorder_L=%d,g=%.2f,h=0.80,k=0,p=1,x=1.dat", L, g) u ($1 / (L+0.0)):($6 / page($1)) w p ps 2 t 'k=1', 1-abs(1-2*x) w l ls 1 notitle
	#exit;
	RANGE = rescale_x_axis? "set xrange[1e-1:2e2]; set yrange[3e-2:6];" : "set xrange[1e-1:200]; set yrange[3e-1:6];"
	RANGE2="set xrange[3e-1:10]; set yrange[1e-1:6];"
	MARGIN = compare_scales? "set lmargin at screen 0.10; set rmargin at screen 0.54; set bmargin at screen 0.10; set tmargin at screen 0.54;"\
					: "set lmargin at screen 0.10; set rmargin at screen 0.98; set bmargin at screen 0.10; set tmargin at screen 0.98;"
	
	if(plot_exponent && compare_scales) { set term qt size 1500, 750 font "Helvetica,18"; }
	else { set term qt size 900, 900 font "Helvetica,18"; }
	if(what_to_plot==1){
		if(compare_scales){
			MARGIN2 = plot_exponent == 0? "set lmargin at screen 0.54; set rmargin at screen 0.98; set bmargin at screen 0.10; set tmargin at screen 0.54;"\
							: "set lmargin at screen 0.06; set rmargin at screen 0.48; set bmargin at screen 0.10; set tmargin at screen 0.98;";
			MARGIN3 = plot_exponent == 0? "set lmargin at screen 0.10; set rmargin at screen 0.54; set bmargin at screen 0.54; set tmargin at screen 0.98;"\
							: "set lmargin at screen 0.56; set rmargin at screen 0.98; set bmargin at screen 0.10; set tmargin at screen 0.98;";
			MARGIN4 = "set lmargin at screen 0.54; set rmargin at screen 0.98; set bmargin at screen 0.54; set tmargin at screen 0.98;";
			LABEL1="unset label 2; set label 1 at 0.3,2.0 sprintf('g=%.2f',g) front"
			LABEL2="unset label 1; set label 2 at 0.4,2.0 sprintf('g=%.2f',g2) front"
			
			if(plot_exponent == 0){
				set key left top
				f_lin_left(x) = 0.3*(x**0.8+0); 	name_lin_left = "\n0.3*x^{0.8}"
				f_log_left(x) = 1.8*log(x)-2.2; 	name_log_left = "\n1.75*log(x)-2.5"
				f_lin_right(x) = 1.05*x-0.1; 		name_lin_right = '1.05*x-0.1'
				f_log_right(x) = 2.8*log(x)+0.0; 	name_log_right = '2.8*log(x)'
				set multiplot
				@LOG_LOG; @MARGIN; @RANGE;  @LABEL1; # LEFT - BOTTOM
					set ylabel 'S(t)'; set xlabel 't'; plot for[i=i0:iend:di] name(i) u 1:2 w lp ps 0.75 pt 5 lw 2 notitle, f_lin_left(x) w l ls 1 notitle, f_log_left(x) w l ls 2 notitle
				@LOG_LOG; @MARGIN2; @RANGE2; @LABEL2; # RIGHT - BOTTOM
					unset ylabel; set format y ''; plot for[i=i0:iend:di] name2(i) u 1:2 w lp ps 0.75 pt 5 lw 2 notitle, f_lin_right(x) w l ls 1 notitle, f_log_right(x) w l ls 2 notitle
				@LIN_LOG; @MARGIN3; @RANGE; @LABEL1;# LEFT - TOP
					unset xlabel; set format x ''; set format y '%g';  plot for[i=i0:iend:di] name(i) u 1:2 w lp ps 0.75 pt 5 lw 2 title key_title(i), f_lin_left(x) w l ls 1 t name_lin_left, f_log_left(x) w l ls 2 t name_log_left
					#if(scaling==1){ @UNSET; @MARGIN3; @RANGE; plot for[L=L0:Lend:dL] ( x > 30? (L * log(2.0) - 0.5)/2. : NaN) w l ls 4 notitle; }; @FORMAT; set border;
				@LIN_LOG; @MARGIN4; @RANGE2; @LABEL2; # RIGHT - TOP
					set key right bottom
					unset xlabel; unset ylabel; set format x ''; set format y '';  plot for[i=i0:iend:di] name2(i) u 1:2 w lp ps 0.75 pt 5 lw 2 notitle, f_lin_right(x) w l ls 1 t name_lin_right, f_log_right(x) w l ls 2 t name_log_right
					#if(scaling==1){ @UNSET; @MARGIN4; @RANGE; plot for[L=L0:Lend:dL] ( x > 30? (L * log(2.0) - 0.5)/2. : NaN) w l ls 4 notitle; }; @FORMAT; set border;
				unset multiplot
			} else {				
				set multiplot
				@LOG_LOG; @MARGIN2; @RANGE; 
					set key left top; set ylabel 'S(t)'; set xlabel 't'; plot for[i=i0:iend:di] name(i) u 1:2 w lp ps 0.75 pt 5 lw 2 t key_title(i)#, sqrt(x) w l ls 0 t "t^{1/2}", x w l ls 1  t 't', log(x) w l ls 2 t 'ln(t)'
				@LIN_LIN; @MARGIN3; set xrange[6e-1:1.5e1]; set yrange[1e-3:2e0];
					set ylabel '{/Symbol g}(t)'; set xlabel 't'; plot for[i=i0:iend:di] name_exp(i) u 1:($2) w lp ps 0.75 pt 5 lw 2 t key_title(i)
				unset multiplot
			}
		} else{
			fbase(x,a,b) = a*log(x) - b
			f(x,a,b) = (fbase(x,a,b) < 0.1 || fbase(x,a,b) > 40)? NaN : fbase(x,a,b);
			set format x "10^{%L}"; 
			unset logscale y;#set format y "10^{%L}"
			set xlabel 't'
			set ylabel 'S(t)'
			set key left top font ",16"
			@LIN_LOG
			if(plot_exponent == 0){
				set multiplot
				@MARGIN; @RANGE;
				STYLE=compare_to_lanczos? "w l" : "w lp ps 0.75 pt 5";
				#plot for[i=i0:iend-di:di] '< paste '.name(i).' '.name(i+di) u 1:($2/$4) w lp ps 0.75 pt 5 lw 2 t key_title(i).' / '.key_title(i+di)
				plot for[i=i0:iend:di] name(i) u (rescale_X($1,i)):(rescale(abs($2),i)) @STYLE lw 2 title key_title(i)#,\
						#f(x,1.8,4.45) w l ls 1 title 'ln(x)', 0.3*x**0.8+0 w l ls 2 t "\n0.3t^{0.8}"#, f(x,0.28,0.4) w l ls 1 notitle#,\
						#f(x,0.68,0.65) w l ls 1 notitle, f(x,0.28,0.2) w l ls 1 notitle, f(x,0.1,-10) w l ls 1 notitle
				if(fit){
					@MARGIN; @UNSET; @RANGE;
					plot for[i=1:(iend-i0)/di + 1] f_plot(a_list[i], b_list[i], t0_list[i], x) w l ls 1 lw 2.5 notitle
				}
				if(compare_to_lanczos){
					@MARGIN; @UNSET; @RANGE;
					plot for[i=i0:iend:di] name_lancz(i) u (rescale_X($1,i)):(rescale($2,i)) w p ps 0.75 pt 5 notitle
					#plot for[i=i0:iend:di] name(i) u (rescale_X($1,i)):(rescale($3,i)) w p ps 0.75 pt 5 notitle
				}
				#plot for[i=i0:iend:di] name(i) u 1:($1 < 10 ? NaN : L/2.*log(2) - 0.5) w l ls 2 notitle 
				unset multiplot
			} else{
				@MARGIN; set xrange[6e-1:100]; set yrange[0:2e0];
				plot for[i=i0:iend:di] name_exp(i) u (rescale_X($1,i)):($2) w lp ps 0.75 pt 5 lw 2 title key_title(i)#, 1/x w l ls 0 lw 2 notitle
			}
		}
	} else {
		if(what_to_plot==0){
			set ylabel '<S_A> - S_{page}'
			set xlabel '1 / L^2'
			set xrange[0.0:*]; set yrange[0:*];
			#set logscale x
			#set format x "10^{%L}"
			@MARGIN;
			plot S2 using (xVal[$1]):(abs(S2[$1] - S2_page[$1])) w lp ps 1.2 pt 6 title 'L_A= L /2',\
				80.0*x w l ls 0 lw 2 title 'linear fit',\
				S3 using (xVal[$1]):(abs(S3[$1] - S3_page[$1])) w lp ps 1.2 pt 7 title 'L_A= L /2 + 1',\
				S1 using (xVal[$1]):(abs(S1[$1] - S1_page[$1])) w lp ps 1.2 pt 6 title 'L_A= L /2 - 1'
			@UNSET;
		}
		else{
			set ylabel '<n|S_A|n>'
			set xlabel 'E_n/L'
			set key right top
			if(scaling == 3){
				plot for[i=i0:iend:di] name(i) u 1:(column(i+3)) w p ps 0.75 pt 5 title key_title(i)
			} else {
				plot for[i=i0:iend:di] name(i) u 1:(column(subsystem_size)) w p ps 0.75 pt 5 title key_title(i)
			}
		}
	}

