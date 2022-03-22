dir_base='../../results/disorder/PBC/'
dir = dir_base.'entropy/'
out_dir = 'Time_Evolution/'

reset 

#------------------------------------ PREAMBLE
set autoscale
set term qt size 900, 900 font "Helvetica,14"
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
LIN_LOG="unset logscale xy; set logscale x; set format x '10^{%L}'";
UNSET = "unset tics; unset xlabel; unset ylabel; unset title; unset border;"

set style line 1 dt (3,5,10,5) lc rgb "black" lw 2.5
set style line 2 dt (3,3) lc rgb "red" lw 2.5
set style line 3 dt (8,8) lc rgb "blue" lw 2.5
set style line 4 dt (1,1) lc rgb "green" lw 2.5
set fit quiet
#------------------------------------ Margins for each row resp. column
TMARGIN = "set tmargin at screen 0.91; set bmargin at screen 0.56"
BMARGIN = "set tmargin at screen 0.46; set bmargin at screen 0.11"
# TICS
NOYTICS = "set format y '';"
YTICS = "set format y '%g';"

#------------------------------------ PARAMETERS
L = 13; 
g = 0.7
g2 = 0.9 
h = 0.8;
J0 = 0.; g_knot = 0.; 
w = 0.01;

subsystem_size=3		# subsystem size
scaling = 1				# size scaling=1 or h-scaling=0 or 	g-scaling=2 or subsystem_size=3
what_to_plot = 1   		# subsystem size=0, time evolution=1  or  eigenstates parabola=2
rescale_x_axis = 0		# rescale x ax0s?
compare_scales = 0		# plot 4 panels, 2g's anf log-log and lin-log for each
compare_to_lanczos = 1	# compare all results to lanczos eovlution (if available)
if(scaling == 2) compare_scales = 0;
	h0 = 20;	hend = 300;		dh = 20;
	g0 = 20;	gend = 70;		dg = 10;
	L0 = 11;	Lend = 14; 		dL = 1;
	
str_name = (what_to_plot==1? "TimeEvolution" : (what_to_plot==0? "SubsystemSize" : "Eigenstates"));

fit = 0
x_min = 3e0; x_max = 1e1;

my_title = "Entalgement entropy ".str_name." for \n";
if(scaling != 1) my_title = my_title.sprintf("L=%d, ",L)
if(scaling != 2) my_title = my_title.sprintf("g=%.2f, ", g)
if(scaling != 0) my_title = my_title.sprintf("h=%.2f, ", h) 
if(scaling != 3) my_title = my_title.sprintf("L_A = %d", L / 2.-3+subsystem_size) 
#set title my_title

fileexist(name)=1#int(system("if exist \"".name."\" ( echo 1) else (echo 0)"))

f(t) = a*log((t-t0)) + b
f_plot(a,b, t0,t) = (t < x_min || t > x_max)? NaN : a*log((t-t0)) + b

	#------------------------------------ DATA, FIT AND PLOT
	name(x) = 0; key_title(x) = 0; rescale_X(x,y) = x;
	name2(x) = 0;
	name_lancz(x) = 0;
	i0 = 0; iend = 0; di = 1;
		if(scaling == 0){
			name(x) = dir.str_name.sprintf("_L=%d,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", L, J0, g, g_knot, 0.01*x, w);	key_title(x) = sprintf("h=%.2f", x/100.)
			name_lancz(x) = dir.'Lanczos/'.str_name.sprintf("_L=%d,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", L, J0, g, g_knot, 0.01*x, w);
			i0 = h0; iend = hend; di = dh;
		} else{
			if(scaling == 1){
				name(x) = dir.str_name.sprintf("_L=%d,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", x, J0, g, g_knot, h, w);	key_title(x) = sprintf("L=%d",x);
				name2(x) = dir.str_name.sprintf("_L=%d,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", x, J0, g2, g_knot, h, w);
				name_lancz(x) = dir.'Lanczos/'.str_name.sprintf("_L=%d,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", x, J0, g, g_knot, h, w);
				i0 = L0; iend = Lend; di = dL;
			} else{
				if(scaling == 2){
					rescale_X(x,y) = rescale_x_axis? x * (1e-2*y)**2 : x;
					name(x) = dir.str_name.sprintf("_L=%d,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", L, J0, 0.01*x, g_knot, h, w); 
					name_lancz(x) = dir.'Lanczos/'.str_name.sprintf("_L=%d,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", L, J0, 0.01*x, g_knot, h, w); 
					key_title(x) = sprintf("g=%.2f",0.01*x)
					i0 = g0; iend = gend; di = dg;
				} else{
					name(x) = dir.str_name.sprintf("_L=%d,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", L, J0, g, g_knot, h, w); 
					name_lancz(x) = dir.'Lanczos/'.str_name.sprintf("_L=%d,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", L, J0, g, g_knot, h, w); 
					key_title(x) = sprintf("L_A=%d", L / 2. + x);
					i0 = -1; iend = 1; di = 1;
				}
			}
		}
		array xVal[(iend-i0)/di+1];
		array S1[(iend-i0)/di+1]; array S1_page[(iend-i0)/di+1];
		array S2[(iend-i0)/di+1]; array S2_page[(iend-i0)/di+1];
		array S3[(iend-i0)/di+1]; array S3_page[(iend-i0)/di+1];
		array a_list[(iend-i0)/di+1]; array b_list[(iend-i0)/di+1]; array t0_list[(iend-i0)/di+1];
		print "g=", g, "\npar\tpage analytical values: L/2-1, L/2, L/2+1"
		
		rescale(x,i) = x;#scaling==1? x * 2.0 / (i*log(2) - 0.5) : x * 2.0 / (L*log(2) - 0.5)
		
		do for[i=i0:iend:di]{
			name = name(i)
			if(fileexist(name)){
				if(fit){
					fit[x_min:x_max][*:*] f(x) name u 1:(rescale($2,i)) via a, b, t0
					a_list[(i-i0)/di+1] = a; b_list[(i-i0)/di+1] = b; t0_list[(i-i0)/di+1] = t0;
				}
				stats name every ::0::0 using 2 nooutput; S1[(i-i0)/di+1] = rescale(STATS_min,i);
				stats name every ::1::1 using 2 nooutput; S2[(i-i0)/di+1] = rescale(STATS_min,i);
				stats name every ::2::2 using 2 nooutput; S3[(i-i0)/di+1] = rescale(STATS_min,i);
				if(what_to_plot==0){
					stats name every ::0::0 using 3 nooutput;	S1_page[(i-i0)/di+1] = rescale(STATS_min,i);
					stats name every ::1::1 using 3 nooutput; 	S2_page[(i-i0)/di+1] = rescale(STATS_min,i);
					stats name every ::2::2 using 3 nooutput; 	S3_page[(i-i0)/di+1] = rescale(STATS_min,i);
					xVal[(i-i0)/di+1] = scaling == 1? 1.0 / (i**3) : 0.01*i;
				}
				print key_title(i),"  ", S1_page[(i-i0)/di+1], S2_page[(i-i0)/di+1], S3_page[(i-i0)/di+1], "\t\t", xVal[(i-i0)/di+1], S1[(i-i0)/di+1], S2[(i-i0)/di+1], S3[(i-i0)/di+1]
			}
			else{
				print key_title(i)," -----------------------------------------------------------------------------"
			}		
		}
	
	page(x) = (L * log(2.0) - 0.5) / 2.
	#set logscale y
	#set xrange[0:1]
	#plot dir.sprintf("compare_to_disorder_L=%d,g=%.2f,h=0.80,k=0,p=1,x=1.dat", L, g) u ($1 / (L+0.0)):($2 / page($1)) w p ps 2 t 'w=1e-2',\
	#	dir.sprintf("compare_to_disorder_L=%d,g=%.2f,h=0.80,k=0,p=1,x=1.dat", L, g) u ($1 / (L+0.0)):($3 / page($1)) w p ps 2 t 'w=1e-3',\
	#	dir.sprintf("compare_to_disorder_L=%d,g=%.2f,h=0.80,k=0,p=1,x=1.dat", L, g) u ($1 / (L+0.0)):($4 / page($1)) w p ps 2 t 'w=1e-4',\
	#	dir.sprintf("compare_to_disorder_L=%d,g=%.2f,h=0.80,k=0,p=1,x=1.dat", L, g) u ($1 / (L+0.0)):($5 / page($1)) w p ps 2 t 'k=0,p=1',\
	#	dir.sprintf("compare_to_disorder_L=%d,g=%.2f,h=0.80,k=0,p=1,x=1.dat", L, g) u ($1 / (L+0.0)):($6 / page($1)) w p ps 2 t 'k=1', 1-abs(1-2*x) w l ls 1 notitle
	#exit;
	RANGE="set xrange[1e-2:1000]; set yrange[1e-3:6];"
	MARGIN = compare_scales? "set lmargin at screen 0.10; set rmargin at screen 0.54; set bmargin at screen 0.10; set tmargin at screen 0.54;"\
					: "set lmargin at screen 0.10; set rmargin at screen 0.98; set bmargin at screen 0.10; set tmargin at screen 0.98;"
	if(what_to_plot==1){
		if(compare_scales){
			MARGIN2 = "set lmargin at screen 0.54; set rmargin at screen 0.98; set bmargin at screen 0.10; set tmargin at screen 0.54;";
			MARGIN3 = "set lmargin at screen 0.10; set rmargin at screen 0.54; set bmargin at screen 0.54; set tmargin at screen 0.98;";
			MARGIN4 = "set lmargin at screen 0.54; set rmargin at screen 0.98; set bmargin at screen 0.54; set tmargin at screen 0.98;";
			LABEL1="unset label 2; set label 1 at 0.3,2.0 sprintf('g=%.2f',g) front"
			LABEL2="unset label 1; set label 2 at 0.3,2.0 sprintf('g=%.2f',g2) front"
			
			set key left top
			set multiplot
			@LOG_LOG; @MARGIN; @RANGE;  @LABEL1;
				set ylabel 'S(t)'; set xlabel 't'; plot for[i=i0:iend:di] name(i) u 1:2 w lp ps 0.75 pt 5 notitle, 0.15*x w l ls 1 notitle, 1.75*log(x)-2.5 w l ls 2 notitle
			@LOG_LOG; @MARGIN2; @RANGE; @LABEL2;
				unset ylabel; set format y ''; plot for[i=i0:iend:di] name2(i) u 1:2 w lp ps 0.75 pt 5 notitle, 1.05*x-0.1 w l ls 1 notitle, 2.5*log(x)+0.4 w l ls 2 notitle
			@LIN_LOG; @MARGIN3; @RANGE; @LABEL1;
				unset xlabel; set format x ''; set format y '%g';  plot for[i=i0:iend:di] name(i) u 1:2 w lp ps 0.75 pt 5 title key_title(i), 0.15*x w l ls 1, 1.75*log(x)-2.5 w l ls 2
				#if(scaling==1){ @UNSET; @MARGIN3; @RANGE; plot for[L=L0:Lend:dL] ( x > 30? (L * log(2.0) - 0.5)/2. : NaN) w l ls 4 notitle; }; @FORMAT; set border;
			@LIN_LOG; @MARGIN4; @RANGE; @LABEL2; set key right bottom
				unset xlabel; unset ylabel; set format x ''; set format y '';  plot for[i=i0:iend:di] name2(i) u 1:2 w lp ps 0.75 pt 5 notitle, 1.05*x-0.1 w l ls 1, 2.5*log(x)+0.4 w l ls 2
				#if(scaling==1){ @UNSET; @MARGIN4; @RANGE; plot for[L=L0:Lend:dL] ( x > 30? (L * log(2.0) - 0.5)/2. : NaN) w l ls 4 notitle; }; @FORMAT; set border;
			unset multiplot
		} else{
			fbase(x,a,b) = a*log(x) - b
			f(x,a,b) = (fbase(x,a,b) < 0.1 || fbase(x,a,b) > 0.9 - 0.5*b)?NaN : fbase(x,a,b);
			set format x "10^{%L}"; #set format y "10^{%L}"
			set xlabel 't'
			set ylabel 'S(t)'
			set key left top font ",16"
			#label_initial = "|s_j{/Symbol \361}=cos({/Symbol \161}_j/2)|{/Symbol \255}_j{/Symbol \361}+e^{i{/Symbol f}_j}sin({/Symbol \161}_j/2)|{/Symbol \257}_j{/Symbol \361}\n\n|{/Symbol \171}(t=0){/Symbol \361}={/*2{/Symbol \304}}_j |s_j{/Symbol \361}"
			#set label 1 at graph 0.02,0.9 sprintf("%s",label_initial) front
			
			set logscale x;
			set multiplot
			@MARGIN; @RANGE;
			plot for[i=i0:iend:di] name(i) u (rescale_X($1,i)):(rescale($2,i)) w l title key_title(i)#,\
					f(x,0.4,0.0) w l ls 1 title 'ln(x)', f(x,0.28,0.0) w l ls 1 notitle, f(x,0.28,0.4) w l ls 1 notitle,\
					f(x,0.28,0.65) w l ls 1 notitle, f(x,0.28,0.2) w l ls 1 notitle, f(x,0.35,-0.2) w l ls 1 notitle
			if(fit){
				@MARGIN; @UNSET; @RANGE;
				plot for[i=1:(iend-i0)/di + 1] f_plot(a_list[i], b_list[i], t0_list[i], x) w l ls 1 lw 2.5 notitle
			}
			if(compare_to_lanczos){
				@MARGIN; @UNSET; @RANGE;
				plot for[i=i0:iend:di] name_lancz(i) u (rescale_X($1,i)):(rescale($2,i)) w p ps 0.75 pt 5 notitle
			}
			#plot for[i=i0:iend:di] name(i) u 1:($1 < 10 ? NaN : L/2.*log(2) - 0.5) w l ls 2 notitle 
			unset multiplot
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

