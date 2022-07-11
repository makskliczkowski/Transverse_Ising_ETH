dir_base='../results/ISING/disorder/PBC/'
dir = dir_base.'IntegratedResponseFunction/'
out_dir = 'Integrated_Spectral_Function/'
reset 

#------------------------------------ PREAMBLE
set autoscale
use_png = 0		# 1 if use png output, and 0 for qt output
if(use_png) { set term pngcairo size 1200, 1200 font sprintf("Helvetica,%d",20); }
else {set term qt size 900, 900 font sprintf("Helvetica,%d",20); }
load './gnuplot-colorbrewer-master/diverging/RdYlGn.plt'

set mxtics
set mytics

FORMAT = "set style line 12 lc rgb '#ddccdd' lt 1 lw 1.5; \
set style line 13 lc rgb '#ddccdd' lt 1 lw 0.5; \
set grid xtics mxtics ytics mytics back ls 12, ls 13;\
set size square;\
set xtics mirror;\
set ytics mirror;"
@FORMAT
x_log = 1; 
y_log = 1;
SCALE = x_log? "unset logscale xy; set logscale x; set format x '10^{%L}'; " : "set format x '%g'"
SCALE = SCALE.(y_log? "set logscale y; set format y '10^{%L}'" : "set format y '%g'")
@SCALE

fileexist(name)=1#"[ -f name ] && echo 1 || echo 0"
#int(system("if exist \"".name."\" ( echo 1) else (echo 0)"))
	
UNSET = "unset tics; unset xlabel; unset ylabel; unset title; unset border;"
set fit quiet
#------------------------------------ Margins for each row resp. column
TMARGIN = "set tmargin at screen 0.91; set bmargin at screen 0.56"
BMARGIN = "set tmargin at screen 0.46; set bmargin at screen 0.11"
# TICS
NOYTICS = "set format y '';"
YTICS = "set format y '%g';"

#------------------------------------ PARAMETERS
L = 10;
J=1.0 
g = 0.9
h = 0.8;
J0 = 0.; g_knot = 0.; 
w = 1.0;

x_range_min=1e-4

integrated_by_hand = 0 #integrated time evolution?
if(integrated_by_hand) cd '.\integrated'
rescale = 0				# rescale the spectral function by f(w, L)?
site = 0				# site at which the operator acts
scaling = 2				# size scaling=1 or h-scaling=0 or 	g-scaling=2	or 	q/j-scaling=3 or realisation=4 or user=5
q_vs_j = 0				# =1 - evolution of Sz_q, else ecol of Sz_j
operator = 1	 		# 1-SigmaZ , 0-Hq :local

two_panels = 1			# plot integrated spectral function next to respons function
smoothed = 1			# smoothed derivative?
#-- IntegratedSpecFun
plot_normalized = 0		# plot renormalized to 1st peak
plot_exponent = 0;	# plot as integrated response function
#-- SpectralFun
plot_derivative = 1		# use derivative as integrated spectal function
substract_LTA = 0

rescale = 0
nu = 2		# power on L
if(scaling != 1) rescale = 0;
LIOM = 0				# plot LIOMs?
local = 0

	h0 = 10;	hend = 150;		dh = 10;
	g0 = 10;	gend = 100;		dg = 10;
	L0 = 10;	Lend = 15; 		dL = 1;

use_fit = 0
which_fit = 1						#=0 - a*exp(-t/b); =1 - a-b*log(t)
x_min = 1e-5; x_max = 3e-2;

#------------------------------------------ SET OPERATOR
op = ""
if(operator == 0) {op = "H"; }
if(operator == 1) {op = "SigmaZ";}
if(operator == 2) {op = "TFIM_LIOM_plus";}
if(operator == 3) {op = "TFIM_LIOM_minus";}

_str(x) = (q_vs_j? "q" : "j").sprintf("=%d",x);
if(operator > 1){ _str(x) = "n".sprintf("=%d",x); }
str(x) = _str(site)
if(scaling == 1) { str(x) = (site == -1)? _str(x/2) : _str(site) }
if(scaling == 3) { str(x) = _str(x) }

	#------------------------------------ DATA
	output_name = ""
	_name(x) = 0; _key_title(x) = 0;
	i0 = 0; iend = 0; di = 1;
		if(scaling == 0){
			_name(x) = op.sprintf("_".str(site)."_L=%d,J=%.2f,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", L, J, J0, g, g_knot, 0.01*x, w);	_key_title(x) = sprintf("h=%.2f", x/100.)
			i0 = h0; iend = hend; di = dh; 	out_dir = out_dir."h_scaling/";
			output_name = output_name.op.sprintf("_".str(site)."_L=%d,J=%.2f,J0=%.2f,g=%.2f,g0=%.2f,w=%.2f", L, J, J0, g, g_knot, w);
		}else{
			if(scaling == 1){
				_name(x) = op.sprintf("_".str(x)."_L=%d,J=%.2f,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", x, J, J0, g, g_knot, h, w);	_key_title(x) = sprintf("L=%d",x);
				i0 = L0; iend = Lend; di = dL; 	out_dir = out_dir."size_scaling/"
				output_name = output_name.op.sprintf("_".str(site)."_J=%.2f,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f", J, J0, g, g_knot, h, w);
			} else{
				if(scaling == 2){
					_name(x) = op.sprintf("_".str(site)."_L=%d,J=%.2f,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", L, J, J0, 0.01*x, g_knot, h, w); 
					_key_title(x) = sprintf("g=%.2f",0.01*x)
					i0 = g0; iend = gend; di = dg; 	out_dir = out_dir."g_scaling/"
					output_name = output_name.op.sprintf("_".str(site)."_L=%d,J=%.2f,J0=%.2f,g0=%.2f,h=%.2f,w=%.2f", L, J, J0, g_knot, h, w);
				} else{
					if(scaling == 3){
						_name(x) = op."_".str(x).sprintf("_L=%d,J=%.2f,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", L, J, J0, g, g_knot, h, w);
						_key_title(x) = q_vs_j && operator < 2? sprintf("q/{/Symbol p}=%.2f", 2*x/(L+0.0))\
								: (operator > 1? sprintf("n=%d", x) :  sprintf("j=%d", x) ) 
						i0 = 0; iend = q_vs_j? L / 2 : L-1; di=1; if(operator > 1){ iend = 6;} 	
						out_dir = out_dir.(q_vs_j? "q" : "j")."_scaling/"
						output_name = output_name.op.sprintf("_L=%d,J=%.2f,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f", L, J, J0, g, g_knot, h, w);
					} else{
						if(scaling == 4){
							_dir(x) = 'realisation='.sprintf("%d/",x);
							_name(x) = _dir(x).op.sprintf("_".str(site)."_L=%d,J=%.2f,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", L, J, J0, g, g_knot, h, w);
							_key_title(x) = sprintf("r=%d",x);
							i0 = 0; iend = 9; di=1; 	out_dir = out_dir."realisation_scaling/"
							output_name = output_name.op.sprintf("_".str(site)."_L=%d,J=%.2f,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f", L, J, J0, g, g_knot, h, w);
						} else{
							_name(x) = op."_".str(site).sprintf("_L=%d,J=%.2f,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", L, J, J0, g, g_knot, h, w);
							i0=1; iend=1; di=1;
							output_name = output_name.op.sprintf("_".str(site)."_L=%d,J=%.2f,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f", L, J, J, J0, g, g_knot, h, w);
							_key_title(x) = sprintf("L=%d, g=%.2f, h=%.2f",L,g,h);
						}
					}
				}
			}
		}
	
	#------------------------------------ FIT
	__size=(iend - i0) / di + 1; 
	_tmp = (iend - i0) + di
	new_end = 0;
	name(x)=0; key_title(x)=0;
	this_dir = plot_exponent? dir_base.'ResponseFunction/DERIVATIVE/' : dir.(plot_normalized? 'NORMALIZED/' : '')
	if(two_panels){
		__size = 2*__size;
		new_end = iend + _tmp;
		str2(x) = str(x > iend? x - _tmp : x)
		that_dir = (plot_derivative? dir.'DERIVATIVE/' : dir_base.'ResponseFunction/')
		
		name(x) = x <= iend? this_dir.str2(x).'/'._name(x) : that_dir.str2(x).'/'.(smoothed? 'smoothed/' : '')._name(x - _tmp);
		key_title(x) = scaling == 5? (x<=iend? "I_A({/Symbol w})" : "S_A({/Symbol w})") : (x<=iend? _key_title(x) : _key_title(x - _tmp));
		
		if(use_png){set term pngcairo size 1900, 900	}
		else{ set term qt size 1500, 700;}
		set key inside bottom left spacing 1.1
	} else{
		new_end = iend;
		name(x) = this_dir.str(x).'/'._name(x); key_title(x) = _key_title(x);
	}
	
	f(w) = plot_normalized? a * atan(alfa * w) : a * atan(alfa * w)+b;
	fun1(w) = C * gamma / ((w - w01) * gamma**2 + 1.0)
	fun2(w) = D * delta / ((w - w02)**2 * delta**2 + 1.0)
	f_plot(a,b,alfa, w) = (w < x_min || w > x_max)? NaN : (plot_normalized? a * atan(alfa * w) : a * atan(alfa * w)+b);
	fun_plot(C, tau, w0, w, which) = which == 1? ( (w < x_min || w > x_max)? NaN : C * tau / ((w-w0) * tau**2 + 1.0) )\
												: ( (w < x_min || w > x_max)? NaN : C * tau / ( (tau * (w-w0))**2 + 1.0) )
	integrated_fun_plot1(aa, tau, w, w00) = 2 * aa / tau * log( (1 + (w-w00)*tau**2) / (1-(w+w00)*tau**2) )
	integrated_fun_plot2(aa, tau, w, w00)	= aa * ( atan( (w - w00) * tau) + atan( (w + w00) * tau));
	der_f_plot(a, alfa, w0, w) = a * alfa / ( alfa**2 * (w-w0)**2 + 1)
	
	#------------------------------------ STATS
		array wH[__size]; 
		array LTA[__size];
		array val_at_wH[__size];
		array val[__size];
		array a_list[__size]; array b_list[__size]; array alfa_list[__size]; array w0_list[__size]; array w01_list[__size]; array w02_list[__size];
		array C_list[__size]; array gamma_list[__size]; array D_list[__size]; array delta_list[__size]; 
		
		sub(x, i) = substract_LTA? abs(x - LTA[(i-i0)/di+1]) : abs(x);
		y_min = 1e5;
		print "par\tI(w=0), \twH, \tLTA, \t I(w=wH)"
		do for[i=i0:new_end:di]{
			idx = (i - i0) / di + 1
					val[idx] = 0; wH[idx] = 0;
					a_list[idx] = 0; b_list[idx] = 0; alfa_list[idx] = 0; w0_list[idx] = 0; w01_list[idx] = 0; w02_list[idx] = 0;
					C_list[idx] = 0; gamma_list[idx] = 0; D_list[idx] = 0; delta_list[idx] = 0;
			a = 1/pi; 	b = 0.0001;	alfa = 200;	w0 = 0.01; w01 = -1e-4; w02=-1e-1
			C = 1/pi;	D = 1/pi;	gamma=100; delta=100;
			name = name(i)
			#print name
			if(fileexist(name)){
				if(use_fit){
					if(!two_panels){
						fit[x_min:x_max][*:*] f(x) name u 1:2 via a, b, alfa; a_list[idx] = a; b_list[idx] = b; alfa_list[idx] = alfa;
					} else {
						if(i <= iend){
							if(plot_normalized) { fit[1e-5:0.2][*:*] f(x) name u 1:2 via a, alfa; a_list[idx] = a; b_list[idx] = 0.0; alfa_list[idx] = alfa;
							} else { fit[1e-5:0.2][*:*] f(x) name u 1:2 via a, b, alfa; a_list[idx] = a; b_list[idx] = b; alfa_list[idx] = alfa; }
						} else{
							#fit[x_min:x_max][*:*] fun1(x) name u 1:2 via C, gamma, w01; C_list[idx] = C; gamma_list[idx] = gamma; w01_list[idx] = w01;
							#fit[x_min:x_max][*:*] fun2(x) name u 1:2 via D, delta, w02; D_list[idx] = D; delta_list[idx] = delta; w02_list[idx] = w02;
							#print C,D,gamma,delta,w01,w02
						}
					}
				}
				stats name every ::0::1 using 2 nooutput prefix "s1";	val[idx] = s1_min; if(two_panels? (s1_min <= y_min && i <= iend) : (s1_min <= y_min)){ y_min = s1_min;}
				if(i <= iend || plot_derivative){
					stats name every ::0::1 using 3 nooutput prefix "s2"; 	wH[idx] = s2_min
					stats name every ::0::1 using 4 nooutput prefix "s3"; 	LTA[idx] = s3_min
					stats name using (abs($1 - s2_min)) nooutput prefix "s4"; 	idx_wH = s4_index_min
					stats name every ::idx_wH::(idx_wH+1) using 2 nooutput prefix "s5";	val_at_wH[idx] = sub(s5_min, i);
				} else {
					wH[idx] = wH[idx - __size / 2];	wH[idx] = LTA[idx - __size / 2];	LTA[idx] = wH[idx - __size / 2]; val_at_wH[idx] = val_at_wH[idx - __size / 2]	
				}
				
				if(use_fit) { print key_title(i),"  ", val[idx], wH[idx], "  ", LTA[idx], "  ", val_at_wH[idx], "\t\tfit:  ".sprintf("%.5f,  %.5f,  %.5f,  %.5f", a, b, alfa, w0)
				} else { print key_title(i),"  ", val[idx], wH[idx], "  ", LTA[idx], "  ", val_at_wH[idx]}
				#label_fit = label_fit.key_title(i).":\t".fstr(a,b,alfa)."\n\n"	  
			}
			else{
				print key_title(i)," -----------------------------------------------------------------------------"
			}		
		}
		if(!two_panels){ set ytics add(0.5); }
		if(y_min > 0.1){ set ytics add(sprintf("%.4f",y_min) y_min);}
		y_min = sub(y_min, iend)
		if(!y_log){ y_min = 0.0; }
		XRANGE = (rescale? "set xrange[15**nu*x_range_min:15**nu*1e2];" : "set xrange[x_range_min:1e1];" );
		YRANGE = plot_exponent? "set yrange[1e-10:1.5e3];"\
					: plot_normalized? "set yrange[0.5*y_min:1.0e0];" : "set yrange[0.9*y_min:1e0];"
		#if(fit) {set label 1 at 0.02,(y_log? 0.6:0.6) sprintf("%s",label_fit) front }
		

	#------------------------------------ PLOT
	if(use_png){
		if(two_panels){ out_dir = out_dir.'compare/'; };
		output_name = out_dir.output_name; 
		command = "mkdir ".out_dir; 
		system command
		set encoding utf8; 
		set output output_name.".png";
	}

		if(plot_normalized) { set key left top;
		} else {set key inside right bottom; }
		set xlabel (rescale? "{/Symbol w}* L^".sprintf("{%.2f}",nu) : '{/Symbol w}')
		set ylabel 'I_A({/Symbol w})' offset 4,2 rotate by 0
		
		MARGIN = !two_panels? "set lmargin at screen 0.10; set rmargin at screen 0.95; set bmargin at screen 0.10; set tmargin at screen 0.95;"\
			: "set lmargin at screen 0.10; set rmargin at screen 0.5; set bmargin at screen 0.12; set tmargin at screen 0.97;"
		MARGIN2= "set lmargin at screen 0.55; set rmargin at screen 0.97; set bmargin at screen 0.12; set tmargin at screen 0.97;"
		XRANGE2 = plot_derivative? "set xrange[1e-4:2e1];" : "set xrange[1e-1:1e1];"; 
		YRANGE2 = plot_derivative? "set yrange[1e-1:1e1];" : "set yrange[1e-8:1e-1];";
		
		set multiplot
		@MARGIN; @XRANGE; @YRANGE;
		if(scaling == 5){ set key at 1e-2,0.92 font ",23";}
		plot for[i=i0:iend:di] name(i) u ($1*(rescale? i**nu : 1.0)):(sub($2, i)) w l ls ((i-i0)/di+1) lw 1.5 title key_title(i),\
			name(iend) using 1:(0.5) w l dt (3,3) lc rgb "red" lw 1.5 notitle
		
		if(two_panels){
			set logscale y;
			set format y '10^{%L}'; @MARGIN2; @XRANGE2; @YRANGE2; unset ylabel; set key left bottom
			plot for[i=iend+di:new_end:di] name(i) u 1:(($2)) w l ls ((i-di-iend)/di+1) lw 1.5 title key_title(i)
		}
		set format y '%g'
		if(use_fit){
			if(scaling == 5){ set key at 1e-2,0.84 font ",23";}
			else{ set key left top font ",25";}
			@MARGIN; @UNSET;@XRANGE; @YRANGE; @SCALE; plot for[i=1:(iend-i0)/di + 1] f_plot(a_list[i], b_list[i], alfa_list[i], x) w l dt (3,5,10,5) lc rgb "black" lw 2.5 notitle
			#@MARGIN; @UNSET;@XRANGE; plot for[i=iend+1:new_end:di] D_list[(i-i0)/di+1] / norm*0.01 * ( atan( (x - w02_list[(i-i0)/di+1]) * delta_list[(i-i0)/di+1]) + atan( (x + w02_list[(i-i0)/di+1]) * delta_list[(i-i0)/di+1])) w l ls 4 lw 2.5 notitle
			if(0 && two_panels){
				norm = 3.5e-4; 
		@MARGIN2; @UNSET;@XRANGE2; @YRANGE2;  set key right top
			plot for[i=iend+di:new_end:di] fun_plot(C_list[(i-i0)/di+1], gamma_list[(i-i0)/di+1], w01_list[(i-i0)/di+1], x, 1) w l dt (8,8) lc rgb "blue" lw 2.5 title "1/{/Symbol w} fit"
		@MARGIN2; @UNSET;@XRANGE2; @YRANGE2; set key at graph 1,0.9
			plot for[i=iend+di:new_end:di] fun_plot(D_list[(i-i0)/di+1], delta_list[(i-i0)/di+1], w02_list[(i-i0)/di+1], x, 2) w l dt (1,1) lc rgb "green" lw 2.5 title "1/{/Symbol w}^2 fit"
		@MARGIN2; @UNSET;@XRANGE2; @YRANGE2; 
			plot for[i=1:(iend-i0)/di + 1] norm * der_f_plot(a_list[i], alfa_list[i], w0_list[i], x) w l dt (3,5,10,5) lc rgb "black" lw 1.5 lw 2.5 notitle
		@MARGIN; @UNSET; @XRANGE; @YRANGE; @SCALE; plot for[i=iend+1:new_end:di] integrated_fun_plot1(C_list[(i-i0)/di+1] / norm, gamma_list[(i-i0)/di+1], w01_list[(i-i0)/di+1], x) w l ls 3 lw 2.5 notitle
		@MARGIN; @UNSET; @XRANGE; @YRANGE; @SCALE; plot for[i=iend+1:new_end:di] integrated_fun_plot2(D_list[(i-i0)/di+1] / norm*0.1, delta_list[(i-i0)/di+1], w02_list[(i-i0)/di+1], x) w l ls 4 lw 2.5 notitle
			}
		}									
		
		if(!rescale){
			if(!substract_LTA && !plot_normalized){
				@MARGIN; @UNSET;@XRANGE; @YRANGE; @SCALE;
				plot for[i=i0:iend:di] name(i) u ($1 > wH[(i-i0)/di + 1]? NaN : $1):(LTA[(i-i0)/di + 1]) w l dt (1,1) lc rgb "green" lw 1.5 notitle
			}
			@MARGIN; @UNSET;@XRANGE; @YRANGE; @SCALE;
			#plot for[i=i0:iend:di] name(i) using (wH[(i-i0)/di + 1]):(val_at_wH[(i-i0)/di + 1]) w p dt (8,8) lc rgb "blue" lw 1.5 pt 4 ps 1.25 notitle
			plot wH using ($1 <= __size / 2? wH[$1] : NaN):(val_at_wH[$1]) w lp dt (8,8) lc rgb "blue" lw 1.5 pt 4 ps 1.25 notitle
		}
		unset multiplot
		if(!two_panels) {set term qt size 900, 800 font "Helvetica,12"}

