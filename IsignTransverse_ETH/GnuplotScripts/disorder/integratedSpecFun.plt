dir_base='../../results/disorder/PBC/'
dir = dir_base.'IntegratedResponseFunction/'
out_dir = 'Integrated_Spectral_Function/'
reset 

#------------------------------------ PREAMBLE
set autoscale
use_png = 0		# 1 if use png output, and 0 for qt output
if(use_png) { set term pngcairo size 1200, 1200 font sprintf("Helvetica,%d",16); }
else {set term qt size 900, 900 font sprintf("Helvetica,%d",14); }

set mxtics
set mytics

FORMAT = "set style line 12 lc rgb '#ddccdd' lt 1 lw 1.5; \
set style line 13 lc rgb '#ddccdd' lt 1 lw 0.5; \
set grid xtics mxtics ytics mytics back ls 12, ls 13;\
set size square;\
set xtics mirror;\
set ytics mirror;"
@FORMAT
set logscale xy
set format x "10^{%L}"

fileexist(name)=1#"[ -f name ] && echo 1 || echo 0"
#int(system("if exist \"".name."\" ( echo 1) else (echo 0)"))
	
UNSET = "unset tics; unset xlabel; unset ylabel; unset title; unset border;"

set style line 1 dt (3,5,10,5) lc rgb "black" lw 1.5
set style line 2 dt (3,3) lc rgb "red" lw 1.5
set style line 3 dt (8,8) lc rgb "blue" lw 1.5
set style line 4 dt (1,1) lc rgb "green" lw 1.5
set fit quiet
#------------------------------------ Margins for each row resp. column
TMARGIN = "set tmargin at screen 0.91; set bmargin at screen 0.56"
BMARGIN = "set tmargin at screen 0.46; set bmargin at screen 0.11"
# TICS
NOYTICS = "set format y '';"
YTICS = "set format y '%g';"

#------------------------------------ PARAMETERS
L = 14; 
g = 0.1
h = 0.3;
J0 = 0.; g_knot = 0.; 
w = 0.01;

x_range_min=1e-5

integrated_by_hand = 0 #integrated time evolution?
if(integrated_by_hand) cd '.\integrated'
rescale = 0				# rescale the spectral function by f(w, L)?
site = 0				# site at which the operator acts
scaling = 2				# size scaling=1 or h-scaling=0 or 	g-scaling=2	or 	q/j-scaling=3 or realisation=4 or user=5
q_vs_j = 0				# =1 - evolution of Sz_q, else ecol of Sz_j
operator = 1	 		# 1-SigmaZ , 0-Hq :local
two_panels = 1

rescale = 0
nu = 2		# power on L
if(scaling != 1) rescale = 0;
LIOM = 0				# plot LIOMs?
local = 0

	h0 = 20;	hend = 160;		dh = 20;
	g0 = 20;	gend = 90;		dg = 10;
	L0 = 10;	Lend = 15; 		dL = 1;

fit = 0
which_fit = 1						#=0 - a*exp(-t/b); =1 - a-b*log(t)
x_min = 1e-3; x_max = 2e-1;

op = operator? "SigmaZ" : "H";
str(x) = (q_vs_j? "q" : "j").sprintf("=%d",x);
#------------------------------------ GRAPHICS
	my_title = "{/*1.1 Integrated Spectral function I_A({/Symbol w}) with"
	if(scaling != 1) {my_title = my_title.sprintf(" L=%d", L);}
	if(scaling != 2) {my_title = my_title.sprintf(" g=%0.2f,", g);}
	if(scaling != 0) {my_title = my_title.sprintf(" h=%0.2f,", h);}
	tmp_title = my_title."}\n{/*1.1 for operator}\n\n{/*1.25 A = ";
	q_str = (site == -1? "L/2}" : ( site == 0? "0" : sprintf("%d{/Symbol p}/L}", 2.*site)))
	if(operator){
		if(scaling != 3) { my_title = tmp_title.(q_vs_j? "L^{-1/2}{/Symbol S}_j e^{iqj} {/Symbol s}^z_j; q=".q_str : sprintf("{/Symbol s}^z_{%d}}", site));}
		else {my_title = tmp_title.(q_vs_j? "L^{-1/2}{/Symbol S}_j e^{iqj} {/Symbol s}^z_j" : "{/Symbol s}^z_{j}");}
	} else{
		if(scaling != 3) { my_title = tmp_title.(q_vs_j? "L^{-1/2}{/Symbol S}_j cos(qj) H^j; q=".q_str : sprintf("H^{%d}}", site));}
		else {my_title = tmp_title.(q_vs_j? "L^{-1/2}{/Symbol S}_j cos(qj) H^j}" : "H_{j}}");}
	}
	#set title my_title

plot_dir = dir_base.'graph/IntegratedSpectralFunction/'

	
	
	#------------------------------------ DATA, FIT AND PLOT
	output_name = ""
	_name(x) = 0; _key_title(x) = 0;
	i0 = 0; iend = 0; di = 1;
		if(scaling == 0){
			_name(x) = str(site).'/'.op.sprintf("_".str(site)."_L=%d,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", L, J0, g, g_knot, 0.01*x, w);	_key_title(x) = sprintf("h=%.2f", x/100.)
			i0 = h0; iend = hend; di = dh; 	out_dir = out_dir."h_scaling/";
			output_name = output_name.op.sprintf("_".str(site)."_L=%d,J0=%.2f,g=%.2f,g0=%.2f,w=%.2f", L, J0, g, g_knot, w);
		}else{
			if(scaling == 1){
				__str(x) = site == -1 ? str(x/2) : str(site)
				_name(x) = __str(x).'/'.op.sprintf("_".__str(x)."_L=%d,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", x, J0, g, g_knot, h, w);	_key_title(x) = sprintf("L=%d",x);
				i0 = L0; iend = Lend; di = dL; 	out_dir = out_dir."size_scaling/"
				output_name = output_name.op.sprintf("_".str(site)."_J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f", J0, g, g_knot, h, w);
			} else{
				if(scaling == 2){
					_name(x) = str(site).'/'.op.sprintf("_".str(site)."_L=%d,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", L, J0, 0.01*x, g_knot, h, w); 
					_key_title(x) = sprintf("g=%.2f",0.01*x)
					i0 = g0; iend = gend; di = dg; 	out_dir = out_dir."g_scaling/"
					output_name = output_name.op.sprintf("_".str(site)."_L=%d,J0=%.2f,g0=%.2f,h=%.2f,w=%.2f", L, J0, g_knot, h, w);
				} else{
					if(scaling == 3){
						_name(x) = str(x).'/'.op."_".str(x).sprintf("_L=%d,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", L, J0, g, g_knot, h, w);
						_key_title(x) = (q_vs_j? sprintf("q/{/Symbol p}=%.2f", 2*x/(L+0.0)): sprintf("j=%d", x) )
						i0 = 0; iend = q_vs_j? L / 2 : L-1; di=1; 	out_dir = out_dir.(q_vs_j? "q" : "j")."_scaling/"
						output_name = output_name.op.sprintf("_L=%d,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f", L, J0, g, g_knot, h, w);
					} else{
						if(scaling == 4){
							_dir(x) = str(site).'/realisation='.sprintf("%d",x).'/';
							_name(x) = _dir(x).op.sprintf("_".str(site)."_L=%d,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", L, J0, g, g_knot, h, w);
							_key_title(x) = sprintf("r=%d",x);
							i0 = 0; iend = 9; di=1; 	out_dir = out_dir."realisation_scaling/"
							output_name = output_name.op.sprintf("_".str(site)."_L=%d,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f", L, J0, g, g_knot, h, w);
						} else{
							_name(x) = str(site).'/'.op."_".str(x).sprintf("_L=%d,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", L, J0, g, g_knot, h, w);
							i0=1; iend=1; di=1;
							output_name = output_name.op.sprintf("_".str(site)."_L=%d,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f", L, J0, g, g_knot, h, w);
						}
					}
				}
			}
		}
		
	__size=(iend - i0) / di + 1; _tmp = (iend - i0) + di
	new_end = 0;
	name(x)=0; key_title(x)=0;
	if(two_panels){
		__size = 2*__size;
		new_end = iend + _tmp;
		name(x) = x <= iend? dir._name(x) : dir_base.'ResponseFunction/'._name(x - _tmp);
		key_title(x) = scaling == 5? (x<=iend? "I_A({/Symbol w})" : "S_A({/Symbol w})") : (x<=iend? _key_title(x) : _key_title(x - _tmp));
		if(use_png){set term pngcairo size 1900, 900	}
		else{ set term qt size 1500, 700;}
		set key inside bottom left spacing 1.1
	} else{
		new_end = iend;
		name(x) = dir._name(x); key_title(x) = _key_title(x);
	}
	
	f(w) = a * atan(alfa * (w-w0))+b
	fun1(w) = C * gamma / ((w - w01) * gamma**2 + 1.0)
	fun2(w) = D * delta / ((w - w02)**2 * delta**2 + 1.0)
	f_plot(a,b,alfa, w0,w) = (w < 1e-5 || w > 0.2)? NaN : a*atan(alfa*(w-w0))+b;
	fun_plot(C, tau, w0, w, which) = which == 1? ( (w < x_min || w > x_max)? NaN : C * tau / ((w-w0) * tau**2 + 1.0) )\
												: ( (w < x_min || w > x_max)? NaN : C * tau / ( (tau * (w-w0))**2 + 1.0) )
	integrated_fun_plot1(aa, tau, w, w00) = 2 * aa / tau * log( (1 + (w-w00)*tau**2) / (1-(w+w00)*tau**2) )
	integrated_fun_plot2(aa, tau, w, w00)	= aa * ( atan( (w - w00) * tau) + atan( (w + w00) * tau));
	der_f_plot(a, alfa, w0, w) = a * alfa / ( alfa**2 * (w-w0)**2 + 1)
	
	y_min = 1e5;
	if(!LIOM){
		array wH[__size]; 
		array LTA[__size];
		array val[__size];
		array a_list[__size]; array b_list[__size]; array alfa_list[__size]; array w0_list[__size]; array w01_list[__size]; array w02_list[__size];
		array C_list[__size]; array gamma_list[__size]; array D_list[__size]; array delta_list[__size]; 
		print "par\tI(w=0),   \twH,	\tLTA,\t\t a\t b\t alfa\t w0"
		do for[i=i0:new_end:di]{
					val[(i-i0)/di+1] = 0; wH[(i-i0)/di+1] = 0;
					a_list[(i-i0)/di+1] = 0; b_list[(i-i0)/di+1] = 0; alfa_list[(i-i0)/di+1] = 0; w0_list[(i-i0)/di+1] = 0; w01_list[(i-i0)/di+1] = 0; w02_list[(i-i0)/di+1] = 0;
					C_list[(i-i0)/di+1] = 0; gamma_list[(i-i0)/di+1] = 0; D_list[(i-i0)/di+1] = 0; delta_list[(i-i0)/di+1] = 0;
			a = 1/pi; 	b = 0.0001;	alfa = 200;	w0 = 0.01; w01 = -1e-4; w02=-1e-1
			C = 1/pi;	D = 1/pi;	gamma=100; delta=100;
			name = name(i)
			print name
			if(fileexist(name)){
				if(fit){
					if(!two_panels){
						fit[1e-5:0.1][*:*] f(x) name u 1:2 via a, b, alfa, w0; a_list[(i-i0)/di+1] = a; b_list[(i-i0)/di+1] = b; alfa_list[(i-i0)/di+1] = alfa; w0_list[(i-i0)/di+1] = w0;
					}else {
						if(i <= iend){
							fit[1e-5:0.2][*:*] f(x) name u 1:2 via a, b, alfa, w0; a_list[(i-i0)/di+1] = a; b_list[(i-i0)/di+1] = b; alfa_list[(i-i0)/di+1] = alfa; w0_list[(i-i0)/di+1] = w0;
						}else{
							fit[x_min:x_max][*:*] fun1(x) name u 1:2 via C, gamma, w01; C_list[(i-i0)/di+1] = C; gamma_list[(i-i0)/di+1] = gamma; w01_list[(i-i0)/di+1] = w01;
							fit[x_min:x_max][*:*] fun2(x) name u 1:2 via D, delta, w02; D_list[(i-i0)/di+1] = D; delta_list[(i-i0)/di+1] = delta; w02_list[(i-i0)/di+1] = w02;
							print C,D,gamma,delta,w01,w02
						}
					}
				}
				stats name nooutput; n_cols = STATS_columns;
				stats name every ::0::1 using 2 nooutput prefix "s1";	val[(i-i0)/di+1] = s1_min; if(two_panels? (s1_min <= y_min && i <= iend) : (s1_min <= y_min)){ y_min = s1_min;}
				if(n_cols >= 3) {stats name every ::0::1 using 3 nooutput prefix "s2"; 	wH[(i-i0)/di+1] = s2_min} 
				else{ wH[(i-i0)/di+1]=0;}
				if(n_cols >= 4) {stats name every ::0::1 using 4 nooutput prefix "s3"; 	LTA[(i-i0)/di+1] = s3_min}
				else {LTA[(i-i0)/di+1] = 0};
				print key_title(i),"  ", val[(i-i0)/di+1], wH[(i-i0)/di+1], "  ", LTA[(i-i0)/di+1], "\tfit:  ".sprintf("%.5f,  %.5f,  %.5f,  %.5f", a, b, alfa, w0)
				#label_fit = label_fit.key_title(i).":\t".fstr(a,b,alfa)."\n\n"	  
			}
			else{
				print key_title(i)," -----------------------------------------------------------------------------"
			}		
		}
		set ytics add(0.5);
		if(y_min > 0.1){ set ytics add(y_min);}
		XRANGE = (rescale? "set xrange[15**nu*x_range_min:15**nu*1e2];" : "set xrange[x_range_min:1e2];" );
		YRANGE = "set yrange[0.9*y_min:1e0];"
		#if(fit) {set label 1 at 0.02,(y_log? 0.6:0.6) sprintf("%s",label_fit) front }
		

	#------------------------------------------ Controling output
	if(use_png){
		if(two_panels){ out_dir = out_dir.'compare/'; };
		output_name = out_dir.output_name; 
		command = "mkdir ".out_dir; 
		system command
		set encoding utf8; 
		set output output_name.".png";
	}

		set key inside right bottom #font ",16"
		set xlabel (rescale? "{/Symbol w}\267 L^".sprintf("{%.2f}",nu) : '{/Symbol w}')
		set ylabel 'I_A({/Symbol w})'
		
		MARGIN = !two_panels? "set lmargin at screen 0.10; set rmargin at screen 0.99; set bmargin at screen 0.10; set tmargin at screen 0.99;"\
			: "set lmargin at screen 0.10; set rmargin at screen 0.5; set bmargin at screen 0.12; set tmargin at screen 0.99;"
		MARGIN2= "set lmargin at screen 0.55; set rmargin at screen 0.99; set bmargin at screen 0.12; set tmargin at screen 0.99;"
		XRANGE2="set xrange[1e-3:1e1];"; YRANGE2="set yrange[1e-6:1e-1];";
		
		set multiplot
		@MARGIN; @XRANGE; @YRANGE;
		plot for[i=i0:iend:di] name(i) u ($1*(rescale? i**nu : 1.0)):2 w l lw 2.5 title key_title(i)
		
		if(two_panels){
			set format y '10^{%L}'; @MARGIN2; @XRANGE2; @YRANGE2; unset ylabel; set key left bottom
			plot for[i=iend+di:new_end:di] name(i) every 10 u ($1*(rescale? i**nu : 1.0)):2 w l lw 2.5 title key_title(i)
		}
		set format y '%g'
		if(fit){
			set key left top
			@MARGIN; @UNSET;@XRANGE; @YRANGE; plot for[i=1:(iend-i0)/di + 1] f_plot(a_list[i], b_list[i], alfa_list[i], w0_list[i], x) w l ls 1 lw 2.5 title "atan({/Symbol w}) fit"
			#@MARGIN; @UNSET;@XRANGE; plot for[i=iend+1:new_end:di] D_list[(i-i0)/di+1] / norm*0.01 * ( atan( (x - w02_list[(i-i0)/di+1]) * delta_list[(i-i0)/di+1]) + atan( (x + w02_list[(i-i0)/di+1]) * delta_list[(i-i0)/di+1])) w l ls 4 lw 2.5 notitle
			if(two_panels){
				norm = 3.5e-4; 
		@MARGIN2; @UNSET;@XRANGE2; @YRANGE2; set key right top
			plot for[i=iend+di:new_end:di] fun_plot(C_list[(i-i0)/di+1], gamma_list[(i-i0)/di+1], w01_list[(i-i0)/di+1], x, 1) w l ls 3 lw 2.5 title "1/{/Symbol w} fit"
		@MARGIN2; @UNSET;@XRANGE2; @YRANGE2; set key at graph 1,0.9
			plot for[i=iend+di:new_end:di] fun_plot(D_list[(i-i0)/di+1], delta_list[(i-i0)/di+1], w02_list[(i-i0)/di+1], x, 2) w l ls 4 lw 2.5 title "1/{/Symbol w}^2 fit"
		@MARGIN2; @UNSET;@XRANGE2; @YRANGE2; 
			plot for[i=1:(iend-i0)/di + 1] norm * der_f_plot(a_list[i], alfa_list[i], w0_list[i], x) w l ls 1 lw 2.5 notitle
		@MARGIN; @UNSET; @XRANGE; @YRANGE; plot for[i=iend+1:new_end:di] integrated_fun_plot1(C_list[(i-i0)/di+1] / norm, gamma_list[(i-i0)/di+1], w01_list[(i-i0)/di+1], x) w l ls 3 lw 2.5 notitle
		@MARGIN; @UNSET; @XRANGE; @YRANGE; plot for[i=iend+1:new_end:di] integrated_fun_plot2(D_list[(i-i0)/di+1] / norm*0.1, delta_list[(i-i0)/di+1], w02_list[(i-i0)/di+1], x) w l ls 4 lw 2.5 notitle
			}
		}									
		
		if(!two_panels && !rescale){
			@MARGIN; @UNSET;@XRANGE; @YRANGE;
			plot for[i=i0:iend:di] name(i) u ($1 > wH[(i-i0)/di + 1]? NaN : $1):(LTA[(i-i0)/di + 1]) w l ls 4 lw 1.5 notitle
		
			@MARGIN; @UNSET;@XRANGE; @YRANGE;
			plot wH using (wH[$1]):(LTA[$1]) w p ls 3 pt 4 ps 1.25 notitle
		}
		@MARGIN; @UNSET;@XRANGE; @YRANGE;
		plot name(iend) using 1:(0.5) w l ls 2 notitle
		unset multiplot
		if(!two_panels) {set term qt size 900, 800 font "Helvetica,12"}
	}
	else{






		set logscale y
		_which = 0 #0-orders, 1-sites
		n = 2
		set key outside right vertical maxcols(1)
		plotname = "timeEvolutionLIOM"
		tail = sprintf("_L=%d,g=%.2f,k=0,p=1,x=1,h=%.5f.dat", L, g, h)
		if(local){
			set title "TFIM LIOMS-densities \n\
				I_j^n = J(S^{xx}_{j,j+n}+S^{yy}_{j,j+n-2})+g(S^{xx}_{j,j+n-1}+S^{yy}_{j,j+n-1}) - n-even\n\n\
				I_j^n = S^{xy}_{j,j+n}-S^{yx}_{j,j+n} - n-odd,\n\n\
				where S^{{/Symbol a}{/Symbol b}}_{j,j+n}={/Symbol s}^{{/Symbol a}}_{j}{/Symbol P}_{k=1}^{n-1}{/Symbol s}^z_{j+k}{/Symbol s}^{{/Symbol b}}_{j+n}"
			s = ''
			do for [i=1:L-1] {s = s.(_which? sprintf(' %d,%d', n, i) : sprintf(' %d,%d', i, site)) } 
			plot for[w in s] plotname.w.tail w l title "n,j=".w
		}
		else{
			set title "TFIM LIOMS A=\n\
				I^n = {/Symbol S}_j J(S^{xx}_{j,j+n}+S^{yy}_{j,j+n-2}) + g(S^{xx}_{j,j+n-1}+S^{yy}_{j,j+n-1}) - n-even\n\n\
				I^n = {/Symbol S}_j S^{xy}_{j,j+n} - S^{yx}_{j,j+n} - n-odd,\n\n\
				where S^{{/Symbol a}{/Symbol b}}_{j,j+n}={/Symbol s}^{{/Symbol a}}_{j}{/Symbol P}_{k=1}^{n-1}{/Symbol s}^z_{j+k}{/Symbol s}^{{/Symbol b}}_{j+n}"
			array norm_arr[L];
			i_end = 10
			do for [i=0:i_end] {
				stats plotname.sprintf("%d", i).tail every ::1::2 using 2 nooutput
				norm_arr[i+1] = (norm ? STATS_min : 1.0) 
			}
			plot for[i=0:i_end] plotname.sprintf("%d", i).tail u 1:($2 / norm_arr[i+1])w l title sprintf("n=%d",i)
		}
	}

