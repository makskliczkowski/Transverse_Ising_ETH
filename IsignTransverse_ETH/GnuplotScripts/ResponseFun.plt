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
SET_LOG = x_log? "unset logscale xy; set logscale x; set format x '10^{%L}'; " : "set format x '%g'; "
SET_LOG = SET_LOG.(y_log? "set logscale y; set format y '10^{%L}'" : "set format y '%g'; ")
@SET_LOG;
SET_LINX = "unset logscale xy; set format x '%g'; set format y '%g';"
fileexist(name)=system("[ -f '".name."' ] && echo '1' || echo '0'") + 0	#int(system("if exist \"".name."\" (echo 1) else (echo 0)"))

UNSET = "unset tics; unset xlabel; unset ylabel; unset title; unset border;"

set fit quiet

# TICS
NOYTICS = "set format y '';"
YTICS = "set format y '%g';"

#------------------------------------ PARAMETERS
L = 10; 
g = 0.9; 
h = 0.8;
J=1.0
J0 = 0.; g_knot = 0.; 
w = 0.8;

SigX_or_SigZ = 1	 	# 0-SigX , 1-SigZ :local
operator_sum = 0		# is the operator a sum
site = 0				# site at which the operator acts
cor = 0					# correlations
scaling = 2				# size scaling=1 or h-scaling=0 or 	g-scaling=2	or 	q/j-scaling=3 or 4-realisations or 5-M scaling or 6-compare
q_vs_j = 0				# =1 - evolution of Sz_q, else ecol of Sz_j
operator = 1	 		# 1-SigmaZ , 0-Hq :local
compare = 0
smoothed = 0;
use_derivative = 0		# use derivative of integrated spectral function
if(use_derivative){ compare = 0; smoothed = 0;};
add_thouless_time = 0	# add thouless time from sff

LIOM = 0				# plot LIOMs?
local = 0

rescale=0				# rescale S_A by power law to find const region
add_line=1				# draw power-law: a/omega^n
a0=8e-5			# value of power-law plot at x=1
	h0 = 10;	hend = 90;		dh = 10;
	g0 = 10;	gend = 90;		dg = 10;
	L0 = 10;	Lend = 15; 		dL = 1;




op = ""
if(operator == 0) {op = "H"; }
if(operator == 1) {op = "SigmaZ";}
if(operator == 2) {op = "TFIM_LIOM_plus";}
if(operator == 3) {op = "TFIM_LIOM_minus";}

str(x) = (q_vs_j? "q" : "j").sprintf("=%d",x);
if(operator > 1){ str(x) = "n".sprintf("=%d",x); }
if(smoothed == 1){ op = 'smoothed/'.op; }
x_min_inset = 1e-2;
x_max_inset = 1e0;
y_min_inset = 1e-5
y_max_inset = 2e-2
RANGE2 = "set xrange[x_min_inset:x_max_inset]; set yrange[y_min_inset:y_max_inset];"
RANGE = (use_derivative? "set xrange[1e-5:3e2]; set yrange[1e-5:1e3];" : "set xrange[1e-3:2e1]; set yrange[1e-10:1e0];")
if(use_derivative){ a0 = 1e-3}
#rescale function
nu=2.0
fun(x, y, i) = scaling==1? x**nu*2**i*y : x**nu*2**L*y
	
dir_base='../results/ISING/disorder/PBC/'
dir = dir_base.(use_derivative? 'IntegratedResponseFunction/DERIVATIVE/' : 'ResponseFunction/DERIVATIVE/')
out_dir = 'Spectral_Function/'
	#------------------------------------ GRAPHICS
	set border back
	set key opaque inside right top 
	set xlabel "{/Symbol w}"
	set ylabel 'S_A({/Symbol w})'
	
	#------------------------------------ DATA, FIT AND PLOT
	output_name = ""
	name(x) = 0; key_title(x) = 0;
	_base(x, dis) = sprintf("_L=%d,J=%.2f,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", L, J, J0, g, g_knot, h, dis);
	name_th(x) = dir_base.'SpectralFormFactor/smoothed/'.sprintf("_L=%d,J=%.2f,J=%.2f,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", L, 1.0, J, J0, g, g_knot, h, 0.1);
	name2(x)=0
	array Mx[9];
	wH_fun(x) = 0;
	Mx[1] = 50; Mx[2] = 100; Mx[3] = 200; Mx[4] = 400; Mx[5] = 700; Mx[6] = 1000; Mx[7] = 2000; Mx[8] = 3000; Mx[9] = 4000;
	i0 = 0; iend = 0; di = 1;
		if(scaling == 0){
			wH_fun(x) = sqrt(L) / (0.341345 * 2**L) * sqrt(J*J + 1e-4*x*x + g*g + w*w / 3.)
			dir = dir.str(site).'/';
			_base(x, dis) = sprintf("_L=%d,J=%.2f,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", L, J, J0, g, g_knot, 0.01*x, dis);
			name_th(x) = dir_base.'SpectralFormFactor/smoothed/'.sprintf("_L=%d,J=%.2f,J=%.2f,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", L, 1.0, J, J0, g, g_knot, 0.01*x, 0.1);
			name(x) = dir.op."_".str(site)._base(x, w);	key_title(x) = sprintf("h=%.2f", x/100.)
			name2(x) = dir_base.'IntegratedResponseFunction/DERIVATIVE/'.str(site).'/'.op.sprintf("_".str(site)."_L=%d,J=%.2f,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", L, J, J0, g, g_knot, 0.01*x, w);
			i0 = h0; iend = hend; di = dh; 	out_dir = out_dir."h_scaling/";
			output_name = output_name.op.sprintf("_".str(site)."_L=%d,J=%.2f,J0=%.2f,g=%.2f,g0=%.2f,w=%.2f", L, J, J0, g, g_knot, w);
		}else{
			if(scaling == 1){
				wH_fun(x) = sqrt(x) / (0.341345 * 2**x) * sqrt(J*J + h*h + g*g + w*w / 3.)
				__str(x) = site == -1 ? str(x/2) : str(site)
				_base(x, dis) = sprintf("_L=%d,J=%.2f,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", x, J, J0, g, g_knot, h, dis);
				name_th(x) = dir_base.'SpectralFormFactor/smoothed/'.sprintf("_L=%d,J=%.2f,J=%.2f,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", x, 1.0, J, J0, g, g_knot, h, 0.1);
				name(x) = dir.__str(x).'/'.op."_".__str(x)._base(x, w);	key_title(x) = sprintf("L=%d",x);
				name2(x) = dir_base.'IntegratedResponseFunction/DERIVATIVE/'.__str(x).'/'.op.sprintf("_".__str(x)."_L=%d,J=%.2f,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", x, J, J0, g, g_knot, h, w);
				i0 = L0; iend = Lend; di = dL; 	out_dir = out_dir."size_scaling/"
				output_name = output_name.op.sprintf("_".str(site)."_J=%.2f,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f", J, J0, g, g_knot, h, w)
			} else{
				if(scaling == 2){
					wH_fun(x) = sqrt(L) / (0.341345 * 2**L) * sqrt(J*J + 1e-4*x*x + h*h + w*w / 3.)
					dir = dir.str(site).'/'
					_base(x, dis) = sprintf("_L=%d,J=%.2f,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", L, J, J0, 0.01*x, g_knot, h, dis);
					name_th(x) = dir_base.'SpectralFormFactor/smoothed/'.sprintf("_L=%d,J=%.2f,J=%.2f,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", L, 1.0, J, J0, 0.01*x, g_knot, h, 0.1);
					name(x) = dir.op."_".str(site)._base(x, w);
					name2(x) = dir_base.'IntegratedResponseFunction/DERIVATIVE/'.str(site).'/'.op.sprintf("_".str(site)."_L=%d,J=%.2f,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", L, J, J0, 0.01*x, g_knot, h, w);
					key_title(x) = sprintf("g=%.2f",0.01*x)
					i0 = g0; iend = gend; di = dg;	out_dir = out_dir."g_scaling/"
					output_name = output_name.op.sprintf("_".str(site)."_L=%d,J=%.2f,J0=%.2f,g0=%.2f,h=%.2f,w=%.2f", L, J, J0, g_knot, h, w);
				} else{
					if(scaling == 3){
						wH_fun(x) = sqrt(L) / (0.341345 * 2**L) * sqrt(J*J + h*h + g*g + w*w / 3.)
						name(x) = dir.str(x).'/'.op."_".str(x)._base(x, w);
						name2(x) = dir_base.'IntegratedResponseFunction/DERIVATIVE/'.str(x).'/'.op.sprintf("_".str(x)."_L=%d,J=%.2f,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", L, J0, g, g_knot, h, w);
						key_title(x) = q_vs_j && operator < 2? sprintf("q/{/Symbol p}=%.2f", 2*x/(L+0.0))\
								: (operator > 1? sprintf("n=%d", x) :  sprintf("j=%d", x) ) 
						i0 = 0; iend = q_vs_j? L / 2 : L-1; di=1; if(operator > 1){ iend = 6;} 
						out_dir = out_dir.(q_vs_j? "q" : "j")."_scaling/"
						output_name = output_name.op.sprintf("_L=%d,J=%.2f,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f", L, J0, g, g_knot, h, w);
					} else{
							if(scaling == 4){
							wH_fun(x) = sqrt(L) / (0.341345 * 2**L) * sqrt(J*J + h*h + g*g + w*w / 3.)
							_dir(x) = dir.str(site).'/realisation='.sprintf("%d",x).'/';
							name(x) = _dir(x).op."_".str(site)._base(x, w);
							key_title(x) = sprintf("r=%d",x);
							i0 = 0; iend = 9; di=1; 	out_dir = out_dir."realisation_scaling/"
							output_name = output_name.op.sprintf("_".str(site)."_L=%d,J=%.2f,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f", L, J0, g, g_knot, h, w);
						} else{
							if(scaling == 5){
								wH_fun(x) = sqrt(L) / (0.341345 * 2**L) * sqrt(J*J + h*h + g*g + w*w / 3.)
								name(x) = dir.str(site).'/'.op.sprintf("_".str(site)."_L=%d,J=%.2f,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f_M=%d.dat", L, J0, g, g_knot, h, w, x);
								key_title(x) = sprintf("M=%d",x);
								i0 = Mx[L-7]; iend = 5*Mx[L-7]; di=0.4*Mx[L-7]; 	out_dir = out_dir."M_scaling/"
								output_name = output_name.op.sprintf("_".str(site)."_L=%d,J=%.2f,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f", L, J, J0, g, g_knot, h, w);
							} else{
								
							}
						}
					}
				}
			}
		}
	__size = (iend-i0)/di+1
	array val1[__size];
	array val2[__size];
	array wH[__size];	array val_at_wH[__size];
	array w_tau[__size];	array val_at_tau[__size];
	
	do for[i=i0:iend:di]{
		idx = (i-i0)/di + 1
		wH[idx] = wH_fun(i);
		val1[idx] = 0.0; val2[idx] = 0.0;
		f(x, y) = abs(x - 0.01) < 1e-2? y : NaN
		if(fileexist(name(i))){ 
			stats name(i) every 5 using (f($1, $2)) nooutput prefix "stat1"; val1[idx] = stat1_max;
			stats name(i) using (abs($1 - wH[idx])) nooutput prefix "X"; 	idx_wH = X_index_min;
			stats name(i) every ::idx_wH::(idx_wH+1) using 2 nooutput prefix "Y";	val_at_wH[idx] = Y_min
		}
		if(add_thouless_time && fileexist(name(i))){ 
			GOE(x) = (x < 1? 2 * x - x*log(1+2*x) : 2-x*log( (2*x+1) / (2*x-1)))
			#stats name_th every ::0::1 using 3:4 nooutput;   w_tau[idx] = 1.0 / ( STATS_min_x * STATS_min_y );
			stats name_th(i) using ($1 > 2.5? NaN : ((log10( $2 / GOE($1) )) - 9e-2)**2) nooutput prefix "Y";       y_min = Y_index_min;
            stats name_th(i) using 1 every ::y_min::(y_min+1) nooutput;   w_tau[idx] = wH[idx] / STATS_min; 
			stats name(i) using (abs($1 - w_tau[idx])) nooutput prefix "X"; 	idx_tau = X_index_min;
			stats name(i) every ::idx_tau::(idx_tau+1) using 2 nooutput prefix "Y";	val_at_tau[idx] = Y_min
		} else { 
			w_tau[idx] = NaN ; val_at_tau[idx] = NaN;
		}
		if(compare){
			if(fileexist(name(i)) && fileexist(name2(i))){ 
				stats name2(i) using (f($1, $2)) nooutput prefix "stat2"; val2[idx] = stat2_max; 
				print val1[idx], val2[idx], val2[idx] / val1[idx], wH[idx]
			} 
		} else {
			print val1[idx], wH[idx], val_at_tau[idx], w_tau[idx]
		}
	}
	a = 1
	b = 0.1
	#------------------------------------------ Controling output
	if(use_png){
		if(compare){ out_dir = out_dir.'compare/'};
		output_name = out_dir.output_name; 
		command = "mkdir ".out_dir; 
		system command
		set encoding utf8; 
		set output output_name.".png";
	}
		MARGIN = "set lmargin at screen 0.10; set rmargin at screen 0.98; set bmargin at screen 0.10; set tmargin at screen 0.98;"
		MARGIN_INSET = "set lmargin at screen 0.20; set rmargin at screen 0.55; set bmargin at screen 0.15; set tmargin at screen 0.65;"
		set multiplot
		@MARGIN; @RANGE; @SET_LOG;
		plot for[i=i0:iend:di] name(i) every 2 u 1:2 w l ls ((i-i0)/di+1) lw 1.5 title key_title(i)
		if(compare){
			@MARGIN; @UNSET; @RANGE; @SET_LOG
			plot for[i=i0:iend:di] name2(i) u 1:($2* val1[(i-i0)/di+1] / val2[(i-i0)/di+1] ) w lp ls ((i-i0)/di+1) lw 1 ps 1.5 pt 6 title key_title(i)
		} else{
			if(rescale){
				@MARGIN_INSET; @RANGE2; unset xlabel; unset title; @SET_LINX;
				set ylabel sprintf("{/Symbol w}^{%.1f} D S_A({/Symbol w})",nu)
				plot for[i=i0:iend:di] name(i) u 1:(fun($1,$2,i)) w l lw 1.5 notitle
			}
			if(add_line){						
				@MARGIN; @UNSET; @RANGE; @SET_LOG
				plot name(i0) u 1:( ($1>x_min && $1 < x_max)? 10**(1-nu)* a0 /$1**nu : NaN) w l dt (3,5,10,5) lc rgb "black" lw 1.5 notitle
				set label 1 at 0.5*(x_max - x_min), 0.5*(y_max-y_min) sprintf("{/Symbol w}^{%.2f}",nu) front
			}
			@MARGIN; @UNSET; @RANGE; @SET_LOG
			plot wH using (wH[$1]):(val_at_wH[$1]) w lp dt (8,8) lc rgb "blue" lw 1.5 pt 4 ps 1.25 notitle
			if(add_thouless_time){
				plot w_tau using (w_tau[$1]):(val_at_tau[$1]) w lp dt (8,8) lc rgb "black" lw 1.5 pt 4 ps 1.25 notitle
			}
		}
		unset multiplot

	
	