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
SET_LOG = "set logscale xy; set format x '10^{%L}'; set format y '10^{%L}';"
@SET_LOG;
SET_LINX = "unset logscale xy; set format x '%g'; set format y '%g';"
fileexist(name)=1#int(system("if exist \"".name."\" ( echo 1) else (echo 0)"))

UNSET = "unset tics; unset xlabel; unset ylabel; unset title; unset border;"

set style line 1 dt (3,5,10,5) lc rgb "black" lw 1.5
set style line 2 dt (3,3) lc rgb "red" lw 1.5
set style line 3 dt (8,8) lc rgb "blue" lw 1.5
set style line 4 dt (1,1) lc rgb "green" lw 1.5
set fit quiet

# TICS
NOYTICS = "set format y '';"
YTICS = "set format y '%g';"

#------------------------------------ PARAMETERS
L = 15; 
g = 0.4; 
h = 0.8;
J0 = 0.; g_knot = 0.; 
w = 0.01;

SigX_or_SigZ = 1	 	# 0-SigX , 1-SigZ :local
operator_sum = 0		# is the operator a sum
site = 1				# site at which the operator acts
cor = 0					# correlations
scaling = 1				# size scaling=1 or h-scaling=0 or 	g-scaling=2	or 	q/j-scaling=3 or 4-realisations or 5-M scaling or 6-compare
q_vs_j = 1				# =1 - evolution of Sz_q, else ecol of Sz_j
operator = 1	 		# 1-SigmaZ , 0-Hq :local
compare = 0
use_derivative = 1		# use derivative of integrated spectral function
if(use_derivative){ compare = 0};

LIOM = 0				# plot LIOMs?
local = 0

rescale=0				# rescale S_A by power law to find const region
add_line=0				# draw power-law: a/omega^n
a0=8e-5					# value of power-law plot at x=1
	h0 = 20;	hend = 120;		dh = 20;
	g0 = 5;	gend = 40;		dg = 5;
	L0 = 11;	Lend = 15; 		dL = 1;



op = operator? "SigmaZ" : "H";
str(x) = (q_vs_j? "q" : "j").sprintf("=%d",x);

x_min = 1e-2;
x_max = 1e-1;
y_min = 1e-1
y_max = 2e0
RANGE2 = "set xrange[x_min:x_max]; set yrange[y_min:y_max];"
RANGE = (use_derivative? "set xrange[1e-3:3e1]; set yrange[1e-4:1e2];" : "set xrange[1e-3:3e1]; set yrange[1e-10:0.05];")
if(use_derivative){ a0 = 1e-3}
#rescale function
nu=2.0
fun(x, y, i) = scaling==1? x**nu*2**i*y : x**nu*2**L*y
	
dir_base='../results/disorder/PBC/'
dir = dir_base.(use_derivative? 'IntegratedResponseFunction/DERIVATIVE/' : 'ResponseFunction/')
out_dir = 'Spectral_Function/'
	#------------------------------------ GRAPHICS
	set key inside right top
	set xlabel "{/Symbol w}"
	set ylabel 'S_A({/Symbol w})'
	
	#------------------------------------ DATA, FIT AND PLOT
	output_name = ""
	name(x) = 0; key_title(x) = 0;
	name2(x)=0
	array Mx[9];
	Mx[1] = 50; Mx[2] = 100; Mx[3] = 200; Mx[4] = 400; Mx[5] = 700; Mx[6] = 1000; Mx[7] = 2000; Mx[8] = 3000; Mx[9] = 4000;
	i0 = 0; iend = 0; di = 1;
		if(scaling == 0){
			dir = dir.str(site).'/';
			name(x) = dir.op.sprintf("_".str(site)."_L=%d,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", L, J0, g, g_knot, 0.01*x, w);	key_title(x) = sprintf("h=%.2f", x/100.)
			name2(x) = dir_base.'IntegratedResponseFunction/DERIVATIVE/'.str(site).'/'.op.sprintf("_".str(site)."_L=%d,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", L, J0, g, g_knot, 0.01*x, w);
			i0 = h0; iend = hend; di = dh; 	out_dir = out_dir."h_scaling/";
			output_name = output_name.op.sprintf("_".str(site)."_L=%d,J0=%.2f,g=%.2f,g0=%.2f,w=%.2f", L, J0, g, g_knot, w);
		}else{
			if(scaling == 1){
				__str(x) = site == -1 ? str(x/2) : str(site)
				name(x) = dir.__str(x).'/'.op.sprintf("_".__str(x)."_L=%d,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", x, J0, g, g_knot, h, w);	key_title(x) = sprintf("L=%d",x);
				name2(x) = dir_base.'IntegratedResponseFunction/DERIVATIVE/'.__str(x).'/'.op.sprintf("_".__str(x)."_L=%d,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", x, J0, g, g_knot, h, w);
				i0 = L0; iend = Lend; di = dL; 	out_dir = out_dir."size_scaling/"
				output_name = output_name.op.sprintf("_".str(site)."_J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f", J0, g, g_knot, h, w)
			} else{
				if(scaling == 2){
					dir = dir.str(site).'/'
					name(x) = dir.op.sprintf("_".str(site)."_L=%d,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", L, J0, 0.01*x, g_knot, h, w);
					name2(x) = dir_base.'IntegratedResponseFunction/DERIVATIVE/'.str(site).'/'.op.sprintf("_".str(site)."_L=%d,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", L, J0, 0.01*x, g_knot, h, w);
					key_title(x) = sprintf("g=%.2f",0.01*x)
					i0 = g0; iend = gend; di = dg;	out_dir = out_dir."g_scaling/"
					output_name = output_name.op.sprintf("_".str(site)."_L=%d,J0=%.2f,g0=%.2f,h=%.2f,w=%.2f", L, J0, g_knot, h, w);
				} else{
					if(scaling == 3){
						name(x) = dir.str(x).'/'.op."_".str(x).sprintf("_L=%d,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", L, J0, g, g_knot, h, w);
						name2(x) = dir_base.'IntegratedResponseFunction/DERIVATIVE/'.str(x).'/'.op.sprintf("_".str(x)."_L=%d,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", L, J0, g, g_knot, h, w);
						key_title(x) = (q_vs_j? sprintf("q/{/Symbol p}=%.2f", 2*x/(L+0.0)): sprintf("j=%d", x) )
						i0 = 0; iend = q_vs_j? L / 2 : L-1; di=1; 	out_dir = out_dir.(q_vs_j? "q" : "j")."_scaling/"
						output_name = output_name.op.sprintf("_L=%d,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f", L, J0, g, g_knot, h, w);
					} else{
							if(scaling == 4){
							_dir(x) = dir.str(site).'/realisation='.sprintf("%d",x).'/';
							name(x) = _dir(x).op.sprintf("_".str(site)."_L=%d,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", L, J0, g, g_knot, h, w);
							key_title(x) = sprintf("r=%d",x);
							i0 = 0; iend = 9; di=1; 	out_dir = out_dir."realisation_scaling/"
							output_name = output_name.op.sprintf("_".str(site)."_L=%d,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f", L, J0, g, g_knot, h, w);
						} else{
							if(scaling == 5){
								name(x) = dir.str(site).'/'.op.sprintf("_".str(site)."_L=%d,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f_M=%d.dat", L, J0, g, g_knot, h, w, x);
								key_title(x) = sprintf("M=%d",x);
								i0 = Mx[L-7]; iend = 5*Mx[L-7]; di=0.4*Mx[L-7]; 	out_dir = out_dir."M_scaling/"
								output_name = output_name.op.sprintf("_".str(site)."_L=%d,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f", L, J0, g, g_knot, h, w);
							} else{
								
							}
						}
					}
				}
			}
		}
	
	array val1[(iend-i0)/di+1];
	array val2[(iend-i0)/di+1];
	do for[i=i0:iend:di]{
		if(fileexist(name(i))){ stats name(i) every ::0::31 using 2 nooutput prefix "stat1"; val1[(i-i0)/di+1] = stat1_max; }
		if(compare){
			if(fileexist(name2(i))){ 
				stats name2(i) every ::0::31 using 2 nooutput prefix "stat2"; val2[(i-i0)/di+1] = stat2_max; 
				print val1[(i-i0)/di+1], val2[(i-i0)/di+1], val2[(i-i0)/di+1] / val1[(i-i0)/di+1] 
			}
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
	if(!LIOM){
		MARGIN = "set lmargin at screen 0.10; set rmargin at screen 0.99; set bmargin at screen 0.10; set tmargin at screen 0.99;"
		MARGIN_INSET = "set lmargin at screen 0.20; set rmargin at screen 0.55; set bmargin at screen 0.15; set tmargin at screen 0.65;"
		set multiplot
		@MARGIN; @RANGE; @SET_LOG;
		plot for[i=i0:iend:di] name(i) u 1:2 w l lw 1.5 title key_title(i)	
		if(compare){
			@MARGIN; @UNSET; @RANGE; @SET_LOG
			plot for[i=i0:iend:di] name2(i) u 1:($2* val1[(i-i0)/di+1] / val2[(i-i0)/di+1] ) w p ps 1.5 pt 6 title key_title(i)
		} else{
			if(rescale){
				@MARGIN_INSET; @RANGE2; unset xlabel; unset title; @SET_LINX;
				set ylabel sprintf("{/Symbol w}^{%.1f}\267D\267S_A({/Symbol w})",nu)
				plot for[i=i0:iend:di] name(i) u 1:(fun($1,$2,i)) w l lw 1.5 notitle
			}
			if(add_line){						
				@MARGIN; @UNSET; @RANGE; @SET_LOG
				plot name(iend) u 1:( ($1>x_min && $1 < x_max)? 10**(1-nu)* a0 /$1**nu : NaN) w l ls 1 notitle
				set label 1 at 0.5*(x_max - x_min), 0.5*(y_max-y_min) sprintf("{/Symbol w}^{%.2f}",nu) front
			}
			#plot sample [i=1:(iend-i0)/di + 1] '+' using (wH[i]):(2*val[i]) w l ls 3 notitle
		}
		unset multiplot
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
	
	
	
	