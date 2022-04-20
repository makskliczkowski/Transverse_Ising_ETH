
dir_base='../results/disorder/PBC/'
dir = dir_base.'TimeEvolution/'
out_dir = 'Time_Evolution/'
reset 

#------------------------------------ PREAMBLE
set autoscale	
use_png = 0		# 1 if use png output, and 0 for qt output
if(use_png) { set term pngcairo size 1200, 1200 font sprintf("Helvetica,%d",22); }
else {set term qt size 900, 900 font sprintf("Helvetica,%d",22); }
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
L = 15; 
g = 0.2;
h = 0.8;
J0 = 0.; g_knot = 0.; 
w = 0.01;
rescale = 0				# rescale the spectral function by f(w, L)?
power = 0.5				# power in scaling with omega
operator = 1	 		# 1-SigmaZ , 0-Hq :local
site = 2				# site at which the operator acts
cor = 0					# correlations
scaling = 3				# size scaling=1 or h-scaling=0 or 	g-scaling=2	or 	q/j-scaling=3 or realisation-scaling=4 or 5-user defined
q_vs_j = 1				# =1 - evolution of Sz_q, else ecol of Sz_j
compare = 0
smoothed_data = 1		# plot smoothed data?
plot_exponent = 1		# plot exponent to find relaxation time

if(plot_exponent) smoothed_data = 0;
substract_LTA = 0
if(scaling == 4){ plot_exponent =0; smoothed_data = 0;}
rescale = 0
nu = 3		# power on L
#if(scaling != 1) rescale = 0;


LIOM = 0				# plot LIOMs?
local = 0

	h0 = 5;	hend = 20;		dh = 5;
	g0 = 20;	gend = 90;		dg = 10;
	L0 = 10;	Lend = 15; 		dL = 1;

use_fit = 0
which_fit = 1		# =1 -power-law || =0 -exp || =2-log

x_min = 1e1; x_max =2e2;  
# fit/rescale range
SHOW_FT_LABEL=1
if(use_fit==0) SHOW_FT_LABEL=0;

rescale_x(x, i) = x
if(scaling == 1){ rescale_x(x, i) = x / (rescale? (i**nu) : 1.0); }
if(scaling == 2){ rescale_x(x, i) = x * (rescale? ((0.01*i)**nu) : 1.0); }
#------------------------------------ GRAPHICS

if(compare){ set key inside bottom left}
set xlabel (rescale? sprintf("t*L^{%.2f}",nu) : 't')
set ylabel '<A(t)A>' rotate by 0 offset 4, 2

op = ""
if(operator == 0) {op = "H"; }
if(operator == 1) {op = "SigmaZ";}
if(operator == 2) {op = "TFIM_LIOM_plus";}
if(operator == 3) {op = "TFIM_LIOM_minus";}
if(smoothed_data) {op = 'smoothed/'.op;}
if(plot_exponent) {op = 'Exponent/'.op;}

str(x) = (q_vs_j? "q" : "j").sprintf("=%d",x);
if(operator > 1){ str(x) = "n".sprintf("=%d",x); }
#------------------------------------ DATA, FIT AND PLOT
#set yrange[0:1]
#plot dir.sprintf("q=1/SigmaZ_q=1_L=12,J0=0.00,g=%.2f,g0=0.00,h=0.80,w=0.01.dat", 0.0) u 1:2 notitle w lp ls 1,\
#	dir.sprintf("q=1/SigmaZ_q=1_L=12,J0=0.00,g=%.2f,g0=0.00,h=0.80,w=0.01.dat", 0.2) u 1:2 notitle w lp ls 2,\
#	dir.sprintf("q=1/SigmaZ_q=1_L=12,J0=0.00,g=%.2f,g0=0.00,h=0.80,w=0.01.dat", 1.1) u 1:2 notitle w lp ls 3
#
#exit;
output_name = ""
_name(x) = 0; _key_title(x) = 0;
i0 = 0; iend = 0; di = 1;
		if(scaling == 0){
			_name(x) = dir.str(site).'/'.op.sprintf("_".str(site)."_L=%d,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", L, J0, g, g_knot, 0.01*x, w);	_key_title(x) = sprintf("h=%.2f", x/100.)
			i0 = h0; iend = hend; di = dh; 	out_dir = out_dir."h_scaling/";
			output_name = output_name.op.sprintf("_".str(site)."_L=%d,J0=%.2f,g=%.2f,g0=%.2f,w=%.2f", L, J0, g, g_knot, w);
		}else{
			if(scaling == 1){
				__str(x) = site == -1 ? str(x/2) : str(site)
				_name(x) = dir.__str(x).'/'.op.sprintf("_".__str(x)."_L=%d,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", x, J0, g, g_knot, h, w);	_key_title(x) = sprintf("L=%d",x);
				i0 = L0; iend = Lend; di = dL; 	out_dir = out_dir."size_scaling/"
				output_name = output_name.op.sprintf("_".str(site)."_J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f", J0, g, g_knot, h, w);
			} else{
				if(scaling == 2){
					_name(x) = dir.str(site).'/'.op.sprintf("_".str(site)."_L=%d,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", L, J0, 0.01*x, g_knot, h, w); 
					_key_title(x) = sprintf("g=%.2f",0.01*x)
					i0 = g0; iend = gend; di = dg; 	out_dir = out_dir."g_scaling/"
					output_name = output_name.op.sprintf("_".str(site)."_L=%d,J0=%.2f,g0=%.2f,h=%.2f,w=%.2f", L, J0, g_knot, h, w);
				} else{
					if(scaling == 3){
						_name(x) = dir.str(x).'/'.op."_".str(x).sprintf("_L=%d,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", L, J0, g, g_knot, h, w);
						_key_title(x) = q_vs_j && operator < 2? sprintf("q/{/Symbol p}=%.2f", 2*x/(L+0.0))\
								: (operator > 1? sprintf("n=%d", x) :  sprintf("j=%d", x) ) 
						i0 = 0; iend = q_vs_j? L / 2 : L-1; di=1; if(operator > 1){ iend = 6;} 
						out_dir = out_dir.(q_vs_j? "q" : "j")."_scaling/"
						output_name = output_name.op.sprintf("_L=%d,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f", L, J0, g, g_knot, h, w);
					} else{
						if(scaling == 4){
							_dir(x) = dir.str(site).'/realisation='.sprintf("%d",x).'/';
							_name(x) = _dir(x).op.sprintf("_".str(site)."_L=%d,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", L, J0, g, g_knot, h, w);
							_key_title(x) = sprintf("r=%d",x);
							i0 = 0; iend = 9; di=1; 	out_dir = out_dir."realisation_scaling/"
							output_name = output_name.op.sprintf("_".str(site)."_L=%d,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f", L, J0, g, g_knot, h, w);
						} else{
							_name(x) = x==i0? dir.str(site).'/'.op.sprintf("_q=%d_L=%d,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", site, L, J0, g, g_knot, h, w)\
									: dir.sprintf("j=%d", x-1).'/'.op.sprintf("_j=%d_L=%d,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", x - 1, L, J0, g, g_knot, h, w);
							_key_title(x) = x==i0? "{/Symbol s}^z_{q=1}" : sprintf("{/Symbol s}_{j=L/2}", x-1)
							i0=0; iend = site; di=site
							out_dir = out_dir."compare_scales/"
							output_name = output_name.op.sprintf("_".str(site)."_L=%d,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f", L, J0, g, g_knot, h, w);
						}
					}
				}
			}
		}
		if(compare){
			set term qt size 1800, 600 font "Helvetica,12"
			if(use_png){set term pngcairo size 2700, 900	}
			else{ set term qt size 1800, 600;}
			set key inside bottom left spacing 1.1
		}
		set key inside bottom left
#------------------------ FIT functions
_a(a) = substract_LTA? 0.0 : a
f(x) = which_fit == 2? a - b*log(x+alfa) :\
		(substract_LTA? (which_fit == 1? b * x**(alfa) : b * exp(alfa*x)):\
							(which_fit == 1? a + b * x**(alfa) : a + b * exp(alfa*x)));

f_plot(a,b,alfa,x) = ((x < x_min || x > x_max)? NaN : (which_fit==2? a-b*log(x+alfa) : (which_fit==1? _a(a) + b*x**(alfa) : _a(a) + b * exp(alfa*x))))

fstr(a,b,alfa) = which_fit==1? sprintf("{\\~t^{%.2e}}",alfa) : sprintf("{\\~exp(-t / {/Symbol t}); {/Symbol t}=%.4f}",1./alfa)
if(which_fit==2) {fstr(a,b,alfa) = "\\~ln[x".(alfa<0?"-":"+").sprintf("%f]",abs(alfa))}

fileexist(name)=system("[ -f '".name."' ] && echo '1' || echo '0'") + 0	#int(system("if exist \"".name."\" (echo 1) else (echo 0)"))
label_fit = ""
a = 0.01
b = 1
alfa = -0.1
size = (iend - i0) / di+1

	array tH[size]
	array LTA[size]
	array val[size]
	array a_list[size]; array b_list[size]; array alfa_list[size]; #array x0_list[size];
	print "\nfit: a + b*x**alfa"
	print "par\tQ(t=0),   tH,    LTA,\t\t a\t b\t alfa"
	do for[i=i0:iend:di]{
		idx = (i-i0)/di+1
		name = _name(i)
		if(fileexist(name) && !plot_exponent){
			stats name every ::0::1 using 2 nooutput;	val[idx] = STATS_min
			stats name every ::0::1 using 3 nooutput; 	tH[idx]  = STATS_min;
			stats name every ::0::1 using 4 nooutput; 	LTA[idx] = STATS_min
			print _key_title(i),"  ", val[idx], tH[idx], LTA[idx], "\tfit:  ".sprintf("%.5f,  %.5f,  %.5f", a, b, alfa)
		}
		else{
			val[idx] = NaN; tH[idx] = NaN; LTA[idx] = NaN;
			if(use_fit){
				a_list[idx] = 0; b_list[idx] = 0; alfa_list[idx] = 0; #x0_list[idx] = 0;
			}
			print _key_title(i)," -----------------------------------------------------------------------------"
		}
	}
	sub(i,i0,di) = substract_LTA? LTA[(i-i0)/di+1] : 0.0;
	
		print "par\t\t a\t b\t alfa"
		y_min = 1e10
		do for[i=i0:iend:di]{
			idx = (i-i0)/di+1
			if( LTA[idx] < y_min){ y_min = LTA[idx]; }
			name = _name(i)	
			if(use_fit && fileexist(name)){
				if(which_fit == 2){
					fit[x_min:x_max][*:*] f(x) name u 1:($2-sub(i,i0,di)) via a, b, alfa;
					a_list[idx] = a; b_list[idx] = b; alfa_list[idx] = alfa;
				} else{
					if(substract_LTA){ fit[x_min:x_max][*:*] f(x) name u 1:($2-LTA[idx]) via b, alfa;}
					else{ fit[x_min:x_max][*:*] f(x) name u 1:2 via a, b, alfa; }
					a_list[idx] = a; b_list[idx] = b; alfa_list[idx] = alfa;
				}
				print _key_title(i), "\tfit:  ".sprintf("%.5f,  %.5f,  %.5f", a_list[idx], b_list[idx], alfa_list[idx])
				label_fit = label_fit."\t"._key_title(i).":\t".fstr(a,b,alfa)."\n\n"
			}
			else{
				a_list[idx] = 0; b_list[idx] = 0; alfa_list[idx] = 0; #x0_list[idx] = 0;
				print _key_title(i)," -----------------------------------------------------------------------------"
			}
		}
		if(y_min > 0.1 && y_log){ set ytics add(y_min);}
	RANGE = substract_LTA? "set xrange[1e0:1e4]; set yrange[1e-4:5e-1];" :\
						 (rescale? "set xrange[1e-3:9e3];" : "set xrange[2e-2:1e4];")."set yrange[".sprintf("%.5f", 0.75*y_min).":1.0];";
	#------------------------------------------ Controling output
	if(use_png){
		if(compare){ out_dir = out_dir.'compare_scales/'; };
		output_name = out_dir.output_name; 
		command = "mkdir ".out_dir; 
		system command
		set encoding utf8; 
		set output output_name.".png";
	}
	
	if(SHOW_FT_LABEL) {set label 1 at 0.01,(y_log? 0.25:0.6) sprintf("%s",label_fit) front }
	MARGIN = !compare? "set lmargin at screen 0.10; set rmargin at screen 0.95; set bmargin at screen 0.10; set tmargin at screen 0.99;"\
			: "set lmargin at screen 0.07; set rmargin at screen 0.37; set bmargin at screen 0.12; set tmargin at screen 0.97;"
	MARGIN2= "set lmargin at screen 0.38; set rmargin at screen 0.68; set bmargin at screen 0.12; set tmargin at screen 0.97;"
	MARGIN3= "set lmargin at screen 0.69; set rmargin at screen 0.99; set bmargin at screen 0.12; set tmargin at screen 0.97;"

	var=0.3
	set multiplot
		if(compare){ set logscale xy;}
		else{ @SCALE; }
		@MARGIN; @RANGE;
		#set yrange[1e-1:1e6] 
		plot for[i=i0:iend:di] _name(i) u (rescale_x($1, i)):(abs($2-sub(i,i0,di))) w l lw 1.5 title _key_title(i)#, (var+exp(-x/85.7258+2))/(1+var) w l ls 4 t "e^{-t/{/Symbol t}}"
		if(compare){
			unset label 1;
			unset logscale x; set format x '%g'; set logscale y; unset ylabel; @MARGIN2; @RANGE; plot for[i=i0:iend:di] _name(i) u (rescale_x($1, i)):($2-sub(i,i0,di)) w l lw 1.5 notitle
			#set label 1 at 0.012,0.4 sprintf("%s",label_fit) front
			unset logscale y; set format y '%g'; set logscale x; unset ylabel; @MARGIN3; @RANGE; plot for[i=i0:iend:di] _name(i) u (rescale_x($1, i)):($2-sub(i,i0,di)) w l lw 1.5 notitle
		}
		#ssf_name(gx,Lx) = sprintf("./Spectracd Pro	lFormFactor/_L=%d,J0=0.00,g=%.2f,g0=0.00,h=0.80,w=0.01.dat", Lx, gx)
		#dim=2.0**L
		#scale=dim / (dim*dim-dim)
		#print scale
		#@MARGIN; @RANGE; plot for[i=i0:iend:di] ssf_name(0.01*i, L) u 1:( scale*($2 - 1.0) + LTA[(i-i0)/di+1] * ( 1.0 + scale - scale*$2) ) w p pt 6 lw 0.5 notitle		
		if(use_fit){
			@UNSET; @MARGIN; @RANGE; 
			if(compare){ set logscale xy;}
			else{ @SCALE; }
			plot for[i=i0:iend:di] f_plot(a_list[(i-i0)/di + 1], b_list[(i-i0)/di + 1], alfa_list[(i-i0)/di + 1], x) w l ls 1 notitle
			if(compare){
				@UNSET; @MARGIN2; @RANGE; unset logscale x; set format x '%g'; set logscale y; plot for[i=i0:iend:di] f_plot(a_list[(i-i0)/di + 1], b_list[(i-i0)/di + 1], alfa_list[(i-i0)/di + 1], x) w l ls 1 notitle
				@UNSET; @MARGIN3; @RANGE; unset logscale y; set format y '%g'; set logscale x; plot for[i=i0:iend:di] f_plot(a_list[(i-i0)/di + 1], b_list[(i-i0)/di + 1], alfa_list[(i-i0)/di + 1], x) w l ls 1 notitle
			}
		}	
		if(!rescale && !substract_LTA && !compare){
			print "am here"
			@UNSET; @MARGIN; @RANGE; @SCALE;
			plot for[i=i0:iend:di] _name(i) u ($1 < tH[(i-i0)/di + 1]? NaN : $1):(LTA[(i-i0)/di + 1]) w l ls 2 notitle
			
			@UNSET; @MARGIN; @RANGE; @SCALE;
			plot tH using (tH[$1]):(LTA[$1]) w lp ls 3 ps 0.75 notitle
		}
	unset multiplot
	