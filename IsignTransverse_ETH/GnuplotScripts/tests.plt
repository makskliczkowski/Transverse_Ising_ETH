reset 
##--PREAMBLE
set autoscale
set term qt size 1050, 900 font "Helvetica,18"
set mxtics
set mytics
set style line 12 lc rgb '#ddccdd' lt 1 lw 1.5
set style line 13 lc rgb '#ddccdd' lt 1 lw 0.5
set grid xtics mxtics ytics mytics back ls 12, ls 13
set size square
set xtics nomirror
set key inside bottom right font ",16" spacing 2

set style line 1 dt (3,5,10,5) lc rgb "black" lw 1.5
set style line 2 dt (3,3) lc rgb "red" lw 1.5
set style line 3 dt (8,8) lc rgb "blue" lw 1.5
set style line 4 dt (1,1) lc rgb "green" lw 1.5

fileexist(name)=system("[ -f '".name."' ] && echo '1' || echo '0'") + 0#int(system("if exist \"".name."\" ( echo 1) else (echo 0)"))

# Margins for each row resp. column
TMARGIN = "set tmargin at screen 0.99; set bmargin at screen 0.56"
BMARGIN = "set tmargin at screen 0.46; set bmargin at screen 0.10"
LMARGIN = "set lmargin at screen 0.10; set rmargin at screen 0.46"
RMARGIN = "set lmargin at screen 0.46; set rmargin at screen 0.82"
RANGE = "set xrange[0:1]; set yrange[0:2.0]"
UNSET = "unset tics; unset xlabel; unset ylabel; unset title; unset key; unset border;"
#-- PARAMETERS
w = 0.01
g = 0.4
L = 14
h = 1.0
h0 = 80
hend = 400
dh = 40
E_or_O = 0  	#energies(0) or operator(1)
var = 1			# variance?
smoothen = 1	# smooth data?
scaling = 2		# size scaling=1 or h-scaling=0 or 	g-scaling=2	or 	q/j-scaling=3 or realisation-scaling=4 or 5-user defined
y_ax = 2		# =0 - ||O_diag||; =1 - ||O_off||; =2 - <r>
model = 0

choose = 2		# 0 - qpt / 1 - magnetization quench / 2 - ssf / 3 - 
#-- GRAPHICS

#set y2tics 0.37, 0.05
#set y2range[0.37:0.54]
set xrange[0.1:5.0]
set ytics nomirror
set autoscale xy
#set autoscale y2

set lmargin at screen 0.13; set rmargin at screen 0.96;
set tmargin at screen 0.96; set bmargin at screen 0.10
dir_base='../results/'.(model? "symmetries" : "disorder").'/PBC/'
if(choose == 0){
	delta = 0.025
	# integral_f(x) takes one variable, the upper limit.  0 is the lower limit. calculate the integral of function f(t) from 0 to x
	# choose a step size no larger than delta such that an integral number of steps will cover the range of integration.

	f(x,gx) = (1. - cos(x)) / (sqrt(1+gx**2+2*gx*cos(x)))
	integral_f(x, gx) = (x>0)?int1a(x,gx,x/ceil(x/delta)):-int1b(x,gx,-x/ceil(-x/delta))
	int1a(x,gx,d) = (x<=d*.1) ? 0 : (int1a(x-d,gx,d)+(f(x-d,gx)+4*f(x-d*.5,gx)+f(x,gx))*d/6.)
	int1b(x,gx,d) = (x>=-d*.1) ? 0 : (int1b(x+d,gx,d)+(f(x+d,gx)+4*f(x+d*.5,gx)+f(x,gx))*d/6.)
	#
	# integral2_f(x,y) takes two variables; x is the lower limit, and y the upper.
	# calculate the integral of function f(t) from x to y

	integral2_f(x,y) = (x<y)?int2(x,y,(y-x)/ceil((y-x)/delta)): \
	                        -int2(y,x,(x-y)/ceil((x-y)/delta))
	int2(x,y,d) = (x>y-d*.5) ? 0 : (int2(x+d,y,d) + (f(x)+4*f(x+d*.5)+f(x+d))*d/6.)
	mx(gx) = gx < 1? 0 : (1. - 1. / gx ** 2)**(1./8.)
	set xlabel 'g / J'
	set ylabel '<{/Symbol s}^z>' rotate by 0
	dir = dir_base.'Magnetization/'
	#_name(Lx) = dir.sprintf("IsingQPT_L=%d,h=%.2f,k=0,p=1,x=1.dat", Lx, h);
	_name(Lx) = dir.sprintf("IsingQPT_L=%d,J0=0.00,g0=0.00,h=%.2f,w=%.2f.dat", Lx, h, w);
	set xrange[0:2]
	set yrange[0:1]
	plot for [Lx=10:26:2] _name(Lx) using 1:3 w lp title sprintf("L=%d", Lx), mx(x)#integral_f(pi, 1. / x) / pi w l ls 1
} else{
	if(choose == 1){
		#_name = dir.sprintf("quench_g_init=0.00,h_init-0.50_L=%d,g=%.2f,h=%.2f,k=0,p=1,x=1.dat", L, g, h)
		eh=0
		nejm = eh? "quench_ferromagnet" : "quench_domain_wall";
		_name = dir.nejm.sprintf("_L=%d,J0=0.00,g=%.2f,g0=0.00,h=%.2f,w=%.2f.dat", L, g, h, w);
		#plot for[gx=20:60:10] dir.sprintf("meson_coverage_L=%d,J0=0.00,g=%.2f,g0=0.00,h=%.2f,w=%.2f.dat", L, gx*0.01, h, w) u 1:2 w lp title sprintf("g=%.2f", 0.01*gx)
		#exit;
		 heatmap = 0
		if(heatmap){
			set xrange[-0.5:L-0.5]
			set yrange[0:50]
			set view map
			set pm3d interpolate 0,0
			#set cbrange[0:0.1]
			set yrange
			splot _name u 1:2:3 with image
		} else {
			set key outside right
			set xrange[0:60];
			set logscale x
			i = 0;
			plot _name u ($1 == 1? $2 : NaN):3 w lp title 'i=1',\
				_name u ($1 == 2? $2 : NaN):3 w lp title 'i=2',\
				_name u ($1 == L / 2? $2 : NaN):3 w lp title 'i=L/2',\
				_name u ($1 == L / 2 + 1? $2 : NaN):3 w lp title 'i=L/2+1'
		}

	} else{
		if(choose == 2){
			set xlabel 't / t_H'
			set logscale xy
			set format x '10^{%L}'
			set format y '10^{%L}'
			dir = dir_base.'SpectralFormFactor/'
			if(smoothen){ dir = dir.'smoothed/';}
			i0=0; iend=1; di=1;
			ssf_name(x) = ""; key_title(x) = "";
			if(scaling == 0){	# h - sclaing
				i0=20; iend=120; di=10;
				ssf_name(x) = dir.sprintf("_L=%d,J0=0.00,g=%.2f,g0=0.00,h=%.2f,w=0.01.dat", L, g, 0.01 * x);
				key_title(x) = sprintf("h = %.2f", 0.01 * x);
			} else {	
				if(scaling == 1){	# L - scaling
					i0=10; iend=14; di=1;
					ssf_name(x) = dir.sprintf("_L=%d,J0=0.00,g=%.2f,g0=0.00,h=%.2f,w=0.01.dat", x, g, h);
					key_title(x) = sprintf("L = %d", x);
				} else {	#	g - scaling
					i0=10; iend=110; di=20;
					ssf_name(x) = dir.sprintf("_L=%d,J0=0.00,g=%.2f,g0=0.00,h=%.2f,w=0.01.dat", L, 0.01 * x, h);
					key_title(x) = sprintf("g = %.2f", 0.01 * x);
				}
			}
			set key right top
			set xrange[1e-5:8];
			set yrange[2e-3:4e4]
				array tH[(iend-i0)/di+1]
				array val[(iend-i0)/di+1]
			do for[i=i0:iend:di]{
					_name = ssf_name(i)
					if(fileexist(_name)){
						stats _name every ::0::0 using 2 nooutput;	val[(i-i0)/di+1] = STATS_min
						stats _name every ::0::0 using 3 nooutput; 	tH[(i-i0)/di+1] = STATS_min;
						print key_title(i),"  ", val[(i-i0)/di+1], tH[(i-i0)/di+1]
					}
					else{
						print key_title(i)," -----------------------------------------------------------------------------"
					}
			}
			plot for[i=i0:iend:di] ssf_name(i) u 1:2 w l lw 1.5 t key_title(i), (x < 1? 2 * x - x*log(1+2*x) : 2-x*log( (2*x+1) / (2*x-1))) t 'GOE'

		} else{
			#set logscale x; set xrange[*:150]
			dir = dir_base.'Entropy/'
			set ylabel 'S(t)'
			set xlabel 't'
			ssf_name(gx,Lx,x,M) = dir.sprintf("TimeEvolution_L=%d,J0=0.00,g=%.2f,g0=0.00,h=0.80,w=0.01,x=%.4f,M=%d.dat", Lx, gx, x, M)
			M=6
			dt=0.16
			L=10
			M_list = '2 3 4 5 6 7 8 9'
			dt_list = '1e-2 2e-2 4e-2 8e-2 16e-2 64e-2 128e-2 256e-2 512e-2'
			MARGIN = "set lmargin at screen 0.10; set rmargin at screen 0.99; set bmargin at screen 0.10; set tmargin at screen 0.99;"
			set xrange[1e-3:1e2]; set yrange[1e-6:4.0]
			set key right bottom
			#plot for[Mx in M_list] ssf_name(0.6, L, dt, 1.*Mx) u 1:3 w lp pt 6 ps 0.5 t sprintf("x=%.3f, M=%d", dt, 1.*Mx), ssf_name(0.6, L, dt, M) u 1:2 w lp ls 1 title 'ED'
			plot for[dt in dt_list] ssf_name(0.6, L, 1.*dt, M) u 1:3 w lp pt 6 ps 1.5 t sprintf("x*{/Symbol w}_{max}=%.3f, M=%d", 1.*dt, M), ssf_name(0.6, L, 1e-2, M) u 1:2 w lp ls 1 title 'ED'

		}
	}
}
