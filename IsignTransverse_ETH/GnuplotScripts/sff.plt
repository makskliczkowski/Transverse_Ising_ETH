reset 
##--PREAMBLE
set autoscale
use_png = 0		# 1 if use png output, and 0 for qt output
if(use_png) { set term pngcairo size 1200, 1200 font sprintf("Helvetica,%d",20); }
else {set term qt size 900, 900 font sprintf("Helvetica,%d",20); }
set mxtics
set mytics
set style line 12 lc rgb '#ddccdd' lt 1 lw 1.5
set style line 13 lc rgb '#ddccdd' lt 1 lw 0.5
set grid xtics mxtics ytics mytics back ls 12, ls 13
set size square
set xtics nomirror
set key inside bottom right font ",16" spacing 2 maxrows 6

#set style line 1 dt (3,5,10,5) lc rgb "black" lw 1.5
#set style line 2 dt (3,3) lc rgb "red" lw 1.5
#set style line 3 dt (8,8) lc rgb "blue" lw 1.5
#set style line 4 dt (1,1) lc rgb "green" lw 1.5

fileexist(name)=system("[ -f '".name."' ] && echo '1' || echo '0'") + 0#int(system("if exist \"".name."\" ( echo 1) else (echo 0)"))

# Margins for each row resp. column
UNSET = "unset tics; unset xlabel; unset ylabel; unset title; unset key; unset border;"

#---------------------------- PARAMETERS
model = 0       # 1=symmetries and 0=disorder
w = 0.1
g = 0.2
L = 13
h = 0.8
J = 1.0
J_knot = 0.; g_knot = 0.; 
scaling = 2		     # 0 - h scaling / 1 - L scaling / 2 - g scaling / 3 - J scaling
smoothed = 1         # smoothed ?
plot_der_GOE = 0     # plot deriviation from GOE value
zoom_in = 0          # zoom in to collapse on GOE
find_Thouless = 1    # find thouless time?
	h0 = 10;     hend = 50;		dh = 5;
	g0 = 10;    gend = 100;		dg = 10;
    J0 = 10;    Jend = 100;     dJ = 20
	L0 = 9;	    Lend = 13; 		dL = 1;

    h_list = '0.20 0.60 1.20 1.40 1.60 1.80 2.40 3.00 3.60'
    g_list = '0.20 0.30 0.70 0.80 1.10 1.40';

GOE(x) = (x < 1? 2 * x - x*log(1+2*x) : 2-x*log( (2*x+1) / (2*x-1)))
eps = 2.5e-1
ADD=plot_der_GOE? sprintf("%f w l ls 1 dt (3,5,10,5) lc rgb 'black' lw 2 notitle", eps)\
         : "GOE(x) w l ls 1 dt (3,5,10,5) lc rgb 'black' lw 2 t 'GOE', (x < 0.2? NaN : 1.0) w l ls 1 dt (3,5,10,5) lc rgb 'black' lw 2 notitle"
dir_base = '../results/'.(model? 'symmetries' : 'disorder').'/PBC/SpectralFormFactor/'
if(smoothed){ dir_base = dir_base.'smoothed/';}

#---------------------------- SET PLOT DATA
load './gnuplot-colorbrewer-master/diverging/RdYlGn.plt'
out_dir = 'SpectralFormFactor/'
_name_long(Lx, Jx, hx, gx) = model? dir_base.sprintf("_L=%d,g=%.2f,h=%.2f,k=0,p=1,x=1.dat", Lx, Jx, gx, hx) :\
                        dir_base.sprintf("_L=%d,J=%.2f,J0=0.00,g=%.2f,g0=0.00,h=%.2f,w=%.2f.dat", Lx, Jx, gx, hx, w);

_name(x) = 0; _key_title(x) = 0;
i0 = 0; iend = 0; di = 1;
		if(scaling == 0){
			_name(x) = _name_long(L, J, 0.01 * x, g);    _key_title(x) = sprintf("h=%.2f", x / 100.)
			i0 = h0; iend = hend; di = dh; 	out_dir = out_dir."h_scaling/";
			output_name = sprintf("_L=%d,J0%.2f,J0=%.2f,g=%.2f,g0=%.2f,w=%.2f", L, J, J_knot, g, g_knot, w);
		}else{
			if(scaling == 1){
				_name(x) = _name_long(x, J, h, g);	_key_title(x) = sprintf("L=%d",x);
				i0 = L0; iend = Lend; di = dL; 	out_dir = out_dir."size_scaling/"
				output_name = sprintf("_J=%.2f,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f", J, J_knot, g, g_knot, h, w);
			} else{
				if(scaling == 2){
					_name(x) = _name_long(L, J, h, 0.01 * x);    _key_title(x) = sprintf("g=%.2f",0.01*x)
					i0 = g0; iend = gend; di = dg; 	out_dir = out_dir."g_scaling/"
					output_name = sprintf("_L=%d,J=%.2f,J0=%.2f,g0=%.2f,h=%.2f,w=%.2f", L, J, J_knot, g_knot, h, w);
				} else{
					if(scaling == 3){
						_name(x) = _name_long(L, 0.01 * x, h, g);    _key_title(x) = sprintf("J=%.2f", x / 100.)
						i0 = J0; iend = Jend; di = dJ; out_dir = out_dir."J_scaling/"
						output_name = sprintf("_L=%d,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f", L, J_knot, g, g_knot, h, w);
					} else{
                    }}}}

#---------------------------- EXTRACT DATA - STATS
    size = (iend - i0) / di+1
    array tau[size]; array y_vals[size]
    x_min = 1e6;    y_min = 1e6
	if(find_Thouless){
    do for[i=i0:iend:di]{
		idx = (i-i0)/di+1
		name = _name(i)
		if(fileexist(name)){
            f(x,y) = x > 2.5? NaN : ((log10( y / GOE(x) )) - eps)**2
			stats name nooutput; n_cols = STATS_columns;
			stats name using (f($1, $2)) nooutput prefix "Y";       y_min = Y_index_min;
            stats name using 1 every ::y_min::(y_min+1) nooutput;   tau[idx] = STATS_min; 
            stats name every ::y_min::(y_min+1) using 2 nooutput;   y_vals[idx] = STATS_min;
            stats name using 1 nooutput;   new_min = STATS_min;     if(new_min < x_min){ x_min = new_min; }
            stats name using 2 nooutput;   new_min = STATS_min;     if(new_min < y_min){ y_min = new_min; }
		}
		else{
			tau[idx] = NaN;
		}
        print _key_title(i),"  ", tau[idx]
	}
    }

#---------------------------- GRAPH VISUALS
set logscale xy;
set format x '10^{%L}'; set format y '10^{%L}';
if(plot_der_GOE){ zoom_in = 0; set key bottom left font ",20";}
else { set arrow from 1, graph 0 to 1,1 nohead ls 1 dt (3,5,10,5) lc rgb 'black' lw 2; set key top right font ",20";}
if(zoom_in) { unset logscale y; set format y '%g'; set key bottom right font ",20";}
RANGE=zoom_in? "set xrange[1e-3:7]; set yrange[0:1.5];"\
                    : sprintf("set xrange[%.6f:7]; set yrange[%.5f:%.2f];", x_min, 0.8 * y_min, 0.1*(scaling == 1? 2**Lend : 2**L))
print RANGE
MARGIN = "set lmargin at screen 0.10; set rmargin at screen 0.95; set bmargin at screen 0.10; set tmargin at screen 0.99;"
#---------------------------- PLOT
set ylabel (plot_der_GOE? '{/Symbol D}K({/Symbol t})' : 'K({/Symbol t})') rotate by 0 offset 2,0;
set xlabel '{/Symbol t}'
data(x, y) = plot_der_GOE? abs(log10( y / GOE(x) )) : y
if(plot_der_GOE){   plot for[i=i0:iend:di] _name(i) u 1:(data($1, $2)) w l ls ((i-i0)/di+1) lw 2 title _key_title(i), @ADD
} else {
    set multiplot
    @RANGE; @MARGIN; plot for[i=i0:iend:di] _name(i) u 1:(data($1, $2)) w l ls ((i-i0)/di+1) lw 2 title _key_title(i), @ADD
    @RANGE; @MARGIN; @UNSET; plot for[i=1:size] '+' using (tau[i]):(y_vals[i]) w lp ls i pt 7 ps 2.0 lw 2 notitle
	#@RANGE; @MARGIN; @UNSET; plot tau using (tau[$1]):(y_vals[$1]) w p pt 6 ps 2.0 lw 2 lc rgb 'black' notitle
    unset multiplot
}
