reset 
##--PREAMBLE
set autoscale
use_png = 0		# 1 if use png output, and 0 for qt output
if(use_png) { set term pngcairo size 1400, 1200 font sprintf("Helvetica,%d",20); }
else {set term qt size 1100, 900 font sprintf("Helvetica,%d",20); }
set mxtics
set mytics
set style line 12 lc rgb '#ddccdd' lt 1 lw 1.5
set style line 13 lc rgb '#ddccdd' lt 1 lw 0.5
set grid xtics mxtics ytics mytics back ls 12, ls 13
set size square
set xtics nomirror
set key outside bottom right font ",20" spacing 1.5
fileexist(name)=system("[ -f '".name."' ] && echo '1' || echo '0'") + 0#int(system("if exist \"".name."\" ( echo 1) else (echo 0)"))

# Margins for each row resp. column
UNSET = "unset tics; unset xlabel; unset ylabel; unset title; unset key; unset border;"

#---------------------------- PARAMETERS
model = 0       # 1=symmetries and 0=disorder
heisenberg=1

w = 0.7
g = 0.9
L = 12
h = 0.0
J = 1.0
Jdis = 0.0
gdis = 0.0
k=1

scaling = 2     # 0 - h scaling / 1 - L scaling / 2 - g scaling / 3 - J scaling / 4 - k scaling (only model=1) : w scaling (only model=0)
smoothed = 0        # smoothed ?
plot_der_GOE = 0	 # plot deriviation from GOE value
zoom_in = 0          # zoom in to collapse on GOE
find_Thouless = 1    # find thouless time?
add_gap_ratio = 1	 # add gap ratio
rescale_times = 0	 # rescale times by size or parameter
nu = -0.5
use_folded = 0		 # use folded spectrum
real_units = 0		 # rescale by Heisenberg time?
compare_folded_to_unfolded = 0



if(scaling < 0 || scaling > 4 || zoom_in == 1) add_gap_ratio = 0;
if(plot_der_GOE){ zoom_in = 0;}

	h0 = 5;     hend = 150;		dh = 5;
	g0 = 10;    gend = 100;		dg = 10;
    J0 = 2;    Jend = 100;     dJ = 2
	L0 = 10;	    Lend = 18; 		dL = 1;
	w0 = 10;		wend = 100;			dw = 1;
	w_num = (wend-w0)/dw + 1;	array w_list[w_num];
	do for[i=1:w_num]{ w_list[i] = 0.01*w0 + 0.01*dw*(i-1);}
    h_list = '0.20 0.60 1.20 1.40 1.60 1.80 2.40 3.00 3.60'
    g_list = '0.20 0.30 0.70 0.80 1.10 1.40';

GOE(x) = (x < 1? 2 * x - x*log(1+2*x) : 2-x*log( (2*x+1) / (2*x-1)))

eps = 8e-2
ADD=plot_der_GOE? sprintf("%f w l ls 1 dt (3,5,10,5) lc rgb 'black' lw 2 notitle", eps)\
         : "GOE(x) w l ls 1 dt (3,5,10,5) lc rgb 'black' lw 2 t 'GOE', (x < 0.2? NaN : 1.0) w l ls 1 dt (3,5,10,5) lc rgb 'black' lw 2 notitle"		 
if(real_units){ ADD = ""; };
dir_base = '../results/'.(heisenberg? 'HEISENBERG/' : 'ISING/').(model? 'symmetries' : 'disorder').'/PBC/SpectralFormFactor/'
if(smoothed){ dir_base = dir_base.'smoothed/';}

#---------------------------- SET PLOT DATA
LOG = (zoom_in? "set logscale x; set format x '10^{%L}'" : "set logscale xy; set format x '10^{%L}'")."; set format y '10^{%L}';"; 
@LOG;
LINE = scaling == 4 && model == 0? "unset logscale y; set format x '10^{%L}'; set format y '%g';" : "unset logscale xy; set format x '%g'; set format y '%g';"

load './gnuplot-colorbrewer-master/diverging/RdYlGn.plt'
out_dir = 'SpectralFormFactor/'
_name_long(Lx, Jx, hx, gx, dis) = dir_base.(use_folded? "folded" : "").(\
								model? sprintf("_L=%d,J=%.2f,g=%.2f,h=%.2f,k=%d,p=1,x=1.dat",Lx, Jx, gx, hx, k) :\
                        			   sprintf("_L=%d,J=%.2f,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f.dat", Lx, Jx, Jdis, gx, gdis, hx, dis));

wvalue(Lx) = Lx == 10? 1.45 : (Lx == 12? 1.55 : (Lx == 14? 1.70 : 1.90))

_name(x) = 0; _key_title(x) = 0;
_rescale_times(x, i) = 0;
i0 = 0; iend = 0; di = 1;
		if(scaling == 0){
			_name(x) = _name_long(L, J, 0.01 * x, g, w);    _key_title(x) = sprintf("h=%.2f", x / 100.)
			if(Jdis == 0){ set key title "w=h/2";}
			i0 = h0; iend = hend; di = dh; 	out_dir = out_dir."h_scaling/";
			output_name = sprintf("_L=%d,J0%.2f,J0=%.2f,g=%.2f,g0=%.2f,w=%.2f", L, J, Jdis, g, gdis, w);
			_rescale_times(x, i) = x * exp(-1. / (0.01*i)**nu)
		}else{
			if(scaling == 1){

				_name(x) = _name_long(x, J, h, g, w);	_key_title(x) = sprintf("L=%d",x);
				i0 = L0; iend = Lend; di = dL; 	out_dir = out_dir."size_scaling/"
				output_name = sprintf("_J=%.2f,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f", J, Jdis, g, gdis, h, w);
				_rescale_times(x, i) = x*exp(0.3*i)
			} else{
				if(scaling == 2){
					_name(x) = _name_long(L, J, h, 0.01 * x, w);    _key_title(x) = sprintf("g=%.2f",0.01*x)
					i0 = g0; iend = gend; di = dg; 	out_dir = out_dir."g_scaling/"
					output_name = sprintf("_L=%d,J=%.2f,J0=%.2f,g0=%.2f,h=%.2f,w=%.2f", L, J, Jdis, gdis, h, w);
					_rescale_times(x, i) = x * exp(10*(0.01*i)**nu)
				} else{
					if(scaling == 3){
						_name(x) = _name_long(L, 0.01 * x, h, g, w);    _key_title(x) = sprintf("J=%.2f", x / 100.)
						i0 = J0; iend = Jend; di = dJ; out_dir = out_dir."J_scaling/"
						output_name = sprintf("_L=%d,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f", L, Jdis, g, gdis, h, w);
						_rescale_times(x, i) = x * exp(-2*(0.01*i)**nu)
					} else{
						if(scaling == 4){
							if(model == 1){
								_name(x) = dir_base.sprintf("_L=%d,J=%.2f,g=%.2f,h=%.2f,k=%d,p=1,x=1.dat", L, J, g, h, x);   
								 _key_title(x) = sprintf("k=%d", x)
								i0 = 1; iend = L; di = 1;
								_rescale_times(x, i) = x;
							} else {
								_name(x) = _name_long(L, J, h, g, w_list[x]);
								_key_title(x) = sprintf("w=%.2f",w_list[x])
								i0 = 1; iend = w_num; di = 1
								_rescale_times(x, i) = x * exp(-1. / (w_list[(i-i0)/di+1])*1)
							}
						} else{
							_name(x) = _name_long(L, J, h, g);    _key_title(x) = sprintf("L=%d,J=%.2f,g=%.2f,h=%.2f", L, J, g, h)
							i0 = 0; iend = 0; di = 1; out_dir = out_dir."single_plots/"
							output_name = sprintf("_L=%d,J=%.2f,J0=%.2f,g=%.2f,g0=%.2f,h=%.2f,w=%.2f", L, J, Jdis, g, gdis, h, w);
					}}}}}
if(!rescale_times){ _rescale_times(x,i)=x;}
#---------------------------- EXTRACT DATA - STATS
    size = (iend - i0) / di+1
    array tau[size]; array y_vals[size];	array tH[size]
    array gap_ratio[size]; array x_vals_gap[size]
    x_min = 1e6;    y_min = 1e6
    do for[i=i0:iend:di]{
		idx = (i-i0)/di+1
		name = _name(i)
		tau[idx] = 0.0; y_vals[idx] = NaN;	gap_ratio[idx] = NaN;
		#print name;
		if(fileexist(name)){
            f(x,y) = x > 2.5? NaN : ((log10( y / GOE(x) )) - eps)**2
			stats name u 1 nooutput; bool = STATS_min;
			if(bool > 0){
				stats name nooutput; n_cols = STATS_columns;
				if(find_Thouless){
					stats name using (f($1, $2)) nooutput prefix "Y";       y_min_idx = Y_index_min;
            		stats name using 1 every ::y_min_idx::(y_min_idx+1) nooutput;   tau[idx] = STATS_min; 
            		stats name every ::y_min_idx::(y_min_idx+1) using 2 nooutput;   y_vals[idx] = STATS_min;
            		stats name using 1 nooutput;   new_min = STATS_min;     if(new_min < x_min){ x_min = new_min; }
            		stats name using 2 nooutput;   new_min = STATS_min;     if(new_min < y_min){ y_min = new_min; }
				} else {
					x_min = 1e-4;
					y_min = 1e-2;
				}
				if(add_gap_ratio && n_cols >= 5){
					stats name every ::0::0 using 5 nooutput; gap_ratio[idx] = STATS_min;
					x_vals_gap[idx] = scaling != 1? 0.01 * i : i;
					if(scaling == 4 && model == 0){ x_vals_gap[idx] = w_list[i]; }
				}
				if(n_cols >= 3){
					stats name every ::0::0 using 3 nooutput; tH[idx] = STATS_min;
				}
			}
		}
        print _key_title(i),"  ", tau[idx]
    }

#---------------------------- GRAPH VISUALS
x_pos = real_units? tau[1] * tH[1] : 1.0
if(!plot_der_GOE){ set arrow from x_pos, graph 0 to x_pos,1 nohead ls 1 dt (3,5,10,5) lc rgb 'black' lw 2;}
if(zoom_in) { unset logscale y; set format y '%g'; set key bottom right font ",20";}
x_max = 98
if(real_units){ x_min = 1e3*x_min;	x_max = x_max * 1e3;	}
print x_min, y_min
#y_min=1e-1
RANGE=zoom_in? sprintf("set xrange[%.6f:20]; set yrange[%.10f:1.2];", x_min, 0.8 * y_min)\
                    : sprintf("set xrange[%.6f:%.6f]; set yrange[%.10f:%.2f];", x_min, x_max, 0.8 * y_min, 0.5*(scaling == 1? 2**Lend : 2**L))
MARGIN = "set lmargin at screen 0.10; set rmargin at screen 0.95; set bmargin at screen 0.10; set tmargin at screen 0.99;"
MARGIN_inset = "set lmargin at screen 0.52; set rmargin at screen 0.92; set bmargin at screen 0.62; set tmargin at screen 0.97;"
#---------------------------- PLOT
set ylabel (plot_der_GOE? '{/Symbol D}K({/Symbol t})' : 'K({/Symbol t})') rotate by 0 offset 2,0;
set xlabel (real_units? '{/Symbol t}*t_H' : '{/Symbol t}')
data(x, y) = plot_der_GOE? abs(log10( y / GOE(x) )) : y
if(plot_der_GOE){   plot for[i=i0:iend:di] _name(i) u 1:(data($1, $2)) w l ls ((i-i0)/di+1) lw 2 title _key_title(i), @ADD
} else {
	_rescale_times_final(x, i) = real_units? _rescale_times(x, i) * tH[ (i-i0)/di+1 ] : _rescale_times(x, i)
    set multiplot
    @RANGE; @MARGIN; plot for[i=i0:iend:di] _name(i) u (_rescale_times_final($1,i)):(data($1, $2)) w l ls ((i-i0)/di+1) lw 2 title _key_title(i), @ADD
	if(add_gap_ratio){
		@LINE; @MARGIN_inset; unset xlabel; unset ylabel; unset title; unset key; unset arrow;
		x_min = scaling != 1? 0.009 * i0 : i0-1;		x_max = scaling != 1? 0.011 * iend : iend+1;
		if(scaling == 4){
			x_min = model == 1? 0 : 0.9 * w_list[1];	x_max = model == 1? L+1 : 1.1 * w_list[w_num];
		}
		GOE(x) = 0.5307;
		Lap(x) = 0.3863;
		set yrange[0.361:0.549];
		set xrange[x_min:x_max];
		set label 1 at 1.5*x_min, 0.54 'GOE' front
		set label 2 at (scaling == 4? 0.25 : 0.75)*x_max, 0.37 'Poisson' front 
		plot gap_ratio using (x_vals_gap[$1]):(gap_ratio[$1]) w p pt 7 ps 2.0 lw 2 lc rgb 'black' notitle,\
			 GOE(x) w l dt (8,8) lc rgb 'black' lw 2 notitle, Lap(x) w l dt (8,8) lc rgb 'black' lw 2 notitle
	}
	if(find_Thouless){
		unset label;
		@LOG; @RANGE; @MARGIN; @UNSET; plot for[i=1:size] '+' using (_rescale_times_final(tau[i],(i-1)*di+i0)):(y_vals[i]) w lp ls i pt 7 ps 2.0 lw 2 notitle
		#@RANGE; @MARGIN; @UNSET; plot tau using (tau[$1]):(y_vals[$1]) w p pt 6 ps 2.0 lw 2 lc rgb 'black' notitle
	}
	if(compare_folded_to_unfolded){
		@RANGE; @MARGIN; plot for[i=i0:iend:di] (use_folded? "" : "folded")._name(i) u (_rescale_times_final($1,i)):(data($1, $2)) w l ls ((i-i0)/di+1) lw 2 dt (8,8) title _key_title(i)
	}
    unset multiplot
}
