reset 
##--PREAMBLE
set autoscale
use_png = 0		# 1 if use png output, and 0 for qt output
if(use_png) { set term pngcairo size 1200, 1200 font sprintf("Helvetica,%d",18); }
else {set term qt size 900, 900 font sprintf("Helvetica,%d",16); }
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

fileexist(name)=1#int(system("if exist \"".name."\" ( echo 1) else (echo 0)"))

# Margins for each row resp. column
TMARGIN = "set tmargin at screen 0.99; set bmargin at screen 0.56"
BMARGIN = "set tmargin at screen 0.46; set bmargin at screen 0.10"
LMARGIN = "set lmargin at screen 0.10; set rmargin at screen 0.46"
RMARGIN = "set lmargin at screen 0.46; set rmargin at screen 0.82"
RANGE = "set xrange[0:1]; set yrange[0:2.0]"
UNSET = "unset tics; unset xlabel; unset ylabel; unset title; unset key; unset border;"
#-- PARAMETERS
model = 1       # 1=symmetries and 0=disorder
w = 0.0
g = 0.2
L = 19
h = 1.0
scaling = 2		# 0 - h scaling / 1 - L scaling / 2 - g scaling
function = 0    # 1 - gap ratio / 0 - prob distribution
h_vs_g = 1      # as function of h
if(scaling == 0) h_vs_g = 0;
if(scaling == 2) h_vs_g = 1;

	h0 = 20;	hend = 120;		dh = 20;
	g0 = 10;	gend = 50;		dg = 5;
	L0 = 14;	Lend = 19; 		dL = 1;

GOE(x) = x < 1.5? 0.5307 : NaN;
Lap(x) = x > 2.? 0.3863 : NaN;
ADD = function?  " GOE(x) w l ls 3 lw 2 notitle, Lap(x) w l ls 3 lw 2 notitle" : " 2/(1+x)**2 w l ls 2 notitle, 27./4 * (x+x**2)/(1+x+x**2)**2.5 w l ls 2 notitle"
dir_base = '../results/'.(model? 'symmetries' : 'disorder').'/PBC/LevelSpacing/'.(function? 'ratio/' : 'distribution/');
set key top right
load './gnuplot-colorbrewer-master/diverging/RdYlGn.plt'

_name_suffix(Lx, hx, gx) = 0;
if(function){
   _name_suffix(Lx, hx, gx) = model? sprintf("_L=%d",Lx).(h_vs_g? sprintf(",g=%.2f", gx) : sprintf(",h=%.2f", hx)).",k=0,p=1,x=1.dat" :\
                        sprintf("_L=%d,J0=0.00",Lx).(h_vs_g? sprintf(",g=%.2f,g0=0.00", gx) : sprintf(",g0=0.00,h=%.2f", hx)).",w=%.2f.dat";
} else {
   _name_suffix(Lx, hx, gx) = model? sprintf("_L=%d,g=%.2f,h=%.2f,k=0,p=1,x=1.dat", Lx, gx, hx) :\
                        sprintf("_L=%d,J0=0.00,g=%.2f,g0=0.00,h=%.2f,w=%.2f.dat", Lx, gx, hx, w);
}
if(scaling == 1){
    _name(x) = dir_base._name_suffix(x, h, g);
    plot for[Lx=L0:Lend:dL] _name(Lx) u 1:2 w lp ls ((Lx-L0)/dL) pt 6 ps 1.5 title sprintf("L=%d", Lx), @ADD
} else {
    if(scaling == 0){
        _name(x) = dir_base._name_suffix(L, x, g);
        h_list = '0.20 0.60 1.20 1.40 1.60 1.80 2.40 3.00 3.60'
        #plot for[hx in h_list] _name(1. * hx) u 1:2 w lp ls ((1.*hx-0.01*h0)/(0.01*dh)) pt 6 ps 1.5 title sprintf("h=%.2f", 1. * hx), @ADD
        plot for[hx=h0:hend:dh] _name(0.01 * hx) u 1:2 w lp ls ((hx-h0)/dh) pt 6 ps 1.5 title sprintf("h=%.2f", 0.01*hx), @ADD
    } else {
        _name(x) = dir_base._name_suffix(L, h, x);
        g_list = '0.20 0.50 0.70 0.80 1.40 2.00'
        plot for[gx in g_list] _name(1. * gx) u 1:2 w lp ls ((1. * gx - 0.01*g0)/ (0.01*dg)) pt 6 ps 1.5 title sprintf("g=%.2f", 1. * gx), @ADD
    }
}