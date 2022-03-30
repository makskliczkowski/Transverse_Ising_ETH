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
model = 0       # 1=symmetries and 0=disorder
w = 0.01
g = 0.6
L = 14
h = 1.0
scaling = 0		# 0 - h scaling / 1 - L scaling / 2 - g scaling
function = 1    # 1 - gap ratio / 0 - prob distribution
h_vs_g = 1      # 1 - as function of h / 0 - as function of g
heatmap = 1
interpolate = 1

if(!heatmap){
    if(scaling == 0) { h_vs_g = 0; }
    if(scaling == 2) { h_vs_g = 1; }
}
	h0 = 20;	hend = 120;		dh = 20;
	g0 = 10;	gend = 100;		dg = 10;
	L0 = 14;	Lend = 19; 		dL = 1;

GOE(x) = x < 1.5? 0.5307 : NaN;
Lap(x) = x > 2.? 0.3863 : NaN;
ADD = function?  " GOE(x) w l ls 3 lw 2 notitle, Lap(x) w l ls 3 lw 2 notitle" : " 2/(1+x)**2 w l ls 2 notitle, 27./4 * (x+x**2)/(1+x+x**2)**2.5 w l ls 2 notitle"
dir_base = '../results/'.(model? 'symmetries' : 'disorder').'/PBC/LevelSpacing/'.(function? 'ratio/' : 'distribution/');
set key top right
load './gnuplot-colorbrewer-master/diverging/RdYlGn.plt'

name_prefix = dir_base.(function? (h_vs_g? "hMap" : "gMap") : "");
_name_ratio(Lx) = name_prefix.sprintf("_L=%d",Lx).(model? ",k=0,p=1,x=1.dat" : sprintf(",J0=0.00,g0=0.00,w=%.2f.dat", w));
_name_dist(Lx, hx, gx) = model? name_prefix.sprintf("_L=%d,g=%.2f,h=%.2f,k=0,p=1,x=1.dat", Lx, gx, hx) :\
                        name_prefix.sprintf("_L=%d,J0=0.00,g=%.2f,g0=0.00,h=%.2f,w=%.2f.dat", Lx, gx, hx, w);
if(function){
    if(heatmap){
        set xrange[0.05:1.2]
        set yrange[0.05:1.2]
        set cbrange[0.39:0.53]
        set ylabel (h_vs_g? "h" : "g")
        set xlabel (h_vs_g? "g" : "h")
        #set cblabel "<r>"
        if(interpolate){
            set pm3d map
            set dgrid3d 24, 24, 1
    	    set pm3d interpolate 0,0
            splot _name_ratio(L) u 1:2:3 with pm3d
        } else{
        	set view map
    	    set pm3d interpolate 0,0
            splot _name_ratio(L) u 1:2:3 with image
        }
    } else {
        if(scaling == 1){
            data(x) = $1 == (h_vs_g? g : h)? x : NaN
            plot for[Lx=L0:Lend:dL] _name_ratio(Lx) u 2:3 w lp ls ((Lx-L0)/dL) pt 6 ps 1.5 title sprintf("L=%d", Lx), @ADD
        } else {
            if(scaling == 0){
                h_list = '0.20 0.60 1.20 1.40 1.60 1.80 2.40 3.00 3.60'
                #plot for[hx in h_list] _name(1. * hx) u 1:2 w lp ls ((1.*hx-0.01*h0)/(0.01*dh)) pt 6 ps 1.5 title sprintf("h=%.2f", 1. * hx), @ADD
                plot for[hx=h0:hend:dh] _name_ratio(L) u 2:($1 == 0.01*hx? $3 : NaN) w lp ls ((hx-h0)/dh) pt 6 ps 1.5 title sprintf("h=%.2f", 0.01*hx), @ADD
            } else {
                g_list = '0.20 0.50 0.70 0.80 1.40 2.00'
                #plot for[gx in g_list] _name u 2:($1 == 1.*gx? $3 : NaN) w lp ls ((1. * gx - 0.01*g0)/ (0.01*dg)) pt 6 ps 1.5 title sprintf("g=%.2f", 1. * gx), @ADD
                plot for[gx=g0:gend:dg] _name_ratio(L) u 2:($1 == 0.01*gx? $3 : NaN) w lp ls ((gx-g0)/dg) pt 6 ps 1.5 title sprintf("g=%.2f", 0.01*gx), @ADD
            }
        }
    }
} else{
    if(scaling == 1){
        plot for[Lx=L0:Lend:dL] _name_dist(Lx, h, g) u 1:2 w lp ls ((Lx-L0)/dL) pt 6 ps 1.5 title sprintf("L=%d", Lx), @ADD
    } else {
        if(scaling == 0){
            h_list = '0.20 0.60 1.20 1.40 1.60 1.80 2.40 3.00 3.60'
            #plot for[hx in h_list] _name(1. * hx) u 1:2 w lp ls ((1.*hx-0.01*h0)/(0.01*dh)) pt 6 ps 1.5 title sprintf("h=%.2f", 1. * hx), @ADD
            plot for[hx=h0:hend:dh] _name_dist(L, 0.01 * hx, g) u 1:2 w lp ls ((hx-h0)/dh) pt 6 ps 1.5 title sprintf("h=%.2f", 0.01*hx), @ADD
        } else {
            g_list = '0.20 0.50 0.70 0.80 1.40 2.00'
            plot for[gx in g_list] _name_dist(L, h, 1. * gx) u 1:2 w lp ls ((1. * gx - 0.01*g0)/ (0.01*dg)) pt 6 ps 1.5 title sprintf("g=%.2f", 1. * gx), @ADD
        }
    }    
}