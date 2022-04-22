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
set key inside bottom right font ",16" spacing 2

#set style line 1 dt (3,5,10,5) lc rgb "black" lw 1.5
#set style line 2 dt (3,3) lc rgb "red" lw 1.5
#set style line 3 dt (8,8) lc rgb "blue" lw 1.5
#set style line 4 dt (1,1) lc rgb "green" lw 1.5

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
w = 0.1
g = 0.05
L = 14
h = 0.1
J=0.1
scaling = 0		     # 0 - h scaling / 1 - L scaling / 2 - g scaling / 3 - J scaling
smoothed = 1         # smoothed ?
plot_der_GOE = 0     # plot deriviation from GOE value

	h0 = 10;     hend = 40;		dh = 5;
	g0 = 5;    gend = 50;		dg = 5;
    J0 = 10;    Jend = 100;     dJ = 10
	L0 = 8;	    Lend = 14; 		dL = 1;

    h_list = '0.20 0.60 1.20 1.40 1.60 1.80 2.40 3.00 3.60'
    g_list = '0.20 0.30 0.70 0.80 1.10 1.40';

GOE(x) = (x < 1? 2 * x - x*log(1+2*x) : 2-x*log( (2*x+1) / (2*x-1)))
eps = 5e-2
ADD=plot_der_GOE? sprintf("%f w l ls 1 dt (3,5,10,5) lc rgb 'black' lw 2 notitle", eps)\
         : "GOE(x) w l ls 1 dt (3,5,10,5) lc rgb 'black' lw 2 t 'GOE'"
dir_base = '../results/'.(model? 'symmetries' : 'disorder').'/PBC/SpectralFormFactor/'
if(smoothed){ dir_base = dir_base.'smoothed/';}

set key top right
load './gnuplot-colorbrewer-master/diverging/RdYlGn.plt'

_name(Lx, Jx, hx, gx) = model? dir_base.sprintf("_L=%d,g=%.2f,h=%.2f,k=0,p=1,x=1.dat", Lx, Jx, gx, hx) :\
                        dir_base.sprintf("_L=%d,J=%.2f,J0=0.00,g=%.2f,g0=0.00,h=%.2f,w=%.2f.dat", Lx, Jx, gx, hx, w);

set logscale xy;
if(plot_der_GOE){ set xrange[1e-3:5]; }
else { set xrange[*:5]; }
set ylabel 'K({/Symbol t})' rotate by 0 offset 2,0;
set xlabel '{/Symbol t}'
data(x, y) = plot_der_GOE? abs(log10( y / GOE(x) )) : y
if(scaling == 0){ plot for[hx=h0:hend:dh] _name(L, J, 0.01 * hx, g) u 1:(data($1, $2)) w l ls ((hx-h0)/dh+1) lw 2 title sprintf("h=%.2f", 0.01*hx), @ADD
} else {
if(scaling == 1){ plot for[Lx=L0:Lend:dL] _name(1.00 * Lx, J, h, g) u 1:(data($1, $2)) w l ls ((Lx-L0)/dL+1) lw 2 title sprintf("L=%d", Lx), @ADD
} else {
if(scaling == 2){ plot for[gx=g0:gend:dg] _name(L, J, h, 0.01 * gx) u 1:(data($1, $2)) w l ls ((gx-g0)/dg+1) lw 2 title sprintf("g=%.2f", 0.01*gx), @ADD
} else {
if(scaling == 3){ plot for[Jx=J0:Jend:dJ] _name(L, 0.01 * Jx, h, g) u 1:(data($1, $2)) w l ls ((Jx-J0)/dJ+1) lw 2 title sprintf("J=%.2f", 0.01*Jx), @ADD
} else{

}}}}  