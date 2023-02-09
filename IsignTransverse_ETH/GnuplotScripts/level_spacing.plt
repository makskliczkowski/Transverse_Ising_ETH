reset 
##--PREAMBLE
set autoscale
use_png = 0		# 1 if use png output, and 0 for qt output
if(use_png) { set term pngcairo transparent size 1200, 1200 font sprintf("Helvetica,%d",20); }
else {set term qt size 900, 900 font sprintf("Helvetica,%d",20); }
load './gnuplot-colorbrewer-master/diverging/RdYlGn.plt'

set mxtics
set mytics
set style line 12 lc rgb '#ddccdd' lt 1 lw 1.5
set style line 13 lc rgb '#ddccdd' lt 1 lw 0.5
set grid xtics mxtics ytics mytics back ls 12, ls 13
set size square
set xtics nomirror
set key inside bottom right font ",16" spacing 2

fileexist(name)=1#int(system("if exist \"".name."\" ( echo 1) else (echo 0)"))

# Margins for each row resp. column
RANGE = "set xrange[0:1]; set yrange[0:2.0]"
UNSET = "unset tics; unset xlabel; unset ylabel; unset title; unset key; unset border;"
#-- PARAMETERS
model = 0       # 1=symmetries and 0=disorder
w = 0.1
g = 0.9
L = 14
h = 0.8
J=1.0
scaling = 2		# 0 - h scaling / 1 - L scaling / 2 - g scaling / 3 - J scaling / 4 - w scaling (k -  scaling for model==1)
function = 0    # 1 - gap ratio / 0 - prob distribution
h_vs_g = 1      # 1 - as function of h / 0 - as function of g
heatmap = 0
interpolate = 1

if(!heatmap){
    if(scaling == 0) { h_vs_g = 0; }
    if(scaling == 2) { h_vs_g = 1; }
}
	h0 = 5;     hend = 55;		dh = 10;
	g0 = 10;    gend = 90;		dg = 20;
    J0 = 10;    Jend = 60;     dJ = 10
	L0 = 9;	    Lend = 13; 		dL = 1;

    h_list = '0.20 0.60 1.20 1.40 1.60 1.80 2.40 3.00 3.60'
    g_list = '0.20 0.30 0.70 0.80 1.10 1.40'; g00=0.20
	w_num = 5;	array w_list[w_num];
	w_list[1] = 0.01;	w_list[2] = 0.05;	w_list[3] = 0.1;	w_list[4] = 0.3;	w_list[5] = 0.5;

GOE(x) = 0.5307;
Lap(x) = 0.3863;
ADD = function?  " GOE(x) w l dt (8,8) lc rgb 'black' lw 2 notitle, Lap(x) w l dt (8,8) lc rgb 'blue' lw 2 notitle" : " 2/(1+x)**2 w l dt (3,5,10,5) lw 2.5 lc rgb 'black' notitle, 27./4 * (x+x**2)/(1+x+x**2)**2.5 w l dt (3,5,10,5) lw 2.5 lc rgb 'blue' notitle"
dir_base = '../results/ISING/'.(model? 'symmetries' : 'disorder').'/PBC/LevelSpacing/'.(function? 'ratio/' : 'distribution/');
set key top right
load './gnuplot-colorbrewer-master/diverging/RdYlGn.plt'

name_prefix = dir_base.(function? (h_vs_g? "hMap" : "gMap") : "");
_name_ratio(Jx, Lx, dis) = name_prefix.sprintf("_L=%d",Lx).(model? ",k=0,p=1,x=1.dat" : sprintf(",J=%.2f,J0=0.00,g0=0.00,w=%.2f.dat", Jx, dis));
_name_dist(Jx, Lx, hx, gx) = model? name_prefix.sprintf("_L=%d,g=%.2f,h=%.2f,k=0,p=1,x=1.dat", Lx, gx, hx) :\
                        name_prefix.sprintf("_L=%d,J=%.2f,J0=0.00,g=%.2f,g0=0.00,h=%.2f,w=%.2f.dat", Lx, Jx, gx, hx, w);
set output "LevelSpacing/level_spacing_J=1.0.png";
if(function){
    if(heatmap){
        set xrange[0.05:1.5]
        set yrange[0.05:1.5]
        set cbrange[0.38:0.53]
        set ylabel (h_vs_g? "h" : "g")
        set xlabel (h_vs_g? "g" : "h")
        #set cblabel "<r>"
        if(interpolate){
            set pm3d map
            set dgrid3d 30, 30, 1
    	    set pm3d interpolate 0,0
            splot _name_ratio(J, L, w) u 1:2:3 with pm3d
        } else{
        	set view map
    	    set pm3d interpolate 0,0
            splot _name_ratio(J, L, w) u 1:2:3 with image
        }
    } else {
        set key right center
        set ylabel '<r>'
        set xlabel (h_vs_g? "h" : "g")
        set yrange [0.3:0.55]
        #plot for[hx in h_list] _name(1. * hx) u 1:2 w lp ls ((1.*hx-0.01*h0)/(0.01*dh)) pt 6 ps 1.5 title sprintf("h=%.2f", 1. * hx), @ADD
        #plot for[gx=g0:gend:dg] _name_ratio(L) u 2:($1 == 0.01*gx? $3 : NaN) w lp ls ((gx-g0)/dg) pt 6 ps 1.5 title sprintf("g=%.2f", 0.01*gx), @ADD
        remove_zeros(x) = abs(x) < 0.1? NaN : x;
        if(scaling == 0){ plot for[hx=h0:hend:dh] _name_ratio(J, L, w) u 2:($1 == 0.01*hx? remove_zeros($3) : NaN) w lp ls ((hx-h0)/dh) pt 6 ps 1.5 title sprintf("h=%.2f", 0.01*hx), @ADD
        } else {
        if(scaling == 1){ plot for[Lx=L0:Lend:dL] _name_ratio(J, Lx, w) u ($1 == (h_vs_g? g : h)? $2 : NaN):(remove_zeros($3)) w lp ls ((Lx-L0)/dL+1) pt 6 ps 1.5 title sprintf("L=%d", Lx), @ADD
        } else { 
        if(scaling == 2){ plot for[gx=g0:gend:dg] _name_ratio(J, L, w) u 2:($1 == 0.01*gx? remove_zeros($3) : NaN) w lp ls ((gx-g0)/dg) pt 6 ps 1.5 title sprintf("g=%.2f", 0.01 * gx), @ADD
        } else {
        if(scaling == 3){ plot for[Jx=J0:Jend:dJ] _name_ratio(0.01*Jx, L, w) u ($1 == (h_vs_g? g : h)? $2 : NaN):(remove_zeros($3)) w lp ls ((Jx-J0)/dJ) pt 6 ps 1.5 title sprintf("J=%.2f", 0.01*Jx), @ADD
        } else {
        if(scaling == 4){ plot for[i=1:w_num] _name_ratio(J, L, w_list[i]) u ($1 == (h_vs_g? g : h)? $2 : NaN):(remove_zeros($3)) w lp ls i pt 6 ps 1.5 title sprintf("w=%.2f", w_list[i]), @ADD
        } else { 

        }}}}}}
} else{
    #plot for[hx in h_list] _name(1. * hx) u 1:2 w lp ls ((1.*hx-0.01*h0)/(0.01*dh)) pt 6 ps 1.5 title sprintf("h=%.2f", 1. * hx), @ADD
    #plot for[gx in g_list] _name_dist(L, h, 1. * gx) u 1:2 w lp ls ((1. * gx - 0.01*g0)/ (0.01*dg)) pt 6 ps 1.5 title sprintf("g=%.2f", 1. * gx), @ADD
    set ylabel 'P(s)' rotate by 0;
    set xlabel 's'
    set label 1 at 0.2, 1.5 "Poisson" tc rgb "black"
    set label 2 at 0.7, 1.2 "Wigner-Dyson" tc rgb "blue"
    if(scaling == 0){ plot for[hx=h0:hend:dh] _name_dist(J, L, 0.01 * hx, g) u 1:2 w lp ls ((hx-h0)/dh) pt 6 ps 1.5 title sprintf("h=%.2f", 0.01*hx), @ADD
    } else {
    if(scaling == 1){ plot for[Lx=L0:Lend:dL] _name_dist(J, Lx, h, g) u 1:2 w lp ls ((Lx-L0)/dL) pt 6 ps 1.5 title sprintf("L=%d", Lx), @ADD
    } else {
    if(scaling == 2){ plot for[gx=g0:gend:dg] _name_dist(J, L, h, 0.01*gx) u 1:2 w lp ls ((gx-g00)/dg) pt 6 ps 1.5 title sprintf("g=%.2f", 0.01*gx), @ADD
    } else {
    if(scaling == 3){ plot for[Jx=J0:Jend:dJ] _name_dist(0.01*Jx, L, h, g) u 1:2 w lp ls ((Jx-J0)/dJ) pt 6 ps 1.5 title sprintf("J=%.2f", 0.01*Jx), @ADD
    } else{

    }}}}    
}