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
load './gnuplot-colorbrewer-master/diverging/RdYlGn.plt'
set fit quiet

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
MARGIN = "set lmargin at screen 0.10; set rmargin at screen 0.99; set bmargin at screen 0.10; set tmargin at screen 0.99;"
UNSET = "unset tics; unset xlabel; unset ylabel; unset title; unset border;"
#---------------------------- PARAMETERS
model = 0       # 1=symmetries and 0=disorder
w = 0.3
g = 0.9
L = 14
h = 0.8
J = 0.1

what = 0                # 0 - DOS / 1 - Lvl-spacing Distribution / 2 - unfolding analysis
plot_typical_wH = 0     # plot only the scaling of wH_typical
scaling = 1             # 0 - J scaling / 1 - w scaling

use_unfolded = 0        # use unfolded energies for DOS
use_logarithmic = 1     # use low(w) distribution
compare_fits = 1        # compare different ploynomial degrees in unfolding
use_fit = 0             # use fit for distribtution
    J0 = 5;    Jend = 100;  dJ = 5
    w0 = 10;    wend = 60;  dw = 10

if(plot_typical_wH){
    J0 = 5; Jend = 100; dJ = 5;
    w0 = 10;    wend = 300;  dw = 10
}
if(scaling){     i0 = w0;    iend = wend;    di = dw;
} else {         i0 = J0;    iend = Jend;    di = dJ; }

_name_long(Lx, Jx, hx, gx, dis) = (model? sprintf("_L=%d,J=%.2f,g=%.2f,h=%.2f.dat", Lx, Jx, gx, hx) :\
                        			   sprintf("_L=%d,J=%.2f,J0=0.00,g=%.2f,g0=0.00,h=%.2f,w=%.2f.dat", Lx, Jx, gx, hx, dis));
_name(x) = scaling? _name_long(L, J, h, g, x) : _name_long(L, x, h, g, w);
_key_title(x) = scaling? sprintf("w=%.2f", x) : sprintf("J=%.2f", x)
#--------- DIRECTORIES
dir_base = '../results/'.(model? 'symmetries' : 'disorder').'/PBC/'
dir_DOS = dir_base.'DensityOfStates/'.(use_unfolded? 'unfolded' : '')
dir_dist = dir_base.'LevelSpacingDistribution/unfolded'.(use_logarithmic? '_log' : '');
dir_unfolding = dir_base.'Unfolding/';

#--------- FIT FUNCITONS for x=log10(w):
fun_fit(x) = a / 10**x * exp(- (x*log(10) - mu)**2 / (2*sig**2))
plot_fit(x, a, mu, sig) = a / 10**x * exp(-(x*log(10) - mu)**2 / (2*sig**2))

E0 = -sqrt(g**2 + (h+w)**2);

#--------- LOOP THROUGH CURVES AND ANALYZE
sizee = (iend-i0)/di + 1;
array wH_typical[sizee];    array vals[sizee];
array wH[sizee];
array a_list[sizee];    array mu_list[sizee];    array sig_list[sizee];
do for[i=i0:iend:di]{
    idx = (i-i0)/di + 1
    vals[idx] = 0.01*i;
    _name_ = dir_dist._name(0.01*i)
    sig = 0.5; mu = -1; a = 1.0
    if(fileexist(_name_)){
        stats _name_ using 3 every ::0::0 nooutput;   wH_typical[idx] = STATS_min;
        stats _name_ using 4 every ::0::0 nooutput;   wH[idx] = STATS_min;
        if(use_fit){ fit fun_fit(x) dir_dist._name(0.01*i) u 1:2 via a, mu, sig; }
        a_list[idx] = a;    mu_list[idx] = mu;  sig_list[idx] = sig;
        print wH_typical[idx], wH[idx], a_list[idx], mu_list[idx], sig_list[idx]
    }
}
array E_noninteracting[L+1];
if(what == 0){
    do for[i=1:L+1]{ E_noninteracting[i] = (L - 2*(i-1)) * E0; }
}

#--------- RMT predictions
PI = 3.141592653
P_GOE(x) = PI / 2.0 * x * exp(-PI * x**2. / 4.0)# * 1.55/1.9;
P_POISSON(x) = exp(-x)

#--------- ANALYZE
if(plot_typical_wH){
    set ylabel 'log {/Symbol w}_H^{typ}' rotate by 0
    set xlabel (scaling? 'w' : 'J')
    plot wH_typical u (vals[$1]):(wH_typical[$1]) w lp pt 7 notitle
    exit;
}

if(what == 0){
    SCALE = use_unfolded? sprintf("set xrange[0:%d]; set yrange[0.0:1e-3];", 2**L) : "set xrange[-15:15]; set yrange[0.0:0.3];"
    set ylabel 'DOS';
    set xlabel '{/Symbol w}';
    set multiplot
    @MARGIN; @SCALE;   plot E_noninteracting u (E_noninteracting[$1]):(0.0) w p ps 3 pt 7 lc rgb "black" notitle
    @UNSET; 
    @MARGIN; @SCALE;   plot for[i=i0:iend:di] dir_DOS._name(0.01*i) w steps lw 2 t _key_title(0.01*i);
    unset multiplot
} else {
    if(what == 1){
        #set logscale y; set yrange[1e-3:2]
        set key inside right top
        SCALE = use_logarithmic? "set xrange[-3:2]; set yrange[1e-3:1.7];" : "set xrange[1e-3:10]; set yrange[1e-3:2.7];" 
        if(use_logarithmic){ set ylabel 'P(log_{10} {/Symbol w})';    set xlabel 'log_{10} {/Symbol w}'}
        else {set ylabel 'P({/Symbol w})';    set xlabel '{/Symbol w}'; };#set logscale x;};
        trueX(x) = use_logarithmic? 10**x : x
        factor = use_logarithmic? log(10) : 1.0
        set multiplot
        @MARGIN; @SCALE;    plot for[i=i0:iend:di] dir_dist._name(0.01*i) u 1:2 w steps lw 2 t _key_title(0.01*i), '+' u (NaN):(NaN) w p ps 2 pt 7 lc rgb "black" t 'log {/Symbol w}_H^{typ}',\
                                     factor*trueX(x)*P_GOE(trueX(x)) w l ls 1 dt (2,2) lc rgb 'black' lw 2 t 'GOE',\
                                     factor*trueX(x)*P_POISSON(trueX(x)) w l ls 1 dt (4,3,2, 2) lc rgb 'red' lw 2 t 'Poisson', '+' u (NaN):(NaN) w l dt (8,8) lw 2 lc 1 t 'log-normal fit'
        @UNSET; 
        @MARGIN; @SCALE;    plot for[i=i0:iend:di] '+' u (wH_typical[(i-i0)/di + 1]):(0.0) w p ps 2 pt 7 notitle
        if(use_fit){
            @UNSET; 
            @MARGIN; @SCALE;    plot for[i=1:sizee] abs(plot_fit(x, a_list[i], mu_list[i], sig_list[i])) w l dt (8,8) lw 2 notitle
        }
        unset multiplot
    } else {
        if(compare_fits){
            plot for[i=3:6] dir_unfolding._name(scaling? w : J) u 1:(column(i)) w l lw 2 t sprintf("n=%d", i), dir_unfolding._name(scaling? w : J) u 1:2 w steps lw 3 lc rgb "black" t 'cdf' 
        } else {
            SCALE = sprintf("set xrange[-18:23]; set yrange[-200:%.2f];", 1.1 * 2**L)
            if(use_logarithmic){ set ylabel 'P(log {/Symbol w})';    set xlabel 'log {/Symbol w}'}
            else {set ylabel 'P({/Symbol w})';    set xlabel '{/Symbol w}'; set logscale x;};
            set multiplot
            @MARGIN; @SCALE;    plot for[i=i0:iend:di] dir_unfolding._name(0.01*i) u 1:2 w steps lw 2 t _key_title(0.01*i)
            @UNSET; 
            #@MARGIN; @SCALE;    plot for[i=J0:Jend:dJ] dir_unfolding._name(0.01*i) u 1:6 w l dt (8,8) lw 2 notitle
            unset multiplot
        }
    }
}