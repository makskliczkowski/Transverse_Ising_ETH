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
w = 0.2
g = 0.9
L = 16
h = 0.8
J = 0.0

what = 1                    # 0 - DOS / 1 - Lvl-spacing Distribution / 2 - unfolding analysis
plot_typical_wH = 0         # plot only the scaling of wH_typical
scaling = 2                 # 0 - J scaling / 1 - L scaling / 2 - w scaling

use_unfolded = 1            # use unfolded energies for DOS
use_logarithmic = 1         # use low(w) distribution
compare_fits = 1            # compare different ploynomial degrees in unfolding
use_fit = 0                 # use fit for distribtution
use_one_band_picture = 0    # use only middle band

    J0 = 5;    Jend = 100;  dJ = 20
    w0 = 20;    wend = 50;  dw = 10

if(plot_typical_wH){
    J0 = 5; Jend = 100; dJ = 5;
    w0 = 20;    wend = 300;  dw = 10
}

_name_long(Lx, Jx, hx, gx, dis) = (model? sprintf("_L=%d,J=%.2f,g=%.2f,h=%.2f.dat", Lx, Jx, gx, hx) :\
                        			   sprintf("_L=%d,J=%.2f,J0=0.00,g=%.2f,g0=0.00,h=%.2f,w=%.2f.dat", Lx, Jx, gx, hx, dis));

if(scaling == 0){    
    i0 = J0;    iend = Jend;    di = dJ;
    _name(x) = _name_long(L, 0.01*x, h, g, w);  _key_title(x) = sprintf("J=%.2f", 0.01*x)
  } else { if(scaling == 1) {   
        i0 = 10;    iend = 18;  di = 2
        _name(x) = _name_long(x, J, h, g, w);  _key_title(x) = sprintf("L=%d", x)
    } else { if(scaling == 2) {  
        i0 = w0;    iend = wend;    di = dw;
        _name(x) = _name_long(L, J, h, g, 0.01*x);  _key_title(x) = sprintf("w=%.2f", 0.01*x)
    } else {

    }}}
#--------- DIRECTORIES
dir_base = '../results/'.(model? 'symmetries' : 'disorder').'/PBC/'
prefix = use_one_band_picture? "oneband" : ""
dir_DOS = dir_base.'DensityOfStates/'.prefix.(use_unfolded? 'unfolded' : '')
dir_dist = dir_base.'LevelSpacingDistribution/'.prefix.(use_unfolded? 'unfolded' : '').(use_logarithmic? '_log' : '');
dir_unfolding = dir_base.'Unfolding/';

#--------- FIT FUNCITONS for x=log10(w):
fun_fit(x) = a / 10**x * exp(- (x*log(10) - mu)**2 / (2*sig**2))
plot_fit(x, a, mu, sig) = a / 10**x * exp(-(x*log(10) - mu)**2 / (2*sig**2))

E0 = -sqrt(g**2 + (h+w)**2);

#--------- LOOP THROUGH CURVES AND ANALYZE
sizee = (iend-i0)/di + 1;
array wH_typical[sizee];    array wH[sizee];    array vals[sizee];
array omega_min[sizee];  array omega_max[sizee];    array dos_max[sizee];
array a_list[sizee];    array mu_list[sizee];    array sig_list[sizee];
do for[i=i0:iend:di]{
    idx = (i-i0)/di + 1
    vals[idx] = i;
    _name_ = (what == 0? dir_DOS : dir_dist)._name(i)
    sig = 0.5; mu = -1; a = 1.0
    if(fileexist(_name_)){
        if(what == 0){
            stats _name_ using 1 nooutput;   omega_min[idx] = STATS_min;
            stats _name_ using 1 nooutput;   omega_max[idx] = STATS_max;
            stats _name_ using 2 nooutput;   dos_max[idx] = STATS_max;
            print omega_min[idx], omega_max[idx], dos_max[idx]
        } else {
            stats _name_ using 4 every ::0::0 nooutput;   wH[idx] = STATS_min;
            if(use_unfolded){ stats _name_ using 3 every ::0::0 nooutput;   wH_typical[idx] = STATS_min; 
            } else { stats _name_ using 5 every ::0::0 nooutput;   wH_typical[idx] = log(exp(STATS_min)/wH[idx]); }

            if(use_fit){ fit fun_fit(x) dir_dist._name(i) u 1:2 via a, mu, sig; }
            a_list[idx] = a;    mu_list[idx] = mu;  sig_list[idx] = sig;
            print wH_typical[idx], wH[idx], a_list[idx], mu_list[idx], sig_list[idx]
        }
    }
}
array E_noninteracting[L+1];
if(what == 0){
    do for[i=1:L+1]{ E_noninteracting[i] = (L - 2*(i-1)) * E0; }
}

#--------- RMT predictions
PI = 3.141592653
P_GOE(x) = PI / 2.0 * x * exp(-PI * x**2. / 4.0)# * 1.55/1.9;
lambda = 1.0
P_POISSON(x) = lambda * exp(-lambda * x)

#--------- ANALYZE
if(plot_typical_wH){
    set ylabel 'log {/Symbol w}_H^{typ}' rotate by 0
    set xlabel (scaling? 'w' : 'J')
    plot wH_typical u (vals[$1]):(wH_typical[$1]) w lp pt 7 notitle
    exit;
}

if(what == 0){
    SCALE = use_unfolded? sprintf("set xrange[0:1]; set yrange[0.0:1];", 2**L) : "set xrange[-5:5]; set yrange[0.0:1];"
    set ylabel 'normalised DOS';
    set xlabel '{/Symbol w}';
    omega(x, idx) = use_unfolded? ( (x-omega_min[idx]) / (omega_max[idx] - omega_min[idx]) ) : x
    if(!use_unfolded){
        xpoint = 0.5 * E0
        set arrow from -xpoint, graph 0 to -xpoint,1 nohead ls 1 dt (3,5,10,5) lc rgb 'black' lw 2;
        set arrow from xpoint, graph 0 to xpoint,1 nohead ls 1 dt (3,5,10,5) lc rgb 'black' lw 2;
    }
    set multiplot
    @MARGIN; @SCALE;   plot for[i=i0:iend:di] dir_DOS._name(i) u (omega($1, (i-i0)/di+1)):($2/dos_max[(i-i0)/di+1]) w lp lw 2 t _key_title(i);
    if(!use_unfolded){
        @UNSET; @MARGIN; @SCALE;   plot E_noninteracting u (E_noninteracting[$1]):(0.0) w p ps 3 pt 7 lc rgb "black" notitle
    }
    unset multiplot
} else {
    if(what == 1){
        #set logscale y; set format y '10^{%L}';
        #set logscale x; set format x '10^{%L}'
        set key inside right top
        SCALE = use_logarithmic? "set xrange[-3:2.5];"."set yrange[1e-4:1.1];" : "set xrange[1e-2:2e1]; set yrange[1e-5:1.5];" 
        if(use_logarithmic){ set ylabel 'P(log_{10} {/Symbol w}/{/Symbol w}_H)';    set xlabel 'log_{10} {/Symbol w}/{/Symbol w}_H'}
        else {set ylabel 'P({/Symbol w})';    set xlabel '{/Symbol w}'; };#set logscale x;};
        trueX(x) = use_logarithmic? 10**x : x
        prefactor(x) = use_logarithmic? 10**x * log(10) : 1.0
        omega(x, idx) = use_unfolded? x : (use_logarithmic? x - log10(wH[idx]) : x / wH[idx])
        set multiplot
        @MARGIN; @SCALE;    plot for[i=i0:iend:di] dir_dist._name(i) u (omega($1, (i-i0)/di+1)):(!use_logarithmic && !use_unfolded? $2*wH[(i-i0)/di+1] : $2) w lp lw 2 t _key_title(i),\
                                     '+' u (NaN):(NaN) w p ps 2 pt 7 lc rgb "black" t 'log {/Symbol w}_H^{typ}',\
                                     prefactor(x)*P_GOE(trueX(x)) w l ls 1 dt (2,2) lc rgb 'black' lw 2 t 'GOE',\
                                     prefactor(x)*P_POISSON(trueX(x)) w l ls 1 dt (4,3,2, 2) lc rgb 'red' lw 2 t 'Poisson', '+' u (NaN):(NaN) w l dt (8,8) lw 2 lc 1 t 'log-normal fit'
        @UNSET; 
        @MARGIN; @SCALE;    plot for[i=i0:iend:di] '+' u (wH_typical[(i-i0)/di + 1]):(1e-3) w p ps 2 pt 7 notitle
        if(use_fit){
            @UNSET; 
            @MARGIN; @SCALE;    plot for[i=1:sizee] abs(plot_fit(x, a_list[i], mu_list[i], sig_list[i])) w l dt (8,8) lw 2 notitle
        }
        unset multiplot
    } else {
        if(compare_fits){
            set yrange[-2**(L-4):1.2*2**L]
            param = scaling == 0? 100*J : (scaling == 1? L : 100*w)
            plot for[i=3:6] dir_unfolding._name(param) u 1:(column(i)) w l lw 2 t (i < 5? sprintf("n=%d", 3*(i-2)) : sprintf("n=L+%d",(5*(i-5)))), dir_unfolding._name(param) u 1:2 w steps lw 3 lc rgb "black" t 'cdf' 
        } else {
            SCALE = sprintf("set xrange[-18:23]; set yrange[-200:%.2f];", 1.1 * 2**L)
            if(use_logarithmic){ set ylabel 'P(log {/Symbol w})';    set xlabel 'log {/Symbol w}'}
            else {set ylabel 'P({/Symbol w})';    set xlabel '{/Symbol w}'; set logscale x;};
            set multiplot
            @MARGIN; @SCALE;    plot for[i=i0:iend:di] dir_unfolding._name(i) u 1:2 w steps lw 2 t _key_title(i)
            @UNSET; 
            #@MARGIN; @SCALE;    plot for[i=J0:Jend:dJ] dir_unfolding._name(i) u 1:6 w l dt (8,8) lw 2 notitle
            unset multiplot
        }
    }
}