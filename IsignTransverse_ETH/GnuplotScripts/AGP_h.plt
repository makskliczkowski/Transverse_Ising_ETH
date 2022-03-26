reset
set term qt size 900, 800 font ",12"
set mxtics
set mytics
set style line 12 lc rgb '#ddccdd' lt 1 lw 1.5
set style line 13 lc rgb '#ddccdd' lt 1 lw 0.5
set grid xtics mxtics ytics mytics back ls 12, ls 13

set logscale y
set size square
set xtics nomirror
set xtics auto
set ytics auto

dir_base='../../results/disorder/PBC/'
dir = dir_base.'AGP/'

L = 14
w = 0.01
g = 0.05
	glist = '0.025 0.10 0.20 0.30 0.40 0.50 0.70'
	glist2 = '0.05 0.10 0.15'
AGP_or_susc = 1						# 0-typ susc; 1-AGP
SigX_or_SigZ = 1 					# 0-SigX as sum, 1-SigZ as sum

L_min = 8							# minimual system size in file
L_m = 10							# minimal size in scaling
num = 7   							# number of system sizes in file
f(i)=(i-2 - AGP_or_susc)/2 + L_min	# system size taken from column id

scaling = 1							# size scaling=1 or g-scaling=0 or site-scaling=2
norm = 0							# normalize by value at h=0?
size_parity = 0						# 0-all sizes; 1-only odd size; 2-only even sizes
q_vs_j = 0
site = 3

# GRAPH DESIGN
my_title = (AGP_or_susc == 1)? sprintf("{/*1.2 AGP with disorder w = %.2f}", w)."\n{/*1.2 Operator: O=}" : \
			 sprintf("{/*1.2 Typical susceptibility with disorder w = %.2f}", w)."\n{/*1.2 Operator: O=}"
str = SigX_or_SigZ == 0? "x" : "z"
my_title = my_title.(q_vs_j? "{/*1.2 L^{-1/2}{/Symbol S}_j e^{iqj} {/Symbol s}^".str."_j}"\
						: "{/*1.2 {/Symbol s}^".str."_".(scaling == 2? "j}" : sprintf("{%d}}", site)))
my_title = my_title.(scaling? sprintf("{/*1.2  for g = %.3f}", g) : sprintf("{/*1.2  for L = %d}", L))
set title my_title

set key outside right

set ylabel (AGP_or_susc == 1)? '{/*1.3|{/Symbol L}_O|^2'.(norm? "/ |{/Symbol L}_O(0)|^2" : " / D")\
 : '{/*1.3 {/Symbol c}_O^{typ}'.(norm? "/ {/Symbol c}_O^{typ}(0)" : " / D")
set xlabel '{/*1.3h}'


## PLOT
name = dir.(SigX_or_SigZ == 0? 'SigmaX' : 'SigmaZ')
plotname(x, s) = name."_".(q_vs_j? "q" : "j").sprintf("=%d",s).sprintf(",J0=%.2f,g=%.2f,g0=%.2f,w=%.2f.dat", 0, x, 0, w);

array val[num]
do for[i=AGP_or_susc+2:2*num+1:2]{
	stats plotname(g,site) every ::0::0 using i nooutput
	val[f(i) - L_min + 1] = (norm ? STATS_min : 2**f(i))
}

#set yrange[0.01:5]
#set xrange[0:3]	
if(scaling==1){
	i0 = 2*(L_m - L_min + 1) + AGP_or_susc
	if(size_parity > 0){
		i0 = i0 + (( (size_parity == 1)^(L_min % 2) )? 2 : 0)
	}
	imax = 2*num+1
	di = 2 + (size_parity > 0? 2 : 0)
	plot for[i=i0:imax:di] plotname(g, site) u 1:(column(i) / val[f(i) - L_min + 1]) w lp lc abs(f(i)-L) pt 7 ps 0.7 title sprintf("L=%d", f(i))
} else{
	if(scaling == 0){
		i = 2*(L - L_min + 1) + AGP_or_susc
		plot for [gx in glist2] plotname(1.*gx, site) u 1:(column(i) / val[f(i) - L_min + 1]) w lp pt 7 ps 0.7 title "g=".gx
	} else{
		i = 2*(L - L_min + 1) + AGP_or_susc
		plot for [s=0:L] plotname(g, s) u 1:(column(i) / val[f(i) - L_min + 1]) w lp pt 7 ps 0.7 title sprintf("j=%d",s)
	}
}
