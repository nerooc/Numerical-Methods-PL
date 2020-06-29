set terminal post enhanced colour solid font 20  

reset

set style line 1 lt rgb "red" lw 2 dt 2
set style line 2 lt rgb "royalblue" lw 2

# PLIK w0.dat powinien zawierac jedna kolumne:
#	wartosc funkcji dla wybranego wedrowca (np. pierwszego) w kazdym z krokow algorytmu.

# PLIK T.dat powinien zawierac trzy serie danych z dwiema kolumnami:
#	- pierwsza seria danych: polozenia x_i, y_i wszystkich wedrowcow dla przypadku i_T = 0
#	- druga seria danych: polozenia x_i, y_i wszystkich wedrowcow dla przypadku i_T = 7
#	- trzecia seria danych: polozenia x_i, y_i wszystkich wedrowcow dla przypadku i_T = 20

set yrange [*:*] # autoskalowanie
set output "w0.eps"
set multiplot
set xlabel "Numer kroku algorytmu"
set ylabel "f(x_0, y_0)"

do for [i=1:20] {
	set arrow i from 100*i, graph 0 to 100*i, graph 1 nohead ls 1 back
}
plot [:2100] "w0.dat" u 1 w l ls 2 notitle


set origin 0.35, 0.4
set ylabel offset 1, 0
set size 0.6, 0.5
set xtics 500
do for [i=1:20] {
	unset arrow i
}
set obj 1 rect from graph 0, graph 0 to graph 1, graph 1 fc rgb "white" fs solid 1.0 noborder 
plot [1000:2100] "w0.dat" u 1 w l ls 2 notitle
unset multiplot


set terminal post enhanced colour solid font 20
f(x,y) = sin(x)*sin(y) - exp(-(x+pi/2.0)**2-(y-pi/2.0)**2)

reset
set xrange [-10:10]
set yrange [-10:10]
set isosample 250, 250
set table "test.dat"
splot f(x,y)
unset table

set contour base
set cntrparam levels 20
unset surface
set table "cont.dat"
splot f(x,y)
unset table

set size square
unset key
set xlabel "x"
set ylabel "y"
set cblabel "f(x,y)" rotate by -90
set xtics out scale 2
set ytics out scale 2

i = 0
do for [i=0:2] {
	set output "T".i.".eps"
	p "test.dat" with image, "cont.dat" w l lt -1 lw 1.5, "T.dat" i i u 1:2 w p lc rgb "black" lw 4 ps 1.5
}

