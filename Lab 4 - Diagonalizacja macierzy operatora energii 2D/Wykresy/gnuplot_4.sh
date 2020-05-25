set term png enhanced font Arial 17
set view map
set pm3d interpolate 4,4
set size ratio -1
set xrange [1:20]
set yrange [1:20]
set xtics out offset 0,0.8
set ytics out offset 1,0
unset key
set xlabel "x" offset 0,0.8
set ylabel "y" offset 1,0

do for [NUMER=3:12] {
	k=NUMER-2
	set cblabel "{/Symbol y_{".k."}}(x,y)" rotate by -90
	set out "vec_".k.".png"
	splot "dane.dat" i 0 u 1:2:NUMER w pm3d
}
