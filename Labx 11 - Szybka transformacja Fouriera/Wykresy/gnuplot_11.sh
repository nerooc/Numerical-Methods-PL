set terminal post enhanced colour solid font 25 

l_width = 5
p_size = 1.3
set style line 1 lt rgb "black" lw 1 pt 1 ps 1
set style line 2 lt rgb "black" lw l_width*1.2 pt 7 ps 0.7
set style line 3 lt rgb "red" lw 1 pt 2 ps 1
set style line 4 lt rgb "red" lw l_width*0.6 pt 7 ps 0.5 
set style line 5 lt rgb "gray50" lw 1 pt 1 ps 1
T = 1.0
omega = 2.*pi/T
f_0(t) = sin(omega*t)+sin(2.*omega*t)+sin(3.*omega*t)
set xlabel "t"
set ylabel "f(t)"
set samples 500
set grid
set key outside top Right horizontal samplen 2 font ",18"

# PETLA DLA WYKRESOW WYNIKOW ODSZUMIANIA wedlug plikow k8.dat, k10.dat, k12.dat.
# W kazdym z nich (dla odpowiedniego k podanego w nazwie pliku) sa dwie serie danych:
#        - 1. seria: oryginalny, zaszumiony sygnal f wedlug wzoru (4) w dwoch kolumnach: t_i, f[2*i].
#        - 2. seria (po dwoch pustych liniach): odszumiony sygnal, 
#              czyli tablica f po zakonczeniu wszystkich krokow programu w dwoch kolumnach: t_i, y[2*i].
# Nie trzeba wypisywac elementow o nieparzystych indeksach (czesci urojone powinny znikac).


do for [k = 8:12:2]{
    set output "FFT_k".k.".eps"
    plot [0:3.][] "k".k.".dat" i 0 u 1:2 w lp ls 5 t "f(t_i) = f_0(t_i)+{/Symbol D}",\
     f_0(x) ls 2 t "f_0(t)",\
     "k".k.".dat" i 1 u 1:2 w lp ls 4 t "f(t_i)*g(t_i)"
}
