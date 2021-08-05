set terminal post enhanced colour solid font 25
set encoding iso_8859_2
 
set style line 1 lw 0.5 ps 0.6 lc rgb "red"              


set xlabel "t"         
set ylabel "x(t)"
set grid
set size square
set yrange [-1:1]
set key vertical opaque outside top center

set output "x1.eps"
plot "out.dat" i 0 w lp ls 1 t "x(t); {/Symbol b} = 0, F_0 = 0, {/Symbol W} = 0.8"

set output "x2.eps"
plot "out.dat" i 1 w lp ls 1 t "x(t); {/Symbol b} = 0.4, F_0 = 0, {/Symbol W} = 0.8"

set output "x3.eps"
plot "out.dat" i 2 w lp ls 1 t "x(t); {/Symbol b} = 0.4, F_0 = 0.1, {/Symbol W} = 0.8"


# INFORMACJE NA TEMAT PLIKU Z DANYMI:
#
# Plik: out.dat zawiera wszystkie wyniki:
#          - W pierwszej kolumnie: t_i -- czas
#          - W drugiej kolumnie:   x_i -- wychylenie w chwili czasowej t_i
# Poszczegolne serie danych sa oddzielone dwiema pustymi liniami. Pierwsza seria danych odpowiada pierwszemu zestawowi parametrow beta, F0, Omega; druga: drugiemu zestawowi, itd.

# WYJASNIENIE OPCJI KOMENDY plot:
#          i 0 -- uzycie pierwszej serii danych z pliku; (i 1 oraz i 2 -- analogicznie przy uzyciu drugiej oraz trzeciej serii danych)
#          w lp  == with linespoints  -- wykres bedzie narysowany przy pomocy punktow polaczonych linia ciagla.
#          ls 1  == linestyle 1       -- wybor stylu linii
#          t "tekst" == title "tekst" -- "tekst" bedzie podpisem na legendzie
