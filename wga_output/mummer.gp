set terminal png tiny size 800,800
set output "wga_output/mummer.png"
set size 1,1
set grid
unset key
set border 0
set tics scale 0
set xlabel "REF"
set ylabel "GL343308.1"
set format "%.0f"
set mouse format "%.0f"
set mouse mouseformat "[%.0f, %.0f]"
set mouse clipboardformat "[%.0f, %.0f]"
set xrange [1:*]
set yrange [1:*]
set style line 1  lt 1 lw 3 pt 6 ps 1
set style line 2  lt 3 lw 3 pt 6 ps 1
set style line 3  lt 2 lw 3 pt 6 ps 1
plot \
 "wga_output/mummer.fplot" title "FWD" w lp ls 1, \
 "wga_output/mummer.rplot" title "REV" w lp ls 2
