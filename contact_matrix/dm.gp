set palette defined ( -1 "white",-0.2 "orange",0.6  "red", 1 "black" )
#set palette defined ( -1 "blue", 0 "white", 1 "red")
#set palette negative
#set palette rgb -34,-35,-36
unset key
#set cbrange [0.05:0.3]
set yrange [0:1000] reverse
set xrange [0:1000]
set logscale cb
#set terminal postscript eps enhanced color font 'Helvetica,12'
set terminal pdfcairo
set output 'cm.pdf'
set size square
p 'mytraj.lammpstrj0_cm.dat' matrix using 1:2:3 with image
#p 'smooth_div_cmap.txt' matrix using 1:2:3 with image
#p 'map_div.txt' matrix using 1:2:3 with image
#p 'diff_hp1-wt.txt' matrix using 1:2:3 with image
