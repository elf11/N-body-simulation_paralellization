set terminal png nocrop enhanced font arial 8 size 500,350 
set output 'histograms.2.png'
set bar 1.000000
set style rectangle back fc lt -3 fillstyle  solid 1.00 border -1
set key inside right top vertical Right noreverse enhanced autotitles columnhead nobox
set datafile missing '-'
set style data linespoints
set xtics border in scale 1,0.5 nomirror rotate by -45  offset character 0, 0, 0
set xtics  norangelimit
set title "Comparatie intre variantele seriala/omp/pthreads/mpi" 
set rrange [ * : * ] noreverse nowriteback  # (currently [0.00000:10.0000] )
set trange [ * : * ] noreverse nowriteback  # (currently [-5.00000:5.00000] )
set urange [ * : * ] noreverse nowriteback  # (currently [-5.00000:5.00000] )
set vrange [ * : * ] noreverse nowriteback  # (currently [-5.00000:5.00000] )
set ylabel  offset character 0, 0, 0 font "" textcolor lt -1 rotate by 90
set y2label  offset character 0, 0, 0 font "" textcolor lt -1 rotate by 90
set cblabel  offset character 0, 0, 0 font "" textcolor lt -1 rotate by 90
set locale "C"
set key outside
plot for [col=1:5] 'in2.dat' using 0:col with lines title columnheader

