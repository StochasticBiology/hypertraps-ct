# plot summarised posterior ordering from verification study
mkdir Plots
cd VerifyData/
gnuplot -e 'set term svg size 640, 480; set xlabel "Feature"; set ylabel "Ordering"; set output "../Plots/plot-cross-0.svg"; plot "synth-cross-samples-0.txt-posterior-1-1-2-5-0.txt.process" u 3:1:(sqrt($5)*5) ps variable pt 7; exit'
gnuplot -e 'set term svg size 640, 480; set xlabel "Feature"; set ylabel "Ordering"; set output "../Plots/plot-cross-1.svg"; plot "synth-cross-samples-1.txt-posterior-1-1-2-5-0.txt.process" u 3:1:(sqrt($5)*5) ps variable pt 7; exit'
gnuplot -e 'set term svg size 640, 480; set xlabel "Feature"; set ylabel "Ordering"; set output "../Plots/plot-cross-2.svg"; plot "synth-cross-samples-2.txt-posterior-1-1-2-5-0.txt.process" u 3:1:(sqrt($5)*5) ps variable pt 7; exit'

