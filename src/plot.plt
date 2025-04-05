set terminal png
set output 'plot.png'
plot 'data.txt' with lines title 'y=sin(x)'
