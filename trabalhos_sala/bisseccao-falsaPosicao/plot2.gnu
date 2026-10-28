set terminal pngcairo size 800,600
set output 'comparacao_metodos.png'

set title 'Comparação de Métodos'
set xlabel 'Iteração'
set ylabel 'Erro'
set grid
set logscale y

plot \
    'saida.dat'  using 1:5 with lines lw 2 title 'Erro Relativo (Bissecção)', \
    'saida.dat'  using 1:6 with lines lw 2 title 'Erro Verdadeiro (Bissecção)', \
    'saida2.dat' using 1:5 with lines lw 2 dt 2 title 'Erro Relativo (Falsa Posição)', \
    'saida2.dat' using 1:6 with lines lw 2 dt 2 title 'Erro Verdadeiro (Falsa Posição)'
