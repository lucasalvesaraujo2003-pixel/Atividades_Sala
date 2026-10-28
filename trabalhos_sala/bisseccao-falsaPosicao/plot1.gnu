set terminal pngcairo size 800,600
set output 'comparacao_erros.png'

set title 'Comparação de Erros - Método da Bissecção'
set xlabel 'Iteração'
set ylabel 'Erro'
set grid

set logscale y

plot \
    'saida.dat' using 1:5 with lines lw 2 title 'Erro Relativo (Bissecção)', \
    'saida.dat' using 1:6 with lines lw 2 title 'Erro Verdadeiro (Bissecção)'
