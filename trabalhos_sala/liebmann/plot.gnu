# ==========================================================
# CONFIGURAÇÃO GERAL E INÍCIO DO MULTIPLOT
# ==========================================================
set terminal pngcairo size 800,1000 enhanced font 'Arial,12'
set output 'analise_81x81.png'

# Ativa o multiplot com 2 linhas e 1 coluna
set multiplot layout 2,1 title "Analise Termica da Aleta (Metodo de Liebmann)" font ",16"

# ==========================================================
# GRÁFICO 1 (Cima): MAPA DE CALOR 2D
# ==========================================================
set title "Distribuicao de Temperatura" font ",14"
set xlabel "Comprimento x"
set ylabel "Altura y"

set view map
set contour base
set cntrparam levels auto 12
unset surface

# Paleta de cores e configuracoes 3D
set palette defined (0 "dark-blue", 1 "blue", 2 "cyan", 3 "green", 4 "yellow", 5 "red", 6 "dark-red")
set pm3d at b interpolate 2,2

splot 'saida.dat' using 1:2:3 with pm3d notitle, \
      'saida.dat' using 1:2:3 with lines nosurface linecolor rgb "black" notitle

# ==========================================================
# GRÁFICO 2 (Baixo): LINHA CENTRAL
# ==========================================================
# Limpa as configuracoes do mapa de calor para nao bugar o grafico 2D
unset view
unset pm3d
unset contour
unset colorbox
set grid

set title "Perfil de Temperatura na Linha Central" font ",14"
set xlabel "Comprimento x"
set ylabel "Temperatura (C)"

# Limites do eixo Y
set yrange [20:200]
set autoscale x

plot 'centro.dat' using 1:2 with linespoints pt 7 lw 2 lc rgb "red" title "Temperatura no Centro"

# Finaliza a montagem da imagem
unset multiplot