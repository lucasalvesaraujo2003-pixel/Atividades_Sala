# Configurações gerais de estilo e saída
set terminal pngcairo size 1200, 500 enhanced font 'Verdana,11'
set output 'comparativo_metodos.png'

# Define que o gráfico será dividido em 2 subplots (1 linha, 2 colunas)
set multiplot layout 1, 2 title "Análise Numérica: Sistema Massa-Mola-Amortecedor Forçado" font 'Verdana,14 B'

# Estilo das linhas para ficar mais visível
set style line 1 lc rgb '#E63946' lt 1 lw 2  # Vermelho para Euler
set style line 2 lc rgb '#1D3557' lt 1 lw 2  # Azul Escuro para Heun
set grid

# -----------------------------------------------------------------
# GRÁFICO 1: Posição vs Tempo (x vs t)
# -----------------------------------------------------------------
set title "Posição ao longo do Tempo (x vs t)"
set xlabel "Tempo (s)"
set ylabel "Posição x (m)"

# No arquivo, Coluna 1 = t, Coluna 2 = x, Coluna 3 = v
plot 'euler.dat' using 1:2 with lines title "Euler" ls 1, \
     'heun.dat'  using 1:2 with lines title "Heun" ls 2

# -----------------------------------------------------------------
# GRÁFICO 2: Espaço de Fase (x vs v)
# -----------------------------------------------------------------
set title "Espaço de Fase (Posição vs Velocidade)"
set xlabel "Posição x (m)"
set ylabel "Velocidade v (m/s)"

# Aqui plotamos a Coluna 2 (x) no eixo X e a Coluna 3 (v) no eixo Y
plot 'euler.dat' using 2:3 with lines title "Euler" ls 1, \
     'heun.dat'  using 2:3 with lines title "Heun" ls 2

unset multiplot