# ==============================================================================
# SCRIPT GNUPLOT - Problema Gravitacional de Dois Corpos
# ==============================================================================

# Configurações gerais de estilo e formato de saída
set terminal pngcairo size 800,600 enhanced font 'Arial,12'
set grid

# ------------------------------------------------------------------------------
# 1. Gráfico da Trajetória no Plano (x, y)
# ------------------------------------------------------------------------------
set output 'trajetoria.png'
set title 'Trajetória dos Corpos no Plano (x, y)'
set xlabel 'Posição x'
set ylabel 'Posição y'

# Plotando as trajetórias usando os resultados do método Leapfrog
plot 'leapfrog.dat' using 2:3 with lines linewidth 1.5 title 'Corpo 1', \
     'leapfrog.dat' using 6:7 with lines linewidth 1.5 title 'Corpo 2'

# ------------------------------------------------------------------------------
# 2. Gráfico do Espaço de Fases (x, vx)
# ------------------------------------------------------------------------------
set output 'espaco_fases.png'
set title 'Espaço de Fases'
set xlabel 'Posição x'
set ylabel 'Velocidade v_x'

# Plotando o espaço de fases x versus vx para cada corpo (método Leapfrog)
plot 'leapfrog.dat' using 2:4 with lines linewidth 1.5 title 'Corpo 1 (x_1 vs v_{x1})', \
     'leapfrog.dat' using 6:8 with lines linewidth 1.5 title 'Corpo 2 (x_2 vs v_{x2})'

# ------------------------------------------------------------------------------
# 3. Gráfico da Conservação de Energia E(t)
# ------------------------------------------------------------------------------
set output 'energia.png'
set title 'Evolução da Energia Total do Sistema'
set xlabel 'Tempo (t)'
set ylabel 'Energia Mecânica (E)'

# Desabilitando o formato científico forçado no eixo Y para melhorar a leitura de variações pequenas
set format y "%g"

# Plotando a energia em função do tempo para comparar os dois métodos
plot 'euler.dat' using 1:10 with lines linewidth 1.5 title 'Método de Euler', \
     'leapfrog.dat' using 1:10 with lines linewidth 1.5 title 'Método Leapfrog'

print "Gráficos gerados com sucesso: trajetoria.png, espaco_fases.png e energia.png"