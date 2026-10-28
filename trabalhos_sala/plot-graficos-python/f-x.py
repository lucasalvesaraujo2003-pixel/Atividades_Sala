# Plot de gráfico para uma função f(x)

import matplotlib.pyplot as plt
import numpy as np

# 1. Definir o intervalo de x (ex: de -10 a 10, gerando 400 pontos no meio)
x = np.linspace(-10, 10, 400)

# 2. Definir a sua função de x (ex: f(x) = x² - 3x + 2)
y = x**2 - 3*x + 2

# 3. Criar o gráfico
plt.figure(figsize=(8, 5))  # Define o tamanho da janela
plt.plot(x, y, label="f(x)", color="blue", linewidth=2)

# 4. Customizar o gráfico (títulos, eixos, etc.)
plt.title("Gráfico da Função f(x)", fontsize=14)
plt.xlabel("Eixo X", fontsize=12)
plt.ylabel("Eixo Y", fontsize=12)

# Adiciona linhas de grade e eixos de referência (0,0)
plt.grid(True, linestyle="--", alpha=0.6)
plt.axhline(0, color="black", linewidth=0.8)
plt.axvline(0, color="black", linewidth=0.8)

# Mostra a legenda que definimos no plt.plot
plt.legend()

# 5. Exibir o gráfico na tela
plt.show()
