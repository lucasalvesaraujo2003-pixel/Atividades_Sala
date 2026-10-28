import matplotlib.pyplot as plt
import numpy as np

# 1. Gerar a grade de pontos (mesmo processo do 3D)
x = np.linspace(-5, 5, 100)
y = np.linspace(-5, 5, 100)
X, Y = np.meshgrid(x, y)
R = np.sqrt(X**2 + Y**2)
Z = np.sin(R)

# 2. Criar uma janela com 2 gráficos (lado a lado)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# --- GRÁFICO 1: Apenas as Linhas de Nível ---
# O parâmetro 'levels=15' define quantas camadas/divisões queremos mostrar
linhas = ax1.contour(X, Y, Z, levels=15, cmap='viridis')
ax1.clabel(linhas, inline=True, fontsize=8) # Adiciona os números da altura nas linhas
ax1.set_title("Linhas de Nível (contour)")
ax1.set_xlabel("Eixo X")
ax1.set_ylabel("Eixo Y")

# --- GRÁFICO 2: Camadas Preenchidas com Cor ---
camadas = ax2.contourf(X, Y, Z, levels=15, cmap='coolwarm')
fig.colorbar(camadas, ax=ax2, label="Valores de f(x, y)") # Barra de cor de referência
ax2.set_title("Camadas Preenchidas (contourf)")
ax2.set_xlabel("Eixo X")
ax2.set_ylabel("Eixo Y")

plt.tight_layout()
plt.show()
