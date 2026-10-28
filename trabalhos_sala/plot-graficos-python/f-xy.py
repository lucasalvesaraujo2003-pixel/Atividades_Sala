# Plot de gráfico para um função f(x, y)

import matplotlib.pyplot as plt
import numpy as np

# 1. Definir os intervalos de x e y
x = np.linspace(-5, 5, 100)
y = np.linspace(-5, 5, 100)

# 2. Criar a grade (mesh) combinando todos os pontos de X e Y
X, Y = np.meshgrid(x, y)

# 3. Definir a função f(x, y) usando as matrizes da grade
# Exemplo: f(x, y) = sen(sqrt(x² + y²))
R = np.sqrt(X**2 + Y**2)
Z = np.sin(R)

# 4. Configurar a figura para projeção 3D
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')

# 5. Plotar a superfície
# 'cmap' define a paleta de cores (viridis, plasma, coolwarm, etc.)
surf = ax.plot_surface(X, Y, Z, cmap='coolwarm', edgecolor='none', alpha=0.9)

# 6. Customizar o gráfico
ax.set_title("Gráfico 3D da Função f(x, y)", fontsize=14)
ax.set_xlabel("Eixo X", fontsize=10)
ax.set_ylabel("Eixo Y", fontsize=10)
ax.set_zlabel("Eixo Z - f(x,y)", fontsize=10)

# Adiciona uma barra lateral de cores para referência dos valores de Z
fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5)

# 7. Exibir o gráfico
plt.show()
