# 🧪 Atividade 1 — Métodos intervalares

## 🎯 Objetivo

Implementar e analisar os métodos da bissecção e falsa posição, aplicados a uma equação não linear.

---

## 📍 Raiz de referência

x_exato = 14.780208593679468

---

## 📝 Atividades

1. Implementar inicialmente no código fornecido o erro verdadeiro;  
2. Imprimir os registros num arquivo saida.dat;  
3. Plotar com gnuplot a comparação dos erros verdadeiro e relativo para a solução pelo método da bissecção usando o script plot1.gnu; 
4. Implementar no código o método da falsa posição modificada e escrever os registros num arquivo saida2.dat;  
5. Comparar graficamente os dois tipos de erros pelos dois métodos usando o gnuplot com o script plot2.gnu;   

---

## ⚙️ Execução (exemplo)

```bash
gfortran atividade1.f90 
./a.out
gnuplot plot1.gnu
xdg-open comparacao_erros.png
```
