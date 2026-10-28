! ==============================================================================
! PROGRAMA: lu (Resolução de Sistemas Lineares)
!
! 1. RESUMO DO MÉTODO NUMÉRICO:
!    O programa implementa o Método de Decomposição LU. A rotina fatora 
!    a matriz de coeficientes [A] no produto de duas matrizes distintas: uma matriz 
!    triangular inferior [L] e uma matriz triangular superior [U]. 
!    A resolução do sistema linear [A]{X} = {B} ocorre em três etapas principais:
!    - Fatoração: Decompõe [A] de forma que [L]*[U] = [A].
!    - Substituição Progressiva: Resolve o sistema triangular inferior [L]{D} = {B} 
!      para encontrar um vetor intermediário {D} (de cima para baixo).
!    - Substituição Regressiva: Resolve o sistema triangular superior [U]{X} = {D} 
!      para determinar o vetor solução final {X} (de baixo para cima).
!    Além disso, o algoritmo calcula o erro associado à decomposição e à solução final.
!
! 2. PARA QUE SERVE (APLICAÇÃO GERAL DO MÉTODO):
!    A Decomposição LU é um método de extrema importância na álgebra linear computacional.
!    Sua principal vantagem sobre a Eliminação de Gauss clássica é a separação das 
!    operações na matriz [A] do vetor independente {B}. Isso a torna extremamente 
!    eficiente em simulações de engenharia onde a mesma matriz de propriedades (como
!    uma matriz de rigidez estrutural, condutividade térmica ou rede elétrica) 
!    permanece constante, mas precisa ser avaliada para múltiplos vetores de carga/força
!    {B} diferentes. Também é o algoritmo padrão para calcular inversas de matrizes e 
!    determinantes de forma computacionalmente eficiente.
!
! 3. INPUTS (ENTRADAS DEFINIDAS NO CÓDIGO):
!    - n : Dimensão do sistema linear (10 variáveis/equações).
!    - A : Matriz quadrada de coeficientes 10x10 explicitamente definida.
!    - B : Vetor de termos independentes com 10 elementos predefinidos.
!
! 4. OUTPUTS (SAÍDAS GERADAS):
!    Os resultados de validação são exibidos no terminal (output padrão):
!    - O erro absoluto total acumulado (somaErro) entre a matriz original [A] e o 
!      produto calculado [L]*[U].
!    - O erro absoluto total acumulado entre o vetor independente original [B] e 
!      o vetor obtido por [A]*{X}.
! ==============================================================================

program lu
implicit none

integer i, j, n, k
real(8) A(10,10), B(10), X(10), D(10), Bcalc(10)
real(8) L(10,10), U(10,10), Acalc(10,10)
real fator, soma, somaErro, erroA, erroB

n = 10

! Matriz [A]
A = reshape((/ &
  3.14, 7.25, 1.83, 9.42, 5.67, 2.11, 8.93, 4.56, 6.78, 0.99, &
  7.44, 2.38, 6.15, 1.27, 8.61, 3.90, 5.72, 9.18, 4.03, 0.56, &
  2.74, 8.36, 5.91, 7.13, 1.49, 6.82, 3.57, 9.65, 4.21, 0.88, &
  6.47, 1.92, 8.24, 3.76, 5.33, 9.07, 2.68, 7.55, 4.89, 0.14, &
  9.31, 4.62, 2.05, 8.79, 6.18, 1.73, 7.96, 3.28, 5.44, 0.67, &
  5.87, 9.14, 3.41, 6.53, 2.96, 8.22, 1.35, 7.69, 4.58, 0.72, &
  1.64, 7.88, 4.11, 9.26, 5.39, 2.57, 8.03, 6.91, 3.75, 0.48, &
  8.54, 3.19, 6.77, 1.95, 9.62, 4.36, 7.08, 2.43, 5.29, 0.81, &
  4.97, 8.16, 2.84, 7.31, 1.58, 6.45, 9.73, 3.66, 5.12, 0.25, &
  7.59, 2.71, 5.48, 9.04, 3.22, 8.67, 1.86, 6.34, 4.95, 0.39  &
/), (/10,10/))

! Definindo [U] como cópia de [A] para não fazer alterações na matriz original
U = A

! Transformando a matriz [L] em identidade
L = 0.0
do i = 1,n
   L(i,i) = 1.0
end do

! Decomposição de [A] em [L] e [U]
do k = 1,n-1					! linha anterior
   do i=k+1,n					! linha atual
      fator=U(i,k)/U(k,k)			! fator
      L(i,k) = fator				! armazenando o fator na matriz [L]
       do j=k+1,n				! colunas
          U(i,j) = U(i,j) - fator*U(k,j)	! eliminando [A]
          U(i,k) = 0.0
       end do
   end do
end do

! Verificando se [L]*[U] = [A]
Acalc = matmul(L,U)		  ! Matriz [A] calculada
somaErro = 0.0
do i=1,n
	do j=1,n
		somaErro = somaErro + abs(A(i,j) - Acalc(i,j))
	end do
end do
erroA = somaErro/(n*n)
write(*,*) 'Erro entre as matrizes [A] e [A*]:', somaErro

! Vetor resultados {B}
B = (/ 4.25, 7.81, 1.93, 9.14, 3.67, 6.52, 2.48, 8.76, 5.39, 0.95 /)

! Substituição para encontrar o vetor {D}
D(1) = B(1)			  ! començando de cima
do i=2,N
   soma = B(i)			  ! descendo de cima para baixo				
    do j= 1, i-1
       soma = soma - L(i,j)*D(j)  ! soma
    end do
       D(i) = soma		  ! substituindo em {D}
end do

! Substituição para encontrar o vetor {X}
X(N) = D(N)/U(N,N)		  ! començando de baixo
do i=N-1,1, -1
   soma = D(i)			  ! subindo de baixo para cima				
    do j= i+1,N
       soma = soma - U(i,j)*X(j)  ! soma
    end do
       X(i) = soma/U(i,i)	  ! substituindo em {X}
end do

Bcalc = matmul(A,X)
somaErro = 0.0
do i=1,n
	somaErro = somaErro + abs(B(i) - Bcalc(i))
end do
erroB = somaErro/n
write(*,*) 'Erro entre as matrizes [B] e [B*]:', somaErro

end program lu
