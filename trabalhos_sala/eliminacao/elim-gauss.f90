! ==============================================================================
! PROGRAMA: eliminacao (Resolução de Sistemas Lineares)
!
! 1. RESUMO DO MÉTODO NUMÉRICO:
!    O programa implementa o Método de Eliminação de Gauss clássico para resolver
!    sistemas de equações lineares do tipo [A]{X} = {B}. O algoritmo opera em
!    duas etapas principais:
!    - Eliminação Progressiva: Transforma a matriz original de coeficientes [A]
!      em uma matriz triangular superior utilizando operações elementares nas 
!      linhas. O vetor de termos independentes {B} é modificado simultaneamente.
!    - Substituição Regressiva: Resolve as equações de baixo para cima (do último
!      termo X(N) até o primeiro X(1)), aproveitando a estrutura triangular da
!      matriz resultante para isolar as incógnitas uma a uma.
!
! 2. PARA QUE SERVE (APLICAÇÃO GERAL DO MÉTODO):
!    A Eliminação de Gauss é um método direto fundamental para resolver sistemas
!    de equações lineares simétricos ou assimétricos. Na engenharia e ciências 
!    computacionais, serve como base para modelar e resolver problemas reais de 
!    grande porte, tais como: análise estrutural (cálculo de tensões em treliças
!    e pórticos), balanço de massa e energia em engenharia química, circuitos 
!    elétricos complexos (leis de Kirchhoff), redes de distribuição de fluidos, 
!    e na resolução numérica de equações diferenciais via Método dos Elementos 
!    Finitos (MEF).
!
! 3. INPUTS (ENTRADAS DEFINIDAS NO CÓDIGO):
!    - n    : Dimensão do sistema linear (Número de equações e incógnitas: 3).
!    - A    : Matriz quadrada de coeficientes de ordem 3x3 (A(1,1) até A(3,3)).
!    - B    : Vetor coluna de termos independentes com 3 elementos (B(1) a B(3)).
!
! 4. OUTPUTS (SAÍDAS GERADAS):
!    Os resultados finais calculados são impressos diretamente na tela (terminal):
!    - Exibição formatada com 3 casas decimais contendo as respostas do vetor {X}:
!      X(1) = [Valor]
!      X(2) = [Valor]
!      X(3) = [Valor]
! ==============================================================================

program eliminacao
implicit none

integer i, j, n, k
real(8) A(3,3), B(3), X(3)
real fator, soma
n = 3

A(1,1) = 1
A(1,2) = 1
A(1,3) = 1
A(2,1) = 1
A(2,2) = 1.0001
A(2,3) = 1
A(3,1) = 1
A(3,2) = 1
A(3,3) = 1.0001

B(1) = 3
B(2) = 3.0001
B(3) = 3.0001

! Eliminação progressiva
do k = 1,n-1					! linha anterior
   do i=k+1,n					! linha atual
      fator=A(i,k)/A(k,k)			! fator
       do j=k+1,n				! colunas
          A(i,j) = A(i,j) - fator*A(k,j)	! eliminando [A]
          A(i,k) = 0.0
       end do
	B(i) = B(i) - fator*B(k)		! alterando {B}
   end do
end do

! Substituição regressiva
X(N) = B(N)/A(N,N)		  ! començando de baixo
do i=N-1,1, -1
   soma = B(i)			  ! subindo de baixo para cima				
    do j= i+1,N
       soma = soma - A(i,j)*X(j)  ! soma
    end do
       X(i) = soma/A(i,i)	  ! substituindo em {X}
end do

do i = 1, n
   write(*,'(A,I1,A,F8.3)') 'X(', i, ') = ', X(i)
end do

end program eliminacao
