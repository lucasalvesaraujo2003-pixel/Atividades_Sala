! ==============================================================================
! PROGRAMA: inversa (Cálculo da Matriz Inversa via Decomposição LU)
!
! 1. RESUMO DO MÉTODO NUMÉRICO:
!    O programa calcula a inversa de uma matriz quadrada [A] utilizando o 
!    método de Decomposição LU. O processo baseia-se no princípio de que 
!    [A]*[A_inv] = [I], onde [I] é a matriz identidade. O algoritmo opera em 3 passos:
!    - Fatoração: A matriz [A] é decomposta nas matrizes triangulares [L] e [U].
!    - Resolução Múltipla: Para cada coluna da matriz identidade (tratada como um
!      vetor independente {e}), o programa resolve o sistema [L]*[U]*{x} = {e}
!      usando substituição progressiva para encontrar {d} e regressiva para {x}.
!    - Construção: Cada vetor solução {x} calculado compõe uma coluna exata da 
!      matriz inversa [A_inv].
!
! 2. PARA QUE SERVE (APLICAÇÃO GERAL DO MÉTODO):
!    Encontrar a inversa de uma matriz é uma das operações mais centrais da álgebra 
!    linear aplicada. Na engenharia e ciências, a matriz inversa é frequentemente 
!    utilizada para encontrar a matriz de flexibilidade em estruturas (que é a 
!    inversa da matriz de rigidez), em regressões estatísticas (Mínimos Quadrados), 
!    em sistemas de controlo multivariáveis, processamento de imagens e computação 
!    gráfica. Utilizar a Decomposição LU para inverter matrizes é o padrão da 
!    indústria, pois é computacionalmente muito mais rápido, estável e eficiente do 
!    que métodos clássicos (como a Regra de Cramer) para sistemas maiores que 3x3.
!
! 3. INPUTS (ENTRADAS DEFINIDAS NO CÓDIGO):
!    - n : Dimensão da matriz (ordem 10, definida no código).
!    - A : Matriz quadrada de coeficientes 10x10 explicitamente declarada.
!
! 4. OUTPUTS (SAÍDAS GERADAS):
!    - A matriz inversa [A_inv] é armazenada na memória.
!    - O programa calcula o produto ident_calc = matmul(A,A_inv) para fins de 
!      validação teórica.
!    - No estado atual, o programa exibe no terminal (output padrão) as linhas 
!      da matriz identidade teórica ('ident'). 
!      (Nota: Há um comentário no código sugerindo a futura implementação do 
!      cálculo de erro entre a identidade calculada e a real).
! ==============================================================================

program inversa
implicit none

integer i, j, n, k
real(8) A(10,10), L(10,10), U(10,10), x(10), d(10), vetor(10)
real(8) A_inv(10,10), ident_calc(10,10), ident(10,10)
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

! Calculando a inversa de [A]

do k = 1, n
	vetor = 0.0
	vetor(k) = 1.0
	
	! Encontrando o vetor {d}
	d(1) = vetor(1)
	do i = 2, n
		soma = vetor(i)
		do j = 1, i-1
			soma = soma - L(i,j)*d(j)
		end do
		d(i) = soma
	end do
	
	! Encontrando o vetor {x}
	x(n) = d(n)/U(n,n)		  ! començando de baixo
	do i = n-1, 1, -1
		soma = d(i)			  ! subindo de baixo para cima				
		do j= i + 1,n
			soma = soma - U(i,j)*x(j)  ! soma
		end do
		x(i) = soma/U(i,i)	  ! substituindo em {X}
	end do
	
	! Armazendo {x} como as colunas de [A_inv]
	do i = 1, n
		A_inv(i,k) = x(i)
	end do
end do

write(*,*) "Matriz inversa encontrada:"
write(*,*) " "
do i = 1, n
	write(*,*) A_inv(i,:)
end do

!calcular o erro entre a identidade calculada e a identidade real
ident_calc = matmul(A,A_inv)

ident = 0.0
do i = 1, n
	ident(i,i) = 1.0
end do

!alterar o código para diferentes ordens de [A]

end program inversa
