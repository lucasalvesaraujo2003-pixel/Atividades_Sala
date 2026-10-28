! ==============================================================================
! PROGRAMA: interpolacao (Otimização Unidimensional)
!
! 1. RESUMO DO MÉTODO NUMÉRICO:
!    O programa implementa o Método da Interpolação Quadrática para otimização
!    unidimensional (encontrar o ponto de máximo ou mínimo de uma função). 
!    O algoritmo não exige o cálculo de derivadas. Ele aproxima a curva da função 
!    original por um polinómio de segundo grau (uma parábola) que passa por três 
!    pontos iniciais (x0, x1 e x2). Em seguida, calcula analiticamente o vértice 
!    dessa parábola (x3). Dependendo do valor da função neste novo ponto, o pior 
!    ponto do conjunto anterior é descartado, e o processo repete-se até que a 
!    distância entre os pontos convirja para um erro aceitável.
!
! 2. PARA QUE SERVE (APLICAÇÃO GERAL DO MÉTODO):
!    Métodos de otimização baseados em interpolação servem para encontrar picos 
!    (máximos) ou vales (mínimos) de funções objetivo contínuas. São amplamente 
!    aplicados na engenharia para maximizar rendimentos (como a eficiência térmica 
!    de um ciclo, ou a máxima transferência de potência em um circuito elétrico, 
!    como é o caso da função P(x) deste código), minimizar custos operacionais, ou 
!    encontrar dimensões ideais em projetos estruturais. É um método preferido 
!    quando a equação é complexa e calcular sua derivada primeira e segunda 
!    (como exigido no método de Newton) é difícil ou computacionalmente custoso.
!
! 3. INPUTS (ENTRADAS DEFINIDAS NO CÓDIGO):
!    - Parâmetros do sistema (Circuito):
!      V  = 80.0 (Tensão)
!      R1 = 8.0, R2 = 12.0, R3 = 10.0 (Resistências fixas)
!    - Chutes Iniciais de busca para a resistência 'Ra' (x0, x1, x2): 
!      1.0, 8.0 e 20.0
!    - tol : Tolerância aceitável para o critério de paragem (Valor: 1.0e-05)
!
! 4. OUTPUTS (SAÍDAS GERADAS):
!    Os resultados finais são impressos diretamente na tela do terminal:
!    - 'Ra_opt'    : O valor ideal da resistência Ra que maximiza a potência.
!    - 'P_max'     : O valor da potência máxima dissipada encontrada.
!    - 'Iteracoes' : O número total de ciclos necessários para a convergência.
! ==============================================================================

program interpolacao
implicit none

real V, R1, R2, R3, Ra		! Parâmetros do sistema
real P_max, Ra_opt		! Valores de potência máxima e Ra ótimo

integer iter			! Quantidade de iterações
real erro, tol			! Erro e tolerância

! Parâmetros para o Método da Interpolação Quadrática
real x0, x1, x2, x3
real Px0, Px1, Px2, Px3

erro = 100.0
tol = 1.0e-5
iter = 0

V = 80.0		! Tensão [V]
R1 = 8.0		! Resistência R1
R2 = 12.0		! Resistência R2
R3 = 10.0		! Resistência R3

x0 = 1.0
Px0 = P(x0)
x1 = 8.0
Px1 = P(x1)
x2 = 20.0
Px2 = P(x2)
x3 = ponto_max()
Px3 = P(x3)

do while (erro.gt.tol)
	if (x3.gt.x1) then
		if (Px3.ge.Px1) then
			x0 = x1
			x1 = x3
			x3 = ponto_max()
		else
			x2 = x3
			x3 = ponto_max()
		end if
	else
		if (Px3.ge.Px1) then
			x2 = x1
			x1 = x3
			x3 = ponto_max()
		else
			x0 = x3
			x3 = ponto_max()
		end if
	end if
	Px0 = P(x0)
	Px1 = P(x1)
	Px2 = P(x2)
	Px3 = P(x3)
	erro = abs(x3 - x1)
	iter = iter + 1
end do

Ra_opt = x3
P_max = Px3

write(*,'(A,F0.2,A)') 'Ra_opt = ', Ra_opt, ' Ω'
write(*,'(A,F0.2,A)') 'P_max  = ', P_max , ' W'
write(*,*) 'Iteracoes = ', iter

contains
	real function P(x)
		implicit none
		real denom, x
		denom = R1*(x + R2 + R3) + R3*x + R3*R2
		P = ((V*R3*x / denom)**2) / x
	end function P
	
	real function ponto_max()
		implicit none
		real numer, denom
		numer = Px0*(x1**2-x2**2) + Px1*(x2**2-x0**2) + Px2*(x0**2-x1**2)
		denom = 2*(Px0*(x1-x2) + Px1*(x2-x0) + Px2*(x0-x1))
		ponto_max = numer / denom
	end function ponto_max
end program
