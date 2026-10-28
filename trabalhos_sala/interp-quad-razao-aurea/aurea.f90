! ==============================================================================
! PROGRAMA: aurea (Otimização Unidimensional)
!
! 1. RESUMO DO MÉTODO NUMÉRICO:
!    O programa implementa o Método da Razão Áurea (Golden Section Search) 
!    para otimização unidimensional (encontrar o ponto de máximo de uma função).
!    O método começa com um intervalo inicial [xl, xu] e escolhe dois pontos 
!    interiores, x1 e x2, baseados na proporção da razão áurea (~0.618). 
!    Avaliando a função nesses pontos, o algoritmo descarta a sub-região que 
!    não contém o extremo (o máximo, neste caso), reduzindo o intervalo de busca. 
!    A grande vantagem da razão áurea é o reaproveitamento de um dos pontos 
!    calculados na iteração anterior, exigindo apenas uma nova avaliação da 
!    função f(x) por ciclo de repetição.
!
! 2. PARA QUE SERVE (APLICAÇÃO GERAL DO MÉTODO):
!    O Método da Razão Áurea serve para encontrar o máximo ou mínimo (extremos) 
!    de funções unimodais contínuas em um determinado intervalo determinado. Assim 
!    como outros métodos de busca direta, não exige o cálculo analítico de derivadas. 
!    É muito utilizado em engenharia, física e computação para otimização de parâmetros de 
!    projetos, ajuste de curvas, maximização de rendimentos (como a potência elétrica 
!    de um circuito mostrada neste código) ou minimização de custos operacionais e erros. 
!    Destaca-se pelo equilíbrio perfeito entre a robustez geométrica e a eficiência 
!    computacional.
!
! 3. INPUTS (ENTRADAS DEFINIDAS NO CÓDIGO):
!    - Parâmetros do sistema elétrico:
!      V  = 80.0 (Tensão)
!      R1 = 8.0, R2 = 12.0, R3 = 10.0 (Resistências fixas)
!    - Limites Iniciais de busca para a resistência 'Ra': 
!      xl = 0.0 (Limite inferior) e xu = 100.0 (Limite superior)
!    - tol : Tolerância aceitável para o critério de paragem (Valor: 1.0e-5)
!
! 4. OUTPUTS (SAÍDAS GERADAS):
!    Os resultados parciais e finais são impressos na tela do terminal:
!    - Log de iterações detalhando as etapas, pontos (x1, x2), potências e erro de cada ciclo.
!    - 'Ra_opt' : O valor ideal da resistência Ra, que calcula-se pela média dos 
!                 limites finais (xu + xl)/2 para minimizar a incerteza final.
!    - 'P_max'  : O valor analítico da potência máxima dissipada encontrada para esse Ra.
! ==============================================================================

program aurea
implicit none

real V, R1, R2, R3, Ra		! Parâmetros do sistema
real P_max, Ra_opt		! Valores de potência máxima e Ra ótimo

integer iter			! Quantidade de iterações
real erro, tol			! Erro e tolerância

! Parâmetros para o Método da Razão Áurea
real razao, xu, xl, x1, x2, l0, l1, Px1, Px2

erro = 100.0
tol = 1.0e-5
iter = 0

V = 80		! Tensão [V]
R1 = 8		! Resistência R1
R2 = 12		! Resistência R2
R3 = 10		! Resistência R3

razao = (sqrt(5.0) - 1.0)/2.0	! Razão áurea
xu = 100.0
xl = 0.0
l0 = xu - xl
l1 = l0 * razao
x1 = xl + l1
x2 = xu - l1
Px1 = P(x1)
Px2 = P(x2)
write(*,*) '----------------------------------------'

do while(erro.gt.tol)
	write(*,*) 'iteração: ', iter
	write(*,*) 'x1: ', x1, 'Px1: ', Px1
	write(*,*) 'x2: ', x2, 'Px2: ', Px2
	write(*,*) 'erro: ', erro
	write(*,*) '----------------------------------------'
	
	if (Px2.gt.Px1) then		! Se P(x2) > P(x1), o domínio de x1 a xu é eliminado
		xu = x1			! xu vira x1
		x1 = x2			! x1 vira x2
		Px1 = Px2		! P(x1) é igual a P(x2)
		l0 = xu - xl		! l0, l1, x2 e P(x2) são recalculados
		l1 = l0 * razao
		x2 = xu - l1
		Px2 = P(x2)
	else				! Se P(x2) < P(x1), o domínio de x2 a xl é eliminado
		xl = x2			! xl vira x2
		x2 = x1			! x2 vira x1
		Px2 = Px1		! P(x2) é igual a P(x1)
		l0 = xu - xl		! l0, l1, x1 e P(x1) são recalculados
		l1 = l0 * razao
		x1 = xl + l1
		Px1 = P(x1)
	end if
	erro = abs(xu - xl)
	iter = iter + 1
end do

Ra_opt = (xu + xl)/2			! Ra_opt é a média de xu e xl com o erro entre eles muito pequeno
P_max = P(Ra_opt)
write(*,'(A,F0.2,A)') 'Ra_opt: ', Ra_opt, ' Ω'
write(*,'(A,F0.2,A)') 'P_max: ', P_max, ' W'
write(*,*) 'Iteracoes = ', iter

contains
	real function P(x)
	implicit none
	real denom, x
	denom = R1*(x + R2 + R3) + R3*x + R3*R2
	P = ((V*R3*x / denom)**2) / x
	end function P
end program
