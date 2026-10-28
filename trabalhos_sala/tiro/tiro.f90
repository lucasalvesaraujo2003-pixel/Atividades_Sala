! ==============================================================================
! RESUMO DO MÉTODO:
! Este algoritmo implementa o Método do Tiro (Shooting Method) para solucionar 
! um Problema de Valor de Contorno (PVC). O código reduz uma Equação Diferencial 
! Ordinária (EDO) de segunda ordem a um sistema de duas EDOs de primeira ordem, 
! que são integradas espacialmente utilizando o Método de Euler. O programa 
! "atira" trajetórias a partir do contorno inicial usando palpites para o 
! gradiente de temperatura (f_old e f_new). O erro no contorno final é avaliado 
! e o palpite inicial é corrigido iterativamente através do Método da Secante, 
! até que o erro (erro_new) seja inferior à tolerância definida (tol).
!
! APLICAÇÃO GERAL:
! Resolução numérica de problemas unidimensionais de valor de contorno. 
! Neste cenário específico, o código simula a condução térmica unidimensional 
! em regime permanente ao longo de uma barra sujeita a uma taxa de geração 
! volumétrica de calor uniforme. Na engenharia, esse modelo é utilizado para 
! analisar a distribuição de temperatura em condutores elétricos (aquecimento 
! por efeito Joule), elementos combustíveis nucleares ou paredes de reatores.
!
! INPUTS (Entradas):
! - Parâmetros da malha geométrica: n, L e delta_x.
! - Propriedades térmicas: k (condutividade) e q_dot (geração de calor).
! - Condições de contorno de temperatura (Dirichlet): T_0 e T_L.
! - Parâmetros numéricos: tol e os palpites iniciais embutidos na lógica 
!   do código (f_old = -20.0d0 e f_new = 20.0d0).
!
! OUTPUTS (Saídas):
! - Impressão direta no terminal contendo o resultado da convergência:
!   f_end    : O gradiente de temperatura inicial convergido.
!   T_end    : A temperatura final calculada (que deve bater com T_L).
!   erro_end : O erro residual na fronteira final.
!   iter     : O número total de iterações até a convergência.
! ==============================================================================

program tiro
    implicit none
    integer, parameter :: n = 100                            ! Número de nós
    real(8), parameter :: tol = 1.0d-6                       ! Tolerância para o erro
    real(8) :: L, q_dot, k, delta_x, pos_x(n)
    real(8) :: T_end, f_end, erro_end
    real(8) :: f_next, f_atual
    real(8) :: T_old(n), f_old, erro_old
    real(8) :: T_new(n), f_new, erro_new
    real(8), parameter :: T_0 = 200.0d0, T_L = 1500.0d0      ! Condições de contorno da temperatura [°C]
    integer :: i, iter

    L = 2.0d0               ! Comprimento da barra [m]
    k = 60.0d0              ! Condutividade térmica [W/m.K]
    q_dot = 1.0d6           ! Geração volumétrica de calor [W/m^3]
    delta_x = L/(n-1)       ! Espaçamento uniforme da malha

    ! Vetor de posições x, indo de 0 até L
    pos_x(1) = 0.0d0
    do i=2, n
        pos_x(i) = pos_x(i-1) + delta_x
    end do

    ! Método de Euler para resolução de EDO de ordem 1
    ! Avaliando o primeiro chute
    T_old(1) = T_0
    f_old = -20.0d0
    f_atual = f_old
    do i=2, n
        T_old(i) = T_old(i-1) + f_atual*delta_x
        f_atual = f_atual - (q_dot/k)*delta_x
    end do
    erro_old = T_old(n)-T_L

    ! Avaliando o segundo chute
    T_new(1) = T_0
    f_new = 20.0d0
    f_atual = f_new
    do i=2, n
        T_new(i) = T_new(i-1) + f_atual*delta_x
        f_atual = f_atual - (q_dot/k)*delta_x
    end do
    erro_new = T_new(n)-T_L

    iter = 1
    do while (abs(erro_new).gt.tol)
        f_next = f_new - erro_new*((f_new-f_old)/(erro_new-erro_old))
        f_old = f_new
        erro_old = erro_new
        f_new = f_next
        T_new(1) = T_0
        f_atual = f_new
        do i=2, n
            T_new(i) = T_new(i-1) + f_atual*delta_x
            f_atual = f_atual - (q_dot/k)*delta_x
        end do
        erro_new = T_new(n)-T_L
        iter = iter + 1
        f_end = f_new
        T_end = T_new(n)
        erro_end = erro_new
    end do

    write(*, '(A, F10.4, A, F10.2, A, E12.4, A, I0)') &
      "Derivada f(0) = ", f_end, &
      " | Temp Final T(L) = ", T_end, &
      " | Erro = ", erro_end, &
      " | Iteracoes = ", iter

end program tiro
