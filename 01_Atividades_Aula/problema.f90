program aleta_transferencia_calor_completo
  ! =========================================================================
  ! REGRA DE OURO: Sempre usar implicit none.
  ! PROGRAMA COMPLETO INTEGRADO PARA PROVA DE CÁLCULO NUMÉRICO APLICADO
  ! Tópicos: Diferenças Finitas, Decomposição LU, Zeros de Função e Otimização
  ! =========================================================================
  implicit none

  ! =========================================================================
  ! BLOCO 1: DECLARAÇÃO DE VARIÁVEIS FISICAS E NUMÉRICAS
  ! =========================================================================
  integer, parameter :: N_nos = 6 ! Numero de nos da discretizacao da aleta
  
  ! Propriedades Fisicas da Aleta
  real*8 :: t_esp, W_larg, k_cond, T_base, T_inf, L_aleta, h_conv
  
  ! Variaveis do Perfil de Temperatura
  real*8, dimension(N_nos) :: T_perfil
  
  ! Variaveis para o Zero de Funcao (Busca de h)
  real*8 :: h_encontrado, T_alvo_ponta
  
  ! Variaveis para Otimizacao (Busca de L otimo)
  real*8 :: L_otimo, custo_minimo
  
  integer :: i

  ! =========================================================================
  ! BLOCO 2: INICIALIZAÇÃO DOS PARÂMETROS DO PROBLEMA
  ! =========================================================================
  t_esp = 0.005d0       ! Espessura da aleta (5 mm = 0.005 m)
  W_larg = 1.000d0      ! Largura unitaria (1 m)
  k_cond = 180.0d0      ! Condutividade termica do Aluminio (W/m.K)
  T_base = 100.0d0      ! Temperatura na base (Celsius)
  T_inf = 25.0d0        ! Temperatura do ar ambiente (Celsius)
  
  L_aleta = 0.15d0      ! Comprimento inicial da aleta (0.15 m)
  h_conv = 40.0d0       ! Coeficiente de conveccao inicial (W/m2.K)

  print *, '================================================================='
  print *, '   CÁLCULO NUMÉRICO APLICADO - PROBLEMA DA ALETA DE CALOR        '
  print *, '================================================================='
  print *, ''

  ! =========================================================================
  ! ETAPA 1: SOLUÇÃO DO SISTEMA LINEAR VIA DECOMPOSIÇÃO LU
  ! =========================================================================
  print *, '--- PARTE 1: Perfil de Temperaturas via Decomposicao LU ---'
  
  call resolve_aleta_lu(L_aleta, h_conv, k_cond, t_esp, T_base, T_inf, N_nos, T_perfil)
  
  print *, 'Perfil de Temperaturas ao longo da aleta (N = 6 nos):'
  do i = 1, N_nos
    print '(A, I2, A, F8.4, A)', 'No ', i, ': T = ', T_perfil(i), ' °C'
  end do
  print *, ''

  ! =========================================================================
  ! ETAPA 2: ZERO DE FUNÇÃO (METODO DA BISSEÇÃO PARA ENCONTRAR h)
  ! Objetivo: Encontrar h tal que T(ponta) = 50 °C
  ! =========================================================================
  print *, '--- PARTE 2: Zero de Funcao - Encontrar h para T(L) = 50 °C ---'
  T_alvo_ponta = 50.0d0
  
  call busca_h_bissecao(T_alvo_ponta, L_aleta, k_cond, t_esp, T_base, T_inf, N_nos, h_encontrado)
  
  print *, 'Coeficiente de Conveccao h necessario: ', h_encontrado, ' W/(m2.K)'
  print *, ''

  ! =========================================================================
  ! ETAPA 3: OTIMIZAÇÃO (MÉTODO DA RAZÃO ÁUREA PARA MINIMIZAR CUSTO)
  ! Objetivo: Encontrar o comprimento L otimo no intervalo [0.05, 0.30] m
  ! =========================================================================
  print *, '--- PARTE 3: Otimizacao da Aleta via Razao Aurea ---'
  
  call otimiza_comprimento_razao_aurea(h_conv, k_cond, t_esp, T_base, T_inf, N_nos, L_otimo, custo_minimo)
  
  print *, 'Comprimento Otimo da Aleta (L) = ', L_otimo, ' m'
  print *, 'Custo Minimo Total            = R$ ', custo_minimo
  print *, '================================================================='

contains

  ! =========================================================================
  ! SUBROTINA A: MONTAGEM DO SISTEMA E RESOLUÇÃO POR DECOMPOSIÇÃO LU
  ! =========================================================================
  subroutine resolve_aleta_lu(L_val, h_val, k_val, t_val, Tb, Tinf, n, T_sol)
    implicit none
    integer, intent(in) :: n
    real*8, intent(in) :: L_val, h_val, k_val, t_val, Tb, Tinf
    real*8, intent(out), dimension(n) :: T_sol
    
    ! Matrizes do Sistema A*T = b e Matrizes L e U
    real*8, dimension(n,n) :: A_mat, L_mat, U_mat
    real*8, dimension(n) :: b_vec, y_vec
    
    real*8 :: dx, m2, soma
    integer :: lin, col, k_pivo, m_idx
    
    ! 1. Discretizacao espacial
    dx = L_val / real(n - 1, 8)
    m2 = (2.0d0 * h_val) / (k_val * t_val) ! Parametro m^2 da aleta
    
    ! 2. Construcao da Matriz A e Vetor b (Zera tudo primeiro)
    A_mat = 0.0d0
    b_vec = 0.0d0
    
    ! No 1: Base com temperatura conhecida T(0) = Tb
    A_mat(1,1) = 1.0d0
    b_vec(1) = Tb
    
    ! Nos Internos (2 ate n-1): Diferencas Finitas Centrais
    ! T_{i-1} - (2 + m^2 * dx^2)*T_i + T_{i+1} = -m^2 * dx^2 * T_inf
    do lin = 2, n - 1
      A_mat(lin, lin - 1) = 1.0d0
      A_mat(lin, lin)     = -(2.0d0 + m2 * dx**2)
      A_mat(lin, lin + 1) = 1.0d0
      b_vec(lin)          = -m2 * dx**2 * Tinf
    end do
    
    ! No n (Ponta da Aleta): Conveccao na extremidade
    ! 2*T_{n-1} - (2 + 2*h*dx/k + m^2*dx^2)*T_n = -(2*h*dx/k + m^2*dx^2)*T_inf
    A_mat(n, n - 1) = 2.0d0
    A_mat(n, n)     = -(2.0d0 + (2.0d0 * h_val * dx / k_val) + m2 * dx**2)
    b_vec(n)        = -((2.0d0 * h_val * dx / k_val) + m2 * dx**2) * Tinf
    
    ! 3. DECOMPOSIÇÃO LU (Algoritmo de Doolittle: L tem 1s na diagonal)
    L_mat = 0.0d0
    U_mat = 0.0d0
    
    do k_pivo = 1, n
      L_mat(k_pivo, k_pivo) = 1.0d0
      
      ! Preenche linha de U
      do col = k_pivo, n
        soma = 0.0d0
        do m_idx = 1, k_pivo - 1
          soma = soma + L_mat(k_pivo, m_idx) * U_mat(m_idx, col)
        end do
        U_mat(k_pivo, col) = A_mat(k_pivo, col) - soma
      end do
      
      ! Preenche coluna de L
      do lin = k_pivo + 1, n
        soma = 0.0d0
        do m_idx = 1, k_pivo - 1
          soma = soma + L_mat(lin, m_idx) * U_mat(m_idx, k_pivo)
        end do
        L_mat(lin, k_pivo) = (A_mat(lin, k_pivo) - soma) / U_mat(k_pivo, k_pivo)
      end do
    end do
    
    ! 4. RESOLUÇÃO DOS SISTEMAS TRIANGULARES
    ! Passo A: Substituicao Progressiva (L * y = b)
    do lin = 1, n
      soma = b_vec(lin)
      do col = 1, lin - 1
        soma = soma - L_mat(lin, col) * y_vec(col)
      end do
      y_vec(lin) = soma / L_mat(lin, lin)
    end do
    
    ! Passo B: Substituicao Regressiva (U * T = y)
    T_sol(n) = y_vec(n) / U_mat(n, n)
    do lin = n - 1, 1, -1
      soma = y_vec(lin)
      do col = lin + 1, n
        soma = soma - U_mat(lin, col) * T_sol(col)
      end do
      T_sol(lin) = soma / U_mat(lin, lin)
    end do
    
  end subroutine resolve_aleta_lu

  ! =========================================================================
  ! SUBROTINA B: BUSCA DO ZERO DE FUNÇÃO (MÉTODO DA BISSEÇÃO)
  ! Encontra h tal que f(h) = T_ponta(h) - T_alvo = 0
  ! =========================================================================
  subroutine busca_h_bissecao(T_target, L_val, k_val, t_val, Tb, Tinf, n, h_sol)
    implicit none
    real*8, intent(in) :: T_target, L_val, k_val, t_val, Tb, Tinf
    integer, intent(in) :: n
    real*8, intent(out) :: h_sol
    
    real*8 :: h_inf, h_sup, h_mid, f_inf, f_mid, tol
    real*8, dimension(n) :: T_temp
    integer :: iter, max_iter
    
    ! Intervalo de busca inicial para h [W/m2.K]
    h_inf = 5.0d0
    h_sup = 300.0d0
    tol = 1.0d-5
    max_iter = 100
    
    do iter = 1, max_iter
      h_mid = (h_inf + h_sup) / 2.0d0
      
      ! Avalia a funcao no ponto inferior: f(h) = T_ponta(h) - T_target
      call resolve_aleta_lu(L_val, h_inf, k_val, t_val, Tb, Tinf, n, T_temp)
      f_inf = T_temp(n) - T_target
      
      ! Avalia a funcao no ponto medio
      call resolve_aleta_lu(L_val, h_mid, k_val, t_val, Tb, Tinf, n, T_temp)
      f_mid = T_temp(n) - T_target
      
      ! Criterio de parada
      if (abs(f_mid) < tol .or. (h_sup - h_inf)/2.0d0 < tol) exit
      
      ! Atualizacao do intervalo da Bissecao
      if (f_inf * f_mid < 0.0d0) then
        h_sup = h_mid
      else
        h_inf = h_mid
      end if
    end do
    
    h_sol = h_mid
  end subroutine busca_h_bissecao

  ! =========================================================================
  ! SUBROTINA C: OTIMIZAÇÃO POR RAZÃO ÁUREA
  ! Minimiza o Custo C(L) = Custo_Material + Custo_Perda_Termica
  ! =========================================================================
  subroutine otimiza_comprimento_razao_aurea(h_val, k_val, t_val, Tb, Tinf, n, L_sol, custo_sol)
    implicit none
    real*8, intent(in) :: h_val, k_val, t_val, Tb, Tinf
    integer, intent(in) :: n
    real*8, intent(out) :: L_sol, custo_sol
    
    real*8 :: a_gold, b_gold, tau, x1, x2, f1, f2, tol
    integer :: iter, max_iter
    
    ! Intervalo de busca para o comprimento L [m]
    a_gold = 0.05d0
    b_gold = 0.30d0
    tol = 1.0d-5
    max_iter = 100
    
    ! Proporcao Aurea
    tau = (sqrt(5.0d0) - 1.0d0) / 2.0d0
    
    ! Pontos internos
    x1 = b_gold - tau * (b_gold - a_gold)
    x2 = a_gold + tau * (b_gold - a_gold)
    
    f1 = avalia_custo_aleta(x1, h_val, k_val, t_val, Tb, Tinf, n)
    f2 = avalia_custo_aleta(x2, h_val, k_val, t_val, Tb, Tinf, n)
    
    do iter = 1, max_iter
      if (abs(b_gold - a_gold) < tol) exit
      
      if (f1 > f2) then
        a_gold = x1
        x1 = x2
        f1 = f2
        x2 = a_gold + tau * (b_gold - a_gold)
        f2 = avalia_custo_aleta(x2, h_val, k_val, t_val, Tb, Tinf, n)
      else
        b_gold = x2
        x2 = x1
        f2 = f1
        x1 = b_gold - tau * (b_gold - a_gold)
        f1 = avalia_custo_aleta(x1, h_val, k_val, t_val, Tb, Tinf, n)
      end if
    end do
    
    L_sol = (a_gold + b_gold) / 2.0d0
    custo_sol = avalia_custo_aleta(L_sol, h_val, k_val, t_val, Tb, Tinf, n)
  end subroutine otimiza_comprimento_razao_aurea

  ! =========================================================================
  ! FUNÇÃO AUXILIAR: CALCULO DA FUNÇÃO OBJETIVO DE CUSTO
  ! =========================================================================
  real*8 function avalia_custo_aleta(L_test, h_val, k_val, t_val, Tb, Tinf, n)
    implicit none
    real*8, intent(in) :: L_test, h_val, k_val, t_val, Tb, Tinf
    integer, intent(in) :: n
    
    real*8, dimension(n) :: T_temp
    real*8 :: q_base, custo_mat, custo_termico
    
    ! 1. Calcula o perfil de temperatura para este L_test
    call resolve_aleta_lu(L_test, h_val, k_val, t_val, Tb, Tinf, n, T_temp)
    
    ! 2. Calor retirado na base q_base = -k * A * (T2 - T1)/dx
    q_base = -k_val * (1.0d0 * t_val) * (T_temp(2) - T_temp(1)) / (L_test / real(n-1, 8))
    
    ! 3. Modelo de Custo: Custo do Material cresce com L; Custo Termico cai com maior calor retirado
    custo_mat = 150.0d0 * L_test                   ! R$ 150 por metro de aleta
    custo_termico = 5000.0d0 / (abs(q_base) + 1.0d0) ! Penalidade por pouca troca termica
    
    avalia_custo_aleta = custo_mat + custo_termico
  end function avalia_custo_aleta

end program aleta_transferencia_calor_completo
