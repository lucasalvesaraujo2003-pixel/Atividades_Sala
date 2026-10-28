program otimizacao_isolamento
  ! =========================================================================
  ! REGRA DE OURO: Sempre usar implicit none.
  ! =========================================================================
  implicit none

  ! =========================================================================
  ! BLOCO 1: DECLARAÇÃO DE VARIÁVEIS GERAIS
  ! =========================================================================
  ! Variaveis para a Funcao Objetivo e Custos[cite: 97]
  real*8 :: r1, L, Ti, Tinf, k, h, t, pE, cI, pi
  real*8 :: e_val, Custo_E, Custo_I, Custo_Total
  
  ! Variaveis para o Metodo da Razao Aurea[cite: 98]
  real*8 :: a, b, tol, tau, x1, x2, f1, f2, e_min
  integer :: iteracao, max_iter, avaliacoes_funcao

  ! =========================================================================
  ! BLOCO 2: INICIALIZAÇÃO E PARÂMETROS
  ! =========================================================================
  pi = acos(-1.0d0)
  
  ! Parametros termicos e economicos do enunciado[cite: 97]
  r1 = 0.040d0      ! m
  L = 20.0d0        ! m
  Ti = 180.0d0      ! C
  Tinf = 40.0d0     ! C
  k = 0.045d0       ! W/(m K)
  h = 12.0d0        ! W/(m2 K)
  t = 30000.0d0     ! h
  pE = 0.70d0       ! R$/kWh
  cI = 9000.0d0     ! R$/m3
  
  ! Parametros do Metodo Numerico (Razao Aurea)[cite: 98]
  a = 0.0d0         ! Limite inferior do intervalo (m)
  b = 0.12d0        ! Limite superior do intervalo (m)
  tol = 1.0d-6      ! Tolerancia do intervalo (m)
  max_iter = 200
  iteracao = 0
  avaliacoes_funcao = 0

  print *, 'Iniciando Minimizacao pelo Metodo da Razao Aurea...'

  ! =========================================================================
  ! BLOCO 3: ALGORITMO DA RAZÃO ÁUREA
  ! =========================================================================
  ! Proporcao aurea[cite: 98]
  tau = (sqrt(5.0d0) - 1.0d0) / 2.0d0
  
  ! Calculo inicial dos pontos internos[cite: 98]
  x1 = b - tau * (b - a)
  x2 = a + tau * (b - a)
  
  ! Avalia a funcao custo nos dois pontos
  f1 = calcula_custo(x1, r1, L, Ti, Tinf, k, h, t, pE, cI, pi)
  f2 = calcula_custo(x2, r1, L, Ti, Tinf, k, h, t, pE, cI, pi)
  avaliacoes_funcao = avaliacoes_funcao + 2

  do while (abs(b - a) > tol .and. iteracao < max_iter)
    iteracao = iteracao + 1
    
    ! Logica de Reducao do Intervalo[cite: 98]
    if (f1 > f2) then
      ! O minimo deve estar mais a direita, o limite 'a' passa a ser 'x1'
      a = x1
      x1 = x2
      f1 = f2
      
      ! Calcula apenas um novo ponto
      x2 = a + tau * (b - a)
      f2 = calcula_custo(x2, r1, L, Ti, Tinf, k, h, t, pE, cI, pi)
      avaliacoes_funcao = avaliacoes_funcao + 1
    else
      ! O minimo deve estar mais a esquerda, o limite 'b' passa a ser 'x2'
      b = x2
      x2 = x1
      f2 = f1
      
      ! Calcula apenas um novo ponto
      x1 = b - tau * (b - a)
      f1 = calcula_custo(x1, r1, L, Ti, Tinf, k, h, t, pE, cI, pi)
      avaliacoes_funcao = avaliacoes_funcao + 1
    end if
  end do

  ! A espessura otima e o meio do intervalo final[cite: 98]
  e_min = (a + b) / 2.0d0
  Custo_Total = calcula_custo(e_min, r1, L, Ti, Tinf, k, h, t, pE, cI, pi)

  ! =========================================================================
  ! BLOCO 4: IMPRESSÃO DOS RESULTADOS FINAIS
  ! =========================================================================
  print *, '---------------------------------------------------'
  print *, 'Convergencia atingida!'
  print *, 'Espessura Otima (e) = ', e_min, ' m'
  print *, 'Custo Minimo Total  = R$ ', Custo_Total
  print *, 'Iteracoes           = ', iteracao
  print *, 'Avaliacoes de f(x)  = ', avaliacoes_funcao

contains

  ! =========================================================================
  ! FUNÇÃO AUXILIAR: CÁLCULO DA FUNÇÃO OBJETIVO C(e)
  ! COPIE SE A QUESTÃO FOR SOBRE AVALIAÇÃO ECONÓMICA / TRANSFERÊNCIA CALOR
  ! =========================================================================
  real*8 function calcula_custo(e, r1_f, L_f, Ti_f, Tinf_f, k_f, h_f, t_f, pE_f, cI_f, pi_f)
    implicit none
    real*8, intent(in) :: e, r1_f, L_f, Ti_f, Tinf_f, k_f, h_f, t_f, pE_f, cI_f, pi_f
    
    real*8 :: R_cond, R_conv, q_linha, C_energia, C_isolante
    
    ! 1. Resistencias termicas (Conducao + Conveccao)[cite: 97]
    R_cond = log((r1_f + e) / r1_f) / (2.0d0 * pi_f * k_f)
    R_conv = 1.0d0 / (2.0d0 * pi_f * h_f * (r1_f + e))
    
    ! 2. Taxa de calor por comprimento (q')[cite: 97]
    q_linha = (Ti_f - Tinf_f) / (R_cond + R_conv)
    
    ! 3. Custo de Energia Acumulado[cite: 97]
    C_energia = pE_f * (L_f * q_linha * t_f) / 1000.0d0
    
    ! 4. Custo Instalado do Isolante[cite: 97]
    C_isolante = cI_f * pi_f * L_f * ((r1_f + e)**2 - r1_f**2)
    
    ! 5. Funcao Objetivo Total (Minimizar)[cite: 97]
    calcula_custo = C_energia + C_isolante
  end function calcula_custo

end program otimizacao_isolamento