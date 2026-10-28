program sedimentacao_esfera_rk4
  ! =========================================================================
  ! REGRA DE OURO: Sempre usar implicit none.
  ! =========================================================================
  implicit none

  ! =========================================================================
  ! BLOCO 1: DECLARAÇÃO DE VARIÁVEIS
  ! =========================================================================
  ! Parametros fisicos adimensionais[cite: 103, 108]
  real*8 :: St, Re_s 
  
  ! Variaveis de tempo e avanco
  real*8 :: t, t_final, dt
  integer :: num_passos, i
  
  ! Variaveis da velocidade (Solucao Numerica e Analitica)[cite: 103, 105, 108]
  real*8 :: v_num, v_analitico
  
  ! Coeficientes do Metodo de Runge-Kutta de 4a Ordem
  real*8 :: k1, k2, k3, k4
  
  ! Flag de controle para escolher qual modelo resolver (1 = Stokes, 2 = Oseen)
  integer :: modelo

  ! =========================================================================
  ! BLOCO 2: ENTRADA DE DADOS E CONFIGURAÇÃO INICIAL
  ! =========================================================================
  ! Escolha o modelo a simular[cite: 108]:
  ! 1 -> Limite de Stokes (Re -> 0)
  ! 2 -> Efeito Inercial de Oseen (Re /= 0)
  modelo = 2
  
  ! Valores Fisicos de Teste[cite: 103, 108]
  St = 1.0d0
  Re_s = 0.5d0 ! Usado apenas se modelo = 2
  
  ! Configuracao do refinamento temporal[cite: 108]
  t = 0.0d0
  t_final = 5.0d0  ! Tempo total de simulacao (adimensional)
  dt = 0.01d0      ! Passo de integracao (h)
  num_passos = int(t_final / dt)
  
  ! Condicao Inicial (A esfera parte do repouso)[cite: 103, 107]
  v_num = 0.0d0

  ! =========================================================================
  ! BLOCO 3: INICIALIZAÇÃO DE ARQUIVO DE SAÍDA (FILE I/O)
  ! Essencial para plotar os graficos exigidos na atividade[cite: 108]
  ! =========================================================================
  open(unit=10, file='resultados_rk4.txt', status='replace')
  write(10, *) 'Tempo_Adimensional | Vel_Numerica | Vel_Analitica'

  print *, 'Iniciando integracao por Runge-Kutta de 4a ordem...'
  print *, 'Tempo   |  V_Num   | V_Exata'

  ! =========================================================================
  ! BLOCO 4: LOOP DE INTEGRAÇÃO TEMPORAL (AVANÇO NO TEMPO)[cite: 107, 108]
  ! =========================================================================
  do i = 1, num_passos
    
    ! -----------------------------------------------------------------------
    ! ETAPA A: Solucao Analitica (Para efeitos de validacao)[cite: 106, 108]
    ! -----------------------------------------------------------------------
    if (modelo == 1) then
       ! Solucao Exata de Stokes (Re -> 0)[cite: 103]
       v_analitico = 1.0d0 - exp(-t / St)
    else
       ! Solucao Exata de Oseen (Re /= 0) via Equacao de Ricatti[cite: 105]
       v_analitico = solucao_analitica_oseen(t, St, Re_s)
    end if
    
    ! Escreve os dados atuais no terminal e no arquivo de saida
    ! Imprime a cada 10 iteracoes para nao poluir a tela
    if (mod(i, 10) == 0 .or. i == 1) then
       print '(F8.4, 2X, F8.4, 2X, F8.4)', t, v_num, v_analitico
    end if
    write(10, '(F12.6, 2X, F12.6, 2X, F12.6)') t, v_num, v_analitico

    ! -----------------------------------------------------------------------
    ! ETAPA B: Metodo Numerico - Calculo dos Coeficientes de Runge-Kutta[cite: 107]
    ! -----------------------------------------------------------------------
    ! k1 avalia a funcao na posicao original[cite: 107]
    k1 = func_edo(t, v_num, modelo, St, Re_s)
    
    ! k2 avalia a funcao no meio do passo de tempo, usando k1[cite: 107]
    k2 = func_edo(t + 0.5d0*dt, v_num + 0.5d0*k1*dt, modelo, St, Re_s)
    
    ! k3 avalia a funcao no meio do passo de tempo, usando k2[cite: 107]
    k3 = func_edo(t + 0.5d0*dt, v_num + 0.5d0*k2*dt, modelo, St, Re_s)
    
    ! k4 avalia a funcao no final do passo de tempo, usando k3[cite: 107]
    k4 = func_edo(t + dt, v_num + k3*dt, modelo, St, Re_s)
    
    ! -----------------------------------------------------------------------
    ! ETAPA C: Avanço do Tempo e da Variável[cite: 107]
    ! -----------------------------------------------------------------------
    ! Previsao da proxima velocidade usando a media ponderada dos 'k'[cite: 107]
    v_num = v_num + (dt / 6.0d0) * (k1 + 2.0d0*k2 + 2.0d0*k3 + k4)
    
    ! Avanca o relogio
    t = t + dt
    
  end do
  
  close(10)
  print *, 'Simulacao concluida. Resultados gravados em "resultados_rk4.txt".'

contains

  ! =========================================================================
  ! FUNÇÃO 1: EQUAÇÃO DIFERENCIAL ORDINÁRIA (EDO)
  ! Contem a funcao f(t, v) = dv/dt[cite: 107]
  ! =========================================================================
  real*8 function func_edo(t_val, v_val, mod_val, St_val, Re_val)
    implicit none
    real*8, intent(in) :: t_val, v_val, St_val, Re_val
    integer, intent(in) :: mod_val
    
    if (mod_val == 1) then
       ! Modelo Linear (Stokes)[cite: 103]: St * dv/dt = 1 - v -> dv/dt = (1 - v) / St
       func_edo = (1.0d0 - v_val) / St_val
    else
       ! Modelo Nao-Linear (Oseen)[cite: 104]: St * dv/dt = 1 - v - (3/8)*Re*v^2
       ! Logo: dv/dt = (1 - v - (3/8)*Re*v^2) / St
       func_edo = (1.0d0 - v_val - (3.0d0/8.0d0) * Re_val * (v_val**2)) / St_val
    end if
  end function func_edo

  ! =========================================================================
  ! FUNÇÃO 2: SOLUÇÃO ANALÍTICA PARA O MODELO DE OSEEN (EQUAÇÃO DE RICATTI)
  ! Calculada estritamente conforme a formula (10) e (11) do documento[cite: 105]
  ! =========================================================================
  real*8 function solucao_analitica_oseen(t_val, St_val, Re_val)
    implicit none
    real*8, intent(in) :: t_val, St_val, Re_val
    real*8 :: v_estrela, Q, P
    
    ! v* (Velocidade assintotica limite)[cite: 105]
    v_estrela = (-1.0d0 + sqrt(1.0d0 + 1.5d0 * Re_val)) / (0.75d0 * Re_val)
    
    ! Parametros Q e P[cite: 105]
    Q = (3.0d0 / (8.0d0 * St_val)) * Re_val
    P = -((0.75d0 * Re_val) / St_val) * v_estrela - (1.0d0 / St_val)
    
    ! Formula final (10)[cite: 105]
    solucao_analitica_oseen = v_estrela + &
         1.0d0 / ( (Q / P) + ((-1.0d0 / v_estrela) - (Q / P)) * exp(-P * t_val) )
  end function solucao_analitica_oseen

end program sedimentacao_esfera_rk4