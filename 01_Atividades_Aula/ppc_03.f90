program conducao_transiente_thomas
  ! =========================================================================
  ! REGRA DE OURO: Sempre usar implicit none.
  ! =========================================================================
  implicit none

  ! =========================================================================
  ! BLOCO 1: DECLARAÇÃO DE VARIÁVEIS
  ! =========================================================================
  integer, parameter :: N = 7      ! Numero de nos da malha (Exemplo do roteiro)
  integer :: p, max_passos_tempo, i
  
  ! Parametros Fisicos e Matematicos[cite: 130, 132, 133]
  real*8 :: L, alpha, h, k, T_inf, q_ponto, rho, Cp
  real*8 :: dx, dt, Fo, Bi, A_fonte
  
  ! Vetores de Temperaturas (p = atual, p+1 = futuro)[cite: 131, 132]
  real*8, dimension(N) :: T_atual, T_futuro
  
  ! Vetores para o Algoritmo de Thomas (Matriz Tridiagonal)[cite: 125]
  real*8, dimension(N) :: e_vec ! Sub-diagonal inferior (esquerda do no)
  real*8, dimension(N) :: f_vec ! Diagonal principal (o proprio no)
  real*8, dimension(N) :: g_vec ! Super-diagonal superior (direita do no)
  real*8, dimension(N) :: r_vec ! Vetor do lado direito (termos independentes e fonte)

  ! =========================================================================
  ! BLOCO 2: ENTRADA DE DADOS E CONFIGURAÇÃO
  ! Utilizando valores genericos/realistas para teste
  ! =========================================================================
  L = 0.1d0               ! Espessura da parede (m)
  alpha = 1.0d-5          ! Difusividade termica (m^2/s)
  k = 15.0d0              ! Condutividade termica (W/m.K)
  h = 100.0d0             ! Coeficiente de conveccao (W/m^2.K)
  T_inf = 25.0d0          ! Temperatura do fluido externo (C)
  
  ! Geracao interna de energia (W/m^3). 
  ! Mude para 0.0d0 para validar o modelo analitico (Passo 1 do PPC)[cite: 134]
  q_ponto = 50000.0d0     
  rho = 8000.0d0          ! Densidade (kg/m^3)
  Cp = 500.0d0            ! Calor especifico (J/kg.K)
  
  ! Malha Espacial e Temporal[cite: 132]
  dx = L / real(N - 1, 8) 
  dt = 1.0d0              ! Passo de tempo em segundos
  max_passos_tempo = 100  ! Numero total de iteracoes no tempo
  
  ! Calculo dos Parametros Adimensionais[cite: 130, 132]
  Fo = (alpha * dt) / (dx**2)
  Bi = (h * dx) / k       ! Biot referenciado ao dx para o contorno
  A_fonte = (q_ponto / (rho * Cp)) * dt
  
  ! Condicao Inicial: Parede uniformemente a 25 graus
  T_atual = 25.0d0
  
  ! Preparar ficheiro de saida
  open(unit=10, file='resultados_transiente.txt', status='replace')
  write(10, *) 'Resultados para T1 (Centro) e T7 (Borda) ao longo do tempo'
  write(10, *) 'Tempo(s) | T_centro | T_borda'

  print *, 'Iniciando simulacao transiente com Esquema Implicito e Algoritmo de Thomas...'

  ! =========================================================================
  ! BLOCO 3: LOOP TEMPORAL (ESQUEMA IMPLÍCITO)
  ! =========================================================================
  do p = 1, max_passos_tempo
    
    ! -----------------------------------------------------------------------
    ! 3.1: Construcao dos vetores do Sistema Tridiagonal (Esquema Implicito)[cite: 133]
    ! -----------------------------------------------------------------------
    
    ! Noh 1: Simetria Adiabatica (Isolamento no centro)[cite: 133]
    ! Equacao: (1+2Fo)T1 - 2FoT2 = T1_velho + A
    f_vec(1) = 1.0d0 + 2.0d0 * Fo
    g_vec(1) = -2.0d0 * Fo
    r_vec(1) = T_atual(1) + A_fonte
    
    ! Nohs 2 a N-1: Nohs Internos da Parede[cite: 132]
    ! Equacao: -Fo*T_{i-1} + (1+2Fo)T_i - Fo*T_{i+1} = T_i_velho + A
    do i = 2, N-1
      e_vec(i) = -Fo
      f_vec(i) = 1.0d0 + 2.0d0 * Fo
      g_vec(i) = -Fo
      r_vec(i) = T_atual(i) + A_fonte
    end do
    
    ! Noh N (7): Conveccao na Borda Externa[cite: 133]
    ! Equacao: -2Fo*T_{N-1} + (1+2Fo+2FoBi)T_N = T_N_velho + 2FoBiT_inf + Fonte_Borda
    e_vec(N) = -2.0d0 * Fo
    f_vec(N) = 1.0d0 + 2.0d0 * Fo + 2.0d0 * Fo * Bi
    r_vec(N) = T_atual(N) + 2.0d0 * Fo * Bi * T_inf + (Fo * q_ponto * dx**2)/k
    
    ! -----------------------------------------------------------------------
    ! 3.2: Resolucao do Sistema Linear usando o ALGORITMO DE THOMAS (TDMA)
    ! Copiado exatamente do pseudocodigo do documento
    ! -----------------------------------------------------------------------
    
    ! Etapa A: Decomposicao LU (Atualizacao dos vetores f e r progressivamente)[cite: 126]
    do i = 2, N
      e_vec(i) = e_vec(i) / f_vec(i-1)
      f_vec(i) = f_vec(i) - e_vec(i) * g_vec(i-1)
      r_vec(i) = r_vec(i) - e_vec(i) * r_vec(i-1)
    end do
    
    ! Etapa B: Substituicao Regressiva (Calcula o vetor de Temperaturas do Futuro)[cite: 126]
    T_futuro(N) = r_vec(N) / f_vec(N)
    
    do i = N-1, 1, -1
      T_futuro(i) = (r_vec(i) - g_vec(i) * T_futuro(i+1)) / f_vec(i)
    end do
    
    ! -----------------------------------------------------------------------
    ! 3.3: Atualizacao do Tempo e Escrita de Dados
    ! -----------------------------------------------------------------------
    ! O futuro de agora eh o passado de amanha!
    T_atual = T_futuro
    
    ! Escreve os resultados no arquivo a cada iteracao
    write(10, '(F10.2, 2X, F12.6, 2X, F12.6)') real(p)*dt, T_atual(1), T_atual(N)
    
  end do

  close(10)
  print *, 'Simulacao concluida. Resultados salvos em "resultados_transiente.txt"'
  print *, 'Temperatura Final no Centro (Noh 1): ', T_atual(1), ' C'
  print *, 'Temperatura Final na Borda (Noh 7): ', T_atual(N), ' C'

end program conducao_transiente_thomas