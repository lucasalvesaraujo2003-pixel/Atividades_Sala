program newton_raphson_sistemas_tubulacoes
  ! =========================================================================
  ! REGRA DE OURO: Sempre usar implicit none.
  ! =========================================================================
  implicit none
  
  ! =========================================================================
  ! BLOCO 1: DECLARAÇÃO DE VARIÁVEIS FÍSICAS E MATEMÁTICAS
  ! =========================================================================
  ! Variáveis geométricas e físicas do problema[cite: 57]
  real*8 :: Q_total, g, pi
  real*8 :: L1, D1, f1, K1
  real*8 :: L2, D2, f2, K2, Kv
  real*8 :: R1, R2, U1, U2, hL1, hL2
  
  ! Variáveis do método de Newton-Raphson[cite: 59, 60]
  real*8 :: Q1, Q2             ! Nossas incógnitas (Vazões nos ramos 1 e 2)
  real*8 :: F1, F2             ! Funções do sistema (Resíduos que queremos zerar)
  real*8 :: J11, J12, J21, J22 ! Elementos da Matriz Jacobiana (Derivadas parciais)
  real*8 :: delta_Q1, delta_Q2 ! Correções de vazão a cada iteração
  
  ! Variáveis para a Regra de Cramer e Controle[cite: 59, 60]
  real*8 :: a, b, c, d, r1_cramer, r2_cramer, Det ! Para resolver o sistema linear 2x2
  real*8 :: erro, tolerancia
  integer :: iteracao, max_iter

  ! =========================================================================
  ! BLOCO 2: ENTRADA DE DADOS E CONFIGURAÇÃO
  ! COPIE SE PRECISAR DECLARAR CONSTANTES FÍSICAS (COMO PI E GRAVIDADE).
  ! =========================================================================
  pi = acos(-1.0d0)            ! Melhor forma de declarar o Pi com precisão máxima no Fortran
  g = 9.81d0
  Q_total = 0.050d0

  ! Dados do ramo 1[cite: 57]
  L1 = 30.0d0; D1 = 0.10d0; f1 = 0.020d0; K1 = 1.5d0
  ! Dados do ramo 2[cite: 57]
  L2 = 20.0d0; D2 = 0.08d0; f2 = 0.022d0; K2 = 1.0d0; Kv = 8.0d0

  ! Cálculo das constantes de resistência (R1 e R2) conforme deduzido no papel[cite: 58]
  R1 = (f1 * L1/D1 + K1) * (8.0d0 / (g * pi**2 * D1**4))
  R2 = (f2 * L2/D2 + K2 + Kv) * (8.0d0 / (g * pi**2 * D2**4))

  ! Estimativas iniciais fornecidas pelo problema[cite: 60]
  Q1 = 0.030d0
  Q2 = 0.020d0

  ! Configuração do critério de parada (loop numérico)[cite: 60]
  tolerancia = 1.0d-8
  max_iter = 100
  erro = 1.0d0
  iteracao = 0

  print *, 'Iniciando Newton-Raphson para Sistemas Nao Lineares...'
  print *, 'Iteracao |       Q1       |       Q2       |      Erro      '
  print *, '--------------------------------------------------------------'

  ! =========================================================================
  ! BLOCO 3: LOOP PRINCIPAL DO MÉTODO DE NEWTON-RAPHSON PARA SISTEMAS
  ! COPIE ESTE BLOCO SE A QUESTÃO PEDIR PARA RESOLVER SISTEMAS NÃO LINEARES.
  ! =========================================================================
  do while (erro > tolerancia .and. iteracao < max_iter)
    iteracao = iteracao + 1
    
    ! 3.1: Avaliar as funções do sistema (F1 e F2)[cite: 58, 60]
    ! F1 = Q1 + Q2 - Q_total = 0 (Conservacao de massa)
    ! F2 = R1*Q1^2 - R2*Q2^2 = 0 (Igualdade da perda de carga)
    F1 = Q1 + Q2 - Q_total
    F2 = R1 * Q1**2 - R2 * Q2**2
    
    ! 3.2: Calcular a Matriz Jacobiana (Derivadas parciais)[cite: 59]
    ! Linha 1: Derivadas de F1 em relação a Q1 e Q2
    J11 = 1.0d0
    J12 = 1.0d0
    ! Linha 2: Derivadas de F2 em relação a Q1 e Q2
    J21 = 2.0d0 * R1 * Q1
    J22 = -2.0d0 * R2 * Q2
    
    ! =======================================================================
    ! BLOCO 4: REGRA DE CRAMER (RESOLUÇÃO DO SISTEMA 2x2)
    ! COPIE APENAS SE PRECISAR RESOLVER UM SISTEMA LINEAR DE 2 VARIÁVEIS[cite: 59].
    ! J * [deltaQ] = [-F]
    ! =======================================================================
    a = J11; b = J12
    c = J21; d = J22
    r1_cramer = -F1
    r2_cramer = -F2
    
    ! Determinante principal (D = ad - bc)
    Det = (a * d) - (b * c)
    
    ! Proteção contra divisão por zero (Matriz singular)
    if (abs(Det) < 1.0d-12) then
      print *, 'ERRO: Determinante nulo. Jacobiana singular!'
      stop
    end if
    
    ! Cálculo das correções (deltas)[cite: 59]
    delta_Q1 = (r1_cramer * d - b * r2_cramer) / Det
    delta_Q2 = (a * r2_cramer - r1_cramer * c) / Det
    
    ! =======================================================================
    ! BLOCO 5: ATUALIZAÇÃO E VERIFICAÇÃO DE ERRO
    ! =======================================================================
    ! Atualizar os valores antigos com as correções[cite: 60]
    Q1 = Q1 + delta_Q1
    Q2 = Q2 + delta_Q2
    
    ! O erro é o maior valor absoluto entre as correções[cite: 60]
    erro = max(abs(delta_Q1), abs(delta_Q2))
    
    ! Imprimir histórico da iteração na tela
    print '(I8, 3X, F14.8, 3X, F14.8, 3X, E14.6)', iteracao, Q1, Q2, erro
    
  end do

  ! =========================================================================
  ! BLOCO 6: PÓS-PROCESSAMENTO E VERIFICAÇÃO FÍSICA
  ! COPIE SE A QUESTÃO PEDIR PARA CONFERIR OS RESULTADOS FINAIS.
  ! =========================================================================
  if (erro <= tolerancia) then
    print *, '--------------------------------------------------------------'
    print *, 'Convergencia atingida com sucesso!'
    
    ! Calcular velocidades médias finais[cite: 61]
    U1 = (4.0d0 * Q1) / (pi * D1**2)
    U2 = (4.0d0 * Q2) / (pi * D2**2)
    
    ! Calcular perdas de carga finais[cite: 61]
    hL1 = R1 * Q1**2
    hL2 = R2 * Q2**2
    
    print *, ''
    print *, 'RESULTADOS FINAIS:'
    print *, 'Vazao Ramo 1 (Q1) = ', Q1, ' m3/s'
    print *, 'Vazao Ramo 2 (Q2) = ', Q2, ' m3/s'
    print *, 'Soma das Vazoes   = ', Q1 + Q2, ' m3/s (Deve ser 0.050)'
    print *, 'Velocidade 1 (U1) = ', U1, ' m/s'
    print *, 'Velocidade 2 (U2) = ', U2, ' m/s'
    print *, 'Perda Carga 1(hL1)= ', hL1, ' m'
    print *, 'Perda Carga 2(hL2)= ', hL2, ' m (Devem ser iguais)'
    
  else
    print *, 'O metodo nao convergiu dentro do limite de iteracoes.'
  end if

end program newton_raphson_sistemas_tubulacoes