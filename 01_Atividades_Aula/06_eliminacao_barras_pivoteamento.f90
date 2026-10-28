program eliminacao_barras_pivoteamento
  ! =========================================================================
  ! REGRA DE OURO: Sempre usar implicit none.
  ! =========================================================================
  implicit none
  
  ! =========================================================================
  ! BLOCO 1: DECLARAÇÃO DE VARIÁVEIS
  ! =========================================================================
  integer, parameter :: n = 5 ! Sistema de 5 barras[cite: 81]
  
  ! Matrizes e vetores do sistema linear[cite: 81]
  real*8, dimension(n,n) :: A, A_orig
  real*8, dimension(n) :: F_calc, b, b_orig, residuo
  real*8, dimension(n) :: F_ref ! Forcas de referencia para gerar o vetor b[cite: 81]
  
  ! Variaveis para geometria e trigonometria[cite: 80, 81]
  real*8 :: theta1_deg, theta2_deg, phi_deg
  real*8 :: th1, th2, phi, pi
  
  ! Variaveis para o algoritmo de Gauss com Pivoteamento[cite: 82, 83]
  integer :: i, j, k, p
  real*8 :: fator, soma, temp, limite_pivo, norma_inf_residuo

  ! =========================================================================
  ! BLOCO 2: ENTRADA DE DADOS E CONVERSÃO DE ÂNGULOS
  ! COPIE SE A QUESTÃO ENVOLVER SENOS/COSSENOS E GRAUS[cite: 80]
  ! =========================================================================
  pi = acos(-1.0d0)
  
  ! Angulos em graus (Cenario A - barras afastadas)[cite: 80]
  theta1_deg = 15.0d0
  theta2_deg = 65.0d0
  phi_deg = 45.0d0
  
  ! Conversao de graus para radianos (Obrigatorio no Fortran)[cite: 80]
  th1 = theta1_deg * pi / 180.0d0
  th2 = theta2_deg * pi / 180.0d0
  phi = phi_deg * pi / 180.0d0

  ! =========================================================================
  ! BLOCO 3: MONTAGEM DA MATRIZ DO SISTEMA
  ! =========================================================================
  A = 0.0d0 ! Zera tudo primeiro
  
  ! Preenchimento das equacoes de equilibrio nos nos A, B e C[cite: 81]
  ! Equacao 1 (Horizontal no no A)
  A(1,1) = sin(th1); A(1,2) = sin(th2); A(1,3) = 1.0d0; A(1,4) = 0.0d0; A(1,5) = 0.0d0
  ! Equacao 2 (Vertical no no A)
  A(2,1) = cos(th1); A(2,2) = cos(th2); A(2,3) = 0.0d0; A(2,4) = 0.0d0; A(2,5) = 0.0d0
  ! Equacao 3 (Horizontal no no B)
  A(3,1) = 0.0d0; A(3,2) = 0.0d0; A(3,3) = -1.0d0; A(3,4) = cos(phi); A(3,5) = 0.0d0
  ! Equacao 4 (Vertical no no B)
  A(4,1) = 0.0d0; A(4,2) = 0.0d0; A(4,3) = 0.0d0; A(4,4) = -sin(phi); A(4,5) = -1.0d0
  ! Equacao 5 (Vertical no no C)
  A(5,1) = 0.0d0; A(5,2) = 0.0d0; A(5,3) = 0.0d0; A(5,4) = 0.0d0; A(5,5) = 1.0d0

  ! Geracao do vetor de cargas b a partir da solucao de referencia[cite: 81]
  F_ref(1) = 100.0d0; F_ref(2) = -80.0d0; F_ref(3) = 20.0d0
  F_ref(4) = 35.0d0;  F_ref(5) = 15.0d0
  
  ! Multiplicacao b = A * F_ref para gerar um cenario consistente[cite: 81]
  do i = 1, n
    soma = 0.0d0
    do j = 1, n
      soma = soma + A(i,j) * F_ref(j)
    end do
    b(i) = soma
  end do

  ! Guardar originais para o residuo[cite: 83]
  A_orig = A
  b_orig = b
  
  limite_pivo = 1.0d-12 ! Criterio para pivo nulo[cite: 82]

  print *, 'Iniciando Eliminacao Gaussiana com Pivoteamento Parcial...'

  ! =========================================================================
  ! BLOCO 4: ELIMINAÇÃO COM PIVOTEAMENTO PARCIAL (O CORAÇÃO DO ALGORITMO)[cite: 82]
  ! =========================================================================
  do k = 1, n-1
    
    ! 4.1 Identificar a linha pivo (maior valor em modulo na coluna k)[cite: 82]
    p = k
    do i = k+1, n
      if (abs(A(i,k)) > abs(A(p,k))) then
        p = i
      end if
    end do
    
    ! 4.2 Parar se o pivo for numericamente zero[cite: 82]
    if (abs(A(p,k)) <= limite_pivo) then
      print *, 'ERRO: Sistema singular. Pivo muito pequeno.'
      stop
    end if
    
    ! 4.3 Troca de linhas em A e b (Pivoteamento efetivo)[cite: 82]
    if (p /= k) then
      print *, 'Iteracao ', k, ' -> Trocando linha ', k, ' com linha ', p
      do j = 1, n
        temp = A(k,j)
        A(k,j) = A(p,j)
        A(p,j) = temp
      end do
      temp = b(k)
      b(k) = b(p)
      b(p) = temp
    end if
    
    ! 4.4 Eliminacao progressiva[cite: 82]
    do i = k+1, n
      fator = A(i,k) / A(k,k)
      do j = k+1, n
        A(i,j) = A(i,j) - fator * A(k,j)
      end do
      b(i) = b(i) - fator * b(k)
      A(i,k) = 0.0d0
    end do
  end do

  ! =========================================================================
  ! BLOCO 5: SUBSTITUIÇÃO REGRESSIVA[cite: 83]
  ! =========================================================================
  F_calc(n) = b(n) / A(n,n)
  
  do i = n-1, 1, -1
    soma = b(i)
    do j = i+1, n
      soma = soma - A(i,j) * F_calc(j)
    end do
    F_calc(i) = soma / A(i,i)
  end do

  ! =========================================================================
  ! BLOCO 6: CÁLCULO DA NORMA INFINITO DO RESÍDUO
  ! COPIE SE PEDIR "NORMA INFINITO" OU "ERRO MAXIMO"[cite: 83]
  ! =========================================================================
  norma_inf_residuo = 0.0d0
  do i = 1, n
    soma = 0.0d0
    do j = 1, n
      soma = soma + A_orig(i,j) * F_calc(j)
    end do
    residuo(i) = soma - b_orig(i)
    
    ! A norma infinito e o maior valor absoluto do vetor residuo[cite: 83]
    if (abs(residuo(i)) > norma_inf_residuo) then
      norma_inf_residuo = abs(residuo(i))
    end if
  end do

  ! =========================================================================
  ! BLOCO 7: IMPRESSÃO FINAL
  ! =========================================================================
  print *, '--------------------------------------------------------------'
  print *, 'Esforcos Calculados nas Barras (F_calc):'
  do i = 1, n
    print *, 'F(', i, ') = ', F_calc(i)
  end do
  print *, '--------------------------------------------------------------'
  print *, 'Norma Infinito do Residuo = ', norma_inf_residuo

end program eliminacao_barras_pivoteamento