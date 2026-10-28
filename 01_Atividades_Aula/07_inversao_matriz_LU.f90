program inversao_matriz_LU
  ! =========================================================================
  ! REGRA DE OURO: Sempre usar implicit none.
  ! =========================================================================
  implicit none

  ! =========================================================================
  ! BLOCO 1: DECLARAÇÃO DE VARIÁVEIS
  ! =========================================================================
  integer, parameter :: n = 5 ! Sistema de 5 graus de liberdade
  
  ! Matrizes
  real*8, dimension(n,n) :: K_mat  ! Matriz de Rigidez Original (Renomeada para evitar conflito com 'k')
  real*8, dimension(n,n) :: L      ! Matriz Triangular Inferior (Lower)
  real*8, dimension(n,n) :: U      ! Matriz Triangular Superior (Upper)
  real*8, dimension(n,n) :: C      ! Matriz Inversa (Flexibilidade)
  real*8, dimension(n,n) :: M      ! Matriz para verificar se K*C = I
  
  ! Vetores
  real*8, dimension(n) :: e_vec    ! Vetor coluna da matriz identidade
  real*8, dimension(n) :: y        ! Vetor intermediario (L*y = e)
  real*8, dimension(n) :: x        ! Vetor solucao (U*x = y)
  real*8, dimension(n) :: F        ! Vetor de forcas de carregamento
  real*8, dimension(n) :: u_vec    ! Vetor de deslocamentos
  
  ! Constantes fisicas das molas (N/m)
  real*8 :: k0, k1, k2, k3, k4, k5
  
  ! Variaveis de laco e calculo
  integer :: i, j, k, coluna
  real*8 :: fator, soma, erro_max, erro_sim

  ! =========================================================================
  ! BLOCO 2: INICIALIZAÇÃO DA MATRIZ DE RIGIDEZ
  ! =========================================================================
  k0 = 1200.0d0; k1 =  800.0d0; k2 = 1500.0d0
  k3 = 1000.0d0; k4 = 1800.0d0; k5 =  900.0d0
  
  ! Preenchimento da matriz tridiagonal de rigidez
  K_mat = 0.0d0
  K_mat(1,1) = k0 + k1; K_mat(1,2) = -k1
  K_mat(2,1) = -k1;     K_mat(2,2) = k1 + k2; K_mat(2,3) = -k2
  K_mat(3,2) = -k2;     K_mat(3,3) = k2 + k3; K_mat(3,4) = -k3
  K_mat(4,3) = -k3;     K_mat(4,4) = k3 + k4; K_mat(4,5) = -k4
  K_mat(5,4) = -k4;     K_mat(5,5) = k4 + k5

  ! Inicializar a matriz Identidade em L e copiar K_mat para U
  L = 0.0d0
  do i = 1, n
    L(i,i) = 1.0d0
  end do
  U = K_mat

  ! =========================================================================
  ! BLOCO 3: DECOMPOSIÇÃO LU (FEITA APENAS UMA VEZ)
  ! =========================================================================
  ! K_mat = L * U
  do k = 1, n-1
    do i = k+1, n
      ! Fator de eliminacao
      fator = U(i,k) / U(k,k)
      ! O fator e guardado na matriz L
      L(i,k) = fator
      
      ! Atualiza a linha da matriz U
      do j = k+1, n
        U(i,j) = U(i,j) - fator * U(k,j)
      end do
      ! Zera os elementos abaixo do pivo
      U(i,k) = 0.0d0
    end do
  end do

  ! =========================================================================
  ! BLOCO 4: CÁLCULO DA MATRIZ INVERSA COLUNA A COLUNA
  ! =========================================================================
  C = 0.0d0
  
  ! Repete o processo para cada coluna da matriz identidade
  do coluna = 1, n
    
    ! 4.1 Cria o vetor 'e' (coluna da matriz identidade)
    e_vec = 0.0d0
    e_vec(coluna) = 1.0d0
    
    ! 4.2 Resolucao do sistema triangular INFERIOR (L*y = e)
    ! Substituicao Progressiva (de cima para baixo)
    do i = 1, n
      soma = e_vec(i)
      do j = 1, i-1
        soma = soma - L(i,j) * y(j)
      end do
      y(i) = soma / L(i,i)
    end do
    
    ! 4.3 Resolucao do sistema triangular SUPERIOR (U*x = y)
    ! Substituicao Regressiva (de baixo para cima)
    x(n) = y(n) / U(n,n)
    do i = n-1, 1, -1
      soma = y(i)
      do j = i+1, n
        soma = soma - U(i,j) * x(j)
      end do
      x(i) = soma / U(i,i)
    end do
    
    ! 4.4 Armazena a solucao 'x' como uma coluna da Matriz Inversa 'C'
    do i = 1, n
      C(i,coluna) = x(i)
    end do
    
  end do

  ! =========================================================================
  ! BLOCO 5: VERIFICAÇÕES NUMÉRICAS
  ! =========================================================================
  
  ! 5.1 Verificar se K_mat*C = Identidade (Erro maximo deve ser zero)
  erro_max = 0.0d0
  do i = 1, n
    do j = 1, n
      soma = 0.0d0
      do k = 1, n
        soma = soma + K_mat(i,k) * C(k,j)
      end do
      M(i,j) = soma
      
      ! Verifica contra a identidade real
      if (i == j) then
        erro_max = max(erro_max, abs(M(i,j) - 1.0d0))
      else
        erro_max = max(erro_max, abs(M(i,j)))
      end if
    end do
  end do
  
  ! 5.2 Verificar Simetria da matriz C (C_ij = C_ji)
  erro_sim = 0.0d0
  do i = 1, n
    do j = 1, n
      erro_sim = max(erro_sim, abs(C(i,j) - C(j,i)))
    end do
  end do

  ! =========================================================================
  ! BLOCO 6: APLICAÇÃO DA MATRIZ DE FLEXIBILIDADE (CÁLCULO DE u)
  ! =========================================================================
  ! Carregamento dado no exercicio
  F(1) = 0.0d0; F(2) = 100.0d0; F(3) = -50.0d0; F(4) = 25.0d0; F(5) = 0.0d0
  
  ! Deslocamentos u = C * F
  do i = 1, n
    soma = 0.0d0
    do j = 1, n
      soma = soma + C(i,j) * F(j)
    end do
    u_vec(i) = soma
  end do

  ! =========================================================================
  ! IMPRESSÃO FINAL
  ! =========================================================================
  print *, '--- VERIFICACOES NUMERICAS ---'
  print *, 'Erro Maximo de Identidade (K*C = I): ', erro_max
  print *, 'Erro de Simetria (C_ij = C_ji): ', erro_sim
  print *, ''
  print *, '--- DESLOCAMENTOS DA ESTRUTURA (u) ---'
  do i = 1, n
    print *, 'u(', i, ') = ', u_vec(i), ' m'
  end do

end program inversao_matriz_LU