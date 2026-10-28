program eliminacao_gaussiana_ingenua
  ! =========================================================================
  ! REGRA DE OURO: Sempre usar implicit none.
  ! =========================================================================
  implicit none
  
  ! =========================================================================
  ! BLOCO 1: DECLARAÇÃO DE VARIÁVEIS
  ! =========================================================================
  integer, parameter :: n = 5 ! O sistema possui 5 equacoes e 5 incognitas[cite: 70]
  
  ! Matrizes e vetores do sistema[cite: 71]
  real*8, dimension(n,n) :: A, A_orig
  real*8, dimension(n) :: b, b_orig, c, residuo
  
  ! Indices dos lacos iterativos[cite: 72]
  integer :: i, j, k
  real*8 :: fator, soma, erro_max
  
  ! Contadores de FLOPs (Custo Computacional)[cite: 75]
  integer :: flops_md ! Contador de multiplicacoes e divisoes
  integer :: flops_ss ! Contador de somas e subtracoes

  ! =========================================================================
  ! BLOCO 2: INICIALIZAÇÃO E MONTAGEM DO SISTEMA
  ! COPIE SE A QUESTÃO PEDIR PARA PREENCHER MATRIZES DIRETAMENTE
  ! =========================================================================
  ! Zerar os contadores antes de iniciar o algoritmo[cite: 75]
  flops_md = 0
  flops_ss = 0
  
  ! Montagem manual da Matriz A com os coeficientes dos balancos de massa[cite: 71]
  ! Inicializando tudo com zero por seguranca
  A = 0.0d0
  A(1,1) =  6.0d0; A(1,2) =  0.0d0; A(1,3) = -1.0d0; A(1,4) =  0.0d0; A(1,5) =  0.0d0
  A(2,1) = -3.0d0; A(2,2) =  3.0d0; A(2,3) =  0.0d0; A(2,4) =  0.0d0; A(2,5) =  0.0d0
  A(3,1) =  0.0d0; A(3,2) = -1.0d0; A(3,3) =  9.0d0; A(3,4) =  0.0d0; A(3,5) =  0.0d0
  A(4,1) =  0.0d0; A(4,2) = -1.0d0; A(4,3) = -8.0d0; A(4,4) = 11.0d0; A(4,5) = -2.0d0
  A(5,1) = -3.0d0; A(5,2) = -1.0d0; A(5,3) =  0.0d0; A(5,4) =  0.0d0; A(5,5) =  4.0d0
  
  ! Montagem do vetor b (Termos independentes baseados nas entradas)[cite: 71]
  b = 0.0d0
  b(1) = 50.0d0   ! Q01 * c01 = 5 * 10
  b(2) = 0.0d0
  b(3) = 160.0d0  ! Q03 * c03 = 8 * 20
  b(4) = 0.0d0
  b(5) = 0.0d0
  
  ! Guardar copias das versoes originais para calcular o residuo final[cite: 76]
  A_orig = A
  b_orig = b

  print *, 'Iniciando Eliminacao Gaussiana Ingenua...'

  ! =========================================================================
  ! BLOCO 3: ELIMINAÇÃO PROGRESSIVA (SEM PIVOTEAMENTO)
  ! COPIE SE A QUESTÃO EXIGIR O MÉTODO INGÊNUO E A CONTAGEM DE FLOPS[cite: 72, 73, 75]
  ! =========================================================================
  ! O laco externo 'k' escolhe a coluna do pivo[cite: 72]
  do k = 1, n-1
  
    ! Se o pivo for zero, o metodo ingenuo falha
    if (abs(A(k,k)) < 1.0d-12) then
       print *, 'ERRO: Pivo nulo detetado na coluna ', k
       print *, 'A eliminacao ingenua nao suporta pivos nulos.'
       stop
    end if
    
    ! O laco 'i' percorre as linhas abaixo do pivo[cite: 73]
    do i = k+1, n
      ! Calcula o multiplicador
      fator = A(i,k) / A(k,k)
      flops_md = flops_md + 1 ! Regista uma divisao[cite: 75]
      
      ! O laco 'j' atualiza as colunas à direita[cite: 73]
      do j = k+1, n
        A(i,j) = A(i,j) - fator * A(k,j)
        flops_md = flops_md + 1 ! Regista uma multiplicacao[cite: 75]
        flops_ss = flops_ss + 1 ! Regista uma subtracao[cite: 75]
      end do
      
      ! Atualiza o vetor do lado direito
      b(i) = b(i) - fator * b(k)
      flops_md = flops_md + 1 ! Regista uma multiplicacao[cite: 75]
      flops_ss = flops_ss + 1 ! Regista uma subtracao[cite: 75]
      
      ! Forca o zero abaixo do pivo
      A(i,k) = 0.0d0
    end do
  end do

  ! =========================================================================
  ! BLOCO 4: SUBSTITUIÇÃO RETROATIVA
  ! COPIE ESTE BLOCO PARA CALCULAR AS INCÓGNITAS A PARTIR DO SISTEMA TRIANGULAR[cite: 75]
  ! =========================================================================
  ! Calcula a ultima incognita
  c(n) = b(n) / A(n,n)
  flops_md = flops_md + 1 ! Regista uma divisao[cite: 75]
  
  ! Caminha de baixo para cima[cite: 75]
  do i = n-1, 1, -1
    soma = b(i)
    
    do j = i+1, n
      soma = soma - A(i,j) * c(j)
      flops_md = flops_md + 1 ! Regista uma multiplicacao[cite: 75]
      flops_ss = flops_ss + 1 ! Regista uma subtracao[cite: 75]
    end do
    
    c(i) = soma / A(i,i)
    flops_md = flops_md + 1 ! Regista uma divisao[cite: 75]
  end do

  ! =========================================================================
  ! BLOCO 5: CÁLCULO DO RESÍDUO E VERIFICAÇÃO FÍSICA
  ! =========================================================================
  ! r = A_orig * c - b_orig[cite: 76]
  erro_max = 0.0d0
  
  do i = 1, n
    soma = 0.0d0
    do j = 1, n
      soma = soma + A_orig(i,j) * c(j)
    end do
    residuo(i) = soma - b_orig(i)
    erro_max = max(erro_max, abs(residuo(i))) ! Medida global do residuo[cite: 76]
  end do

  ! =========================================================================
  ! BLOCO 6: IMPRESSÃO DOS RESULTADOS FINAIS E DE FLOPS[cite: 75, 77]
  ! =========================================================================
  print *, '--------------------------------------------------------------'
  print *, 'Resultados das Concentracoes (mg/m3):'
  do i = 1, n
    print *, 'c(', i, ') = ', c(i)
  end do
  print *, '--------------------------------------------------------------'
  print *, 'Verificacao de Qualidade Numérica:'
  print *, 'Residuo Maximo = ', erro_max, ' (Deve ser proximo de zero)'
  print *, '--------------------------------------------------------------'
  print *, 'Custo Computacional (FLOPs):'
  print *, 'Multiplicacoes/Divisoes = ', flops_md
  print *, 'Somas/Subtracoes        = ', flops_ss
  print *, 'Total de FLOPs          = ', flops_md + flops_ss

end program eliminacao_gaussiana_ingenua