program calor_gauss_seidel
  ! =========================================================================
  ! REGRA DE OURO: Sempre usar implicit none.
  ! =========================================================================
  implicit none

  ! =========================================================================
  ! BLOCO 1: DECLARAÇÃO DE VARIÁVEIS
  ! =========================================================================
  ! n = numero de nos INTERNOS em cada direcao[cite: 94]
  ! Se a questao pedir malha mais fina, mude para 10 ou 20 aqui[cite: 95]
  integer, parameter :: n = 5 
  
  ! Matriz de Temperaturas. Tamanho (n+2)x(n+2) para incluir as 4 bordas[cite: 94].
  real*8, dimension(n+2, n+2) :: T 
  
  ! Variaveis de laco e controle iterativo
  integer :: i, j, iteracao, max_iter
  real*8 :: erro, tolerancia, antigo, diferenca
  
  ! Variaveis para calculo do residuo final[cite: 95]
  real*8 :: residuo, max_residuo

  ! =========================================================================
  ! BLOCO 2: CONDIÇÕES INICIAIS E DE CONTORNO (FRONTEIRAS)
  ! COPIE SE A QUESTÃO PEDIR TEMPERATURAS FIXAS NAS BORDAS
  ! =========================================================================
  ! 2.1 Inicializar TODOS os elementos da placa com 20 graus[cite: 94]
  T = 20.0d0
  
  ! 2.2 Fixar a borda inferior (j=1) a 100 graus (apenas nos internos i=2 ate n+1)[cite: 94]
  ! O eixo 'i' representa a horizontal (x) e o 'j' a vertical (y)[cite: 94]
  do i = 2, n+1
    T(i, 1) = 100.0d0
  end do
  
  ! Configuração do critério de parada[cite: 94]
  tolerancia = 1.0d-6
  max_iter = 10000
  erro = 1.0d0
  iteracao = 0

  print *, 'Iniciando Metodo de Gauss-Seidel...'

  ! =========================================================================
  ! BLOCO 3: LOOP PRINCIPAL DE GAUSS-SEIDEL
  ! COPIE SE A QUESTÃO PEDIR MÉTODO ITERATIVO POR DIFERENÇAS FINITAS
  ! =========================================================================
  do while (erro > tolerancia .and. iteracao < max_iter)
    erro = 0.0d0 ! Zera o erro a cada nova varredura da malha
    
    ! Percorre apenas os NOS INTERNOS (de 2 ate n+1) nas duas direcoes[cite: 94]
    do j = 2, n+1
      do i = 2, n+1
        
        ! Salva a temperatura antiga para calcular a diferenca depois
        antigo = T(i,j)
        
        ! Equacao de Laplace discretizada (Media dos 4 vizinhos)[cite: 93]
        ! Como eh Gauss-Seidel, T(i-1,j) e T(i,j-1) ja usam valores atualizados desta mesma iteracao!
        T(i,j) = (T(i-1,j) + T(i+1,j) + T(i,j-1) + T(i,j+1)) / 4.0d0
        
        ! Calcula a diferenca absoluta
        diferenca = abs(T(i,j) - antigo)
        
        ! O erro da iteracao eh a MAIOR mudanca que ocorreu em toda a placa[cite: 94]
        if (diferenca > erro) then
          erro = diferenca
        end if
        
      end do
    end do
    
    iteracao = iteracao + 1
  end do

  ! =========================================================================
  ! BLOCO 4: CONFERÊNCIA ADICIONAL (CÁLCULO DO RESÍDUO)[cite: 95]
  ! =========================================================================
  max_residuo = 0.0d0
  do j = 2, n+1
    do i = 2, n+1
      ! Formula do residuo pedida na atividade[cite: 95]: r = 4*T - (soma dos vizinhos)
      residuo = abs(4.0d0 * T(i,j) - (T(i-1,j) + T(i+1,j) + T(i,j-1) + T(i,j+1)))
      if (residuo > max_residuo) then
        max_residuo = residuo
      end if
    end do
  end do

  ! =========================================================================
  ! BLOCO 5: IMPRESSÃO DOS RESULTADOS[cite: 95]
  ! =========================================================================
  if (erro <= tolerancia) then
    print *, 'Convergencia atingida com sucesso!'
    print *, 'Numero de Iteracoes = ', iteracao
    print *, 'Erro Final Maximo   = ', erro
    print *, 'Residuo Maximo      = ', max_residuo
    print *, '---------------------------------------------------'
    print *, 'Temperatura no centro da placa T(4,4) = ', T(4,4)
  else
    print *, 'Metodo parou sem convergir pelo limite de iteracoes.'
  end if

  ! =========================================================================
  ! BLOCO 6: RESPOSTAS ÀS QUESTÕES TEÓRICAS (COMO COMENTÁRIOS)[cite: 95]
  ! =========================================================================
  ! a) Por que as temperaturas das bordas nao sao atualizadas?
  ! R: Porque elas representam as "Condicoes de Contorno" do problema (Dirichlet). 
  ! Sao temperaturas fixas mantidas pelo ambiente e que forcam o calor a fluir pela placa.
  !
  ! b) O que diferencia a atualizacao de Gauss-Seidel da de Jacobi?
  ! R: No metodo de Jacobi, usamos apenas os valores da iteracao anterior (k) para 
  ! calcular todos os nos da iteracao nova (k+1). No Gauss-Seidel, mal calculamos 
  ! um valor e ja o utilizamos instantaneamente no calculo do proximo no vizinho, 
  ! o que acelera a convergencia.
  !
  ! c) Como o numero de iteracoes muda quando a malha eh refinada?
  ! R: Quanto mais fina a malha (maior 'n'), MAIOR o numero de iteracoes necessario. 
  ! Isso acontece porque a informacao termica (condicoes de contorno) demora mais 
  ! "passos" para se propagar das bordas ate o centro do dominio.
  
end program calor_gauss_seidel