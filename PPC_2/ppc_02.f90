program raizes_bairstow
  ! =========================================================================
  ! REGRA DE OURO: Sempre usar implicit none.
  ! =========================================================================
  implicit none

  ! =========================================================================
  ! BLOCO 1: DECLARAÇÃO DE VARIÁVEIS[cite: 114, 115]
  ! =========================================================================
  integer, parameter :: max_grau = 20
  integer :: grau_atual, iteracao, max_iter, i, num_raizes
  
  ! Vetores de coeficientes do polinomio[cite: 114, 115]
  real*8, dimension(0:max_grau) :: a, b, c
  
  ! Variaveis para o Metodo de Bairstow[cite: 114, 115]
  real*8 :: r, s, delta_r, delta_s
  real*8 :: det, erro_r, erro_s, tol
  
  ! Variaveis para extrair as raizes (Equacao Quadratica)[cite: 115]
  real*8 :: discriminante
  complex*16, dimension(max_grau) :: raizes ! As raizes podem ser complexas
  
  ! =========================================================================
  ! BLOCO 2: ENTRADA DE DADOS E CONFIGURAÇÃO
  ! =========================================================================
  ! Exemplo de validacao pedido na atividade (Polinomio de ordem 7)[cite: 118]
  ! Raizes conhecidas para este teste: x = -1, -2, 1, 2, 3, -3, 0 (Grau 7)
  ! f(x) = x^7 - 14x^5 + 49x^3 - 36x
  ! Se quiser testar o Polinomio Caracteristico do Massa-Mola, altere estes coeficientes!
  grau_atual = 7
  a = 0.0d0
  a(7) =  1.0d0
  a(6) =  0.0d0
  a(5) = -14.0d0
  a(4) =  0.0d0
  a(3) =  49.0d0
  a(2) =  0.0d0
  a(1) = -36.0d0
  a(0) =  0.0d0
  
  ! Parametros iterativos iniciais sugeridos[cite: 118]
  r = 1.0d0
  s = 1.0d0
  tol = 1.0d-6
  max_iter = 1000
  num_raizes = 0

  print *, 'Iniciando Metodo de Bairstow para Polinomio de grau ', grau_atual

  ! =========================================================================
  ! BLOCO 3: LOOP PRINCIPAL DE DEFLAÇÃO (EXTRAI AS RAÍZES 2 A 2)[cite: 115, 116]
  ! =========================================================================
  do while (grau_atual > 2)
    iteracao = 0
    erro_r = 1.0d0
    erro_s = 1.0d0
    
    ! 3.1: Loop de Convergencia Newton-Raphson para achar (r,s)[cite: 114, 115]
    do while ((erro_r > tol .or. erro_s > tol) .and. iteracao < max_iter)
      iteracao = iteracao + 1
      
      ! 3.1.1: Divisao Sintetica para calcular os b's (Restos)[cite: 114, 116]
      b(grau_atual) = a(grau_atual)
      b(grau_atual-1) = a(grau_atual-1) + r * b(grau_atual)
      do i = grau_atual-2, 0, -1
        b(i) = a(i) + r * b(i+1) + s * b(i+2)
      end do
      
      ! 3.1.2: Divisao Sintetica para calcular os c's (Derivadas Parciais)[cite: 115, 116]
      c(grau_atual) = b(grau_atual)
      c(grau_atual-1) = b(grau_atual-1) + r * c(grau_atual)
      do i = grau_atual-2, 1, -1
        c(i) = b(i) + r * c(i+1) + s * c(i+2)
      end do
      
      ! 3.1.3: Montagem e Solucao do Sistema Linear 2x2 por Cramer[cite: 114, 115, 116]
      ! [c2  c3] * [delta_r] = [-b1]
      ! [c1  c2]   [delta_s]   [-b0]
      det = c(2) * c(2) - c(3) * c(1)
      
      if (abs(det) < 1.0d-12) then
         print *, 'Erro: Determinante nulo. Tente outros valores iniciais de r e s.'
         stop
      end if
      
      ! Regra de Cramer para calcular os passos delta_r e delta_s[cite: 115]
      delta_r = (-b(1) * c(2) - (-b(0)) * c(3)) / det
      delta_s = (c(2) * (-b(0)) - c(1) * (-b(1))) / det
      
      ! Atualizacao das estimativas r e s[cite: 115, 116]
      r = r + delta_r
      s = s + delta_s
      
      ! Calculo dos erros percentuais para criterio de parada[cite: 116]
      if (abs(r) > 0.0d0) erro_r = abs(delta_r / r)
      if (abs(s) > 0.0d0) erro_s = abs(delta_s / s)
      
    end do ! Fim do laco Newton-Raphson
    
    ! 3.2: Extrair o par de raizes usando a Formula Quadratica[cite: 115, 116]
    ! x = (r +- sqrt(r^2 + 4s)) / 2
    discriminante = r**2 + 4.0d0 * s
    if (discriminante >= 0.0d0) then
       ! Raizes Reais
       num_raizes = num_raizes + 1
       raizes(num_raizes) = cmplx((r + sqrt(discriminante))/2.0d0, 0.0d0, 8)
       num_raizes = num_raizes + 1
       raizes(num_raizes) = cmplx((r - sqrt(discriminante))/2.0d0, 0.0d0, 8)
    else
       ! Raizes Complexas Conjugadas
       num_raizes = num_raizes + 1
       raizes(num_raizes) = cmplx(r/2.0d0, sqrt(abs(discriminante))/2.0d0, 8)
       num_raizes = num_raizes + 1
       raizes(num_raizes) = cmplx(r/2.0d0, -sqrt(abs(discriminante))/2.0d0, 8)
    end if
    
    ! 3.3: Deflacao - Reduzir a ordem do polinomio substituindo 'a' por 'b'[cite: 114, 115]
    grau_atual = grau_atual - 2
    do i = 0, grau_atual
       a(i) = b(i+2)
    end do
    
  end do ! Fim do laco de deflacao

  ! =========================================================================
  ! BLOCO 4: TRATAMENTO DO POLINÓMIO RESTANTE (GRAU 1 OU 2)[cite: 115, 116]
  ! =========================================================================
  if (grau_atual == 2) then
    ! Restou uma equacao do 2o grau pura: a(2)x^2 + a(1)x + a(0) = 0
    discriminante = a(1)**2 - 4.0d0 * a(2) * a(0)
    if (discriminante >= 0.0d0) then
       num_raizes = num_raizes + 1
       raizes(num_raizes) = cmplx((-a(1) + sqrt(discriminante))/(2.0d0*a(2)), 0.0d0, 8)
       num_raizes = num_raizes + 1
       raizes(num_raizes) = cmplx((-a(1) - sqrt(discriminante))/(2.0d0*a(2)), 0.0d0, 8)
    else
       num_raizes = num_raizes + 1
       raizes(num_raizes) = cmplx(-a(1)/(2.0d0*a(2)), sqrt(abs(discriminante))/(2.0d0*a(2)), 8)
       num_raizes = num_raizes + 1
       raizes(num_raizes) = cmplx(-a(1)/(2.0d0*a(2)), -sqrt(abs(discriminante))/(2.0d0*a(2)), 8)
    end if
  else if (grau_atual == 1) then
    ! Restou uma equacao linear: a(1)x + a(0) = 0[cite: 115, 116]
    num_raizes = num_raizes + 1
    raizes(num_raizes) = cmplx(-a(0)/a(1), 0.0d0, 8)
  end if

  ! =========================================================================
  ! BLOCO 5: IMPRESSÃO FINAL DAS RAÍZES
  ! =========================================================================
  print *, '---------------------------------------------------'
  print *, 'Raizes Encontradas:'
  do i = 1, num_raizes
    print *, 'Raiz ', i, ': ', raizes(i)
  end do

end program raizes_bairstow