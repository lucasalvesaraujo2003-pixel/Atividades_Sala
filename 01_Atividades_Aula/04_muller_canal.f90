program muller_canal_critico
  ! =========================================================================
  ! REGRA DE OURO: 'implicit none' para evitar bugs com variáveis.
  ! =========================================================================
  implicit none
  
  ! =========================================================================
  ! BLOCO 1: DECLARAÇÃO DE VARIÁVEIS
  ! =========================================================================
  ! Constantes físicas do canal[cite: 63]
  real*8 :: b, z, g, Q, Ac, Tc, Fr2
  
  ! Variáveis iterativas do Método de Müller (As estimativas)[cite: 64, 65]
  real*8 :: H0, H1, H2, H3
  real*8 :: f0, f1, f2, f3
  real*8 :: erro, tolerancia
  integer :: iteracao, max_iter
  
  ! Variáveis internas do Método de Müller[cite: 65, 66]
  ! NOTA: No método de Müller, o discriminante (raiz quadrada) pode ficar negativo
  ! transitoriamente. Por isso, as variáveis de cálculo interno devem ser COMPLEXAS!
  complex*16 :: h0_c, h1_c, delta0, delta1, a_c, b_M, c_c
  complex*16 :: discriminante, E_plus, E_minus, E_den, delta_H
  complex*16 :: H0_c, H1_c, H2_c
  
  ! =========================================================================
  ! BLOCO 2: ENTRADA DE DADOS E CONFIGURAÇÃO[cite: 63, 64]
  ! =========================================================================
  b = 4.0d0
  z = 1.5d0
  g = 9.81d0
  Q = 20.0d0  ! Altere a vazão aqui para a Parte 5 (varredura de 5 a 40 m3/s)[cite: 67]
  
  ! Estimativas iniciais fornecidas pelo problema[cite: 64]
  H0 = 0.5d0
  H1 = 1.0d0
  H2 = 1.5d0
  
  ! Critério de convergência[cite: 66]
  tolerancia = 1.0d-8
  max_iter = 100
  erro = 1.0d0
  iteracao = 0
  
  print *, 'Iniciando o Metodo de Muller...'
  print *, ' Iter |      H0      |      H1      |      H2      |      H3      |     f(H3)    |     Erro     '
  print *, '------------------------------------------------------------------------------------------------'

  ! =========================================================================
  ! BLOCO 3: LOOP PRINCIPAL DO MÉTODO DE MÜLLER
  ! COPIE ESTE BLOCO SE A QUESTÃO PEDIR PARA IMPLEMENTAR MÜLLER[cite: 64, 65, 66]
  ! =========================================================================
  do while (erro > tolerancia .and. iteracao < max_iter)
    iteracao = iteracao + 1
    
    ! 3.1 Converter as estimativas reais para complexas (necessário para o cálculo seguro)
    H0_c = cmplx(H0, 0.0d0, 8)
    H1_c = cmplx(H1, 0.0d0, 8)
    H2_c = cmplx(H2, 0.0d0, 8)
    
    ! 3.2 Avaliar a função nos três pontos[cite: 64]
    f0 = f_polinomio(H0, Q, b, z, g)
    f1 = f_polinomio(H1, Q, b, z, g)
    f2 = f_polinomio(H2, Q, b, z, g)
    
    ! 3.3 Calcular as diferenças e inclinações (Secantes)[cite: 64, 65]
    h0_c = H1_c - H0_c
    h1_c = H2_c - H1_c
    
    delta0 = (f1 - f0) / h0_c
    delta1 = (f2 - f1) / h1_c
    
    ! 3.4 Calcular os coeficientes da parábola (a, b_M, c)[cite: 65]
    a_c = (delta1 - delta0) / (h1_c + h0_c)
    b_M = a_c * h1_c + delta1
    c_c = cmplx(f2, 0.0d0, 8)
    
    ! 3.5 Calcular a raiz (Discriminante)[cite: 65]
    discriminante = sqrt(b_M**2 - 4.0d0 * a_c * c_c)
    
    ! 3.6 Escolher o denominador de maior módulo para evitar cancelamento catastrófico[cite: 65]
    E_plus = b_M + discriminante
    E_minus = b_M - discriminante
    
    if (abs(E_plus) > abs(E_minus)) then
      E_den = E_plus
    else
      E_den = E_minus
    end if
    
    ! 3.7 Calcular a correção (delta_H) e a nova raiz (H3)[cite: 65]
    delta_H = -2.0d0 * c_c / E_den
    
    ! A nova aproximação (Consideramos apenas a parte real para o fluxo físico)
    H3 = H2 + real(delta_H)
    
    ! 3.8 Avaliar o erro e atualizar pontos[cite: 65, 66]
    erro = abs(H3 - H2)
    f3 = f_polinomio(H3, Q, b, z, g)
    
    ! Imprimir o estado atual da iteração[cite: 66]
    print '(I5, 2X, F12.6, 2X, F12.6, 2X, F12.6, 2X, F12.6, 2X, E12.4, 2X, E12.4)', &
          iteracao, H0, H1, H2, H3, f3, erro
          
    ! Atualiza a "janela" de pontos para a próxima iteração: H0 vira H1, H1 vira H2, etc.[cite: 65]
    H0 = H1
    H1 = H2
    H2 = H3
    
  end do

  ! =========================================================================
  ! BLOCO 4: VERIFICAÇÃO FÍSICA (NÚMERO DE FROUDE)
  ! COPIE SE A QUESTÃO PEDIR PARA VALIDAR AS CONDIÇÕES FÍSICAS DA RAIZ[cite: 66, 67]
  ! =========================================================================
  if (erro <= tolerancia) then
    print *, '------------------------------------------------------------------------------------------------'
    print *, 'Convergencia atingida!'
    print *, 'Profundidade Critica (Hc) = ', H3, ' m'
    
    ! Verifica se a raiz faz sentido físico (deve ser maior que zero)[cite: 67]
    if (H3 > 0.0d0) then
       ! Calcular Propriedades da Secção[cite: 62, 63, 66]
       Ac = b*H3 + z*(H3**2)
       Tc = b + 2.0d0*z*H3
       
       ! Calcular Quadrado do Número de Froude (Fr^2 = Q^2 * Tc / (g * Ac^3))[cite: 66]
       Fr2 = (Q**2 * Tc) / (g * Ac**3)
       
       print *, ''
       print *, '=== VERIFICACAO FISICA ==='
       print *, 'Area Molhada (Ac) = ', Ac, ' m2'
       print *, 'Largura Sup. (Tc) = ', Tc, ' m'
       print *, 'Froude ao quadrado(Fr^2) = ', Fr2, ' (Deve ser aprox. 1.0)'
    else
       print *, 'ATENCAO: A raiz encontrada eh negativa. Fisicamente inadmissivel!'
    end if
  else
    print *, 'O metodo nao convergiu dentro do limite de iteracoes.'
  end if

contains

  ! =========================================================================
  ! FUNÇÃO AUXILIAR: POLINÓMIO DE SEXTA ORDEM
  ! Esta secção fica *dentro* do programa, antes do 'end program'[cite: 63, 66].
  ! =========================================================================
  real*8 function f_polinomio(H, Q_val, b_val, z_val, g_val)
    implicit none
    real*8, intent(in) :: H, Q_val, b_val, z_val, g_val
    
    ! Polinómio de grau 6[cite: 64]
    f_polinomio = g_val * (z_val**3) * (H**6) &
                + 3.0d0 * g_val * b_val * (z_val**2) * (H**5) &
                + 3.0d0 * g_val * (b_val**2) * z_val * (H**4) &
                + g_val * (b_val**3) * (H**3) &
                - 2.0d0 * (Q_val**2) * z_val * H &
                - (Q_val**2) * b_val
  end function f_polinomio

end program muller_canal_critico