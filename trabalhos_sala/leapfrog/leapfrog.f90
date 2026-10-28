! ==============================================================================
! ANÁLISE OBRIGATÓRIA - CONSERVAÇÃO DE ENERGIA
! ==============================================================================
!
! 1. O método de Euler conserva a energia total do sistema?
! Não. Ao analisarmos o gráfico gerado (energia.png), fica evidente que a linha 
! correspondente ao método de Euler não se mantém reta. A energia do sistema vai 
! se alterando (geralmente aumentando) à medida que o tempo passa, o que prova 
! que ele falha em conservar a energia mecânica.
!
! 2. O que acontece com as órbitas quando o método de Euler é usado por muitos 
! passos de tempo?
! Como o método de Euler "injeta" erros na energia a cada passo de tempo, as 
! órbitas perdem a estabilidade. Em vez de formarem ciclos fechados (como elipses 
! perfeitas que se repetem), os corpos começam a espiralar e a se afastar. Se a 
! simulação rodar por muito tempo, eles simplesmente se perdem no espaço.
!
! 3. O método Leapfrog apresenta melhor comportamento na conservação de energia?
! Sim, o comportamento é infinitamente melhor. Olhando novamente o gráfico de 
! energia, a curva do Leapfrog permanece o tempo todo oscilando minimamente em 
! torno do valor inicial. Ele não sofre daquele erro cumulativo do Euler, 
! mantendo a energia praticamente constante do começo ao fim.
!
! 4. Qual método é mais adequado para esse problema conservativo? Justifique.
! O método Leapfrog é o mais adequado. A principal justificativa é que problemas 
! gravitacionais são sistemas conservativos ideais, ou seja, não há perda de 
! energia. O Leapfrog foi matematicamente construído para respeitar a física de 
! sistemas que conservam energia (ele pertence a uma classe chamada de 
! integradores "simpléticos"). Isso garante que os dois corpos mantenham suas 
! órbitas estáveis e fechadas, independentemente de quantos passos de tempo 
! você simule.
!
! ==============================================================================

! ==============================================================================
! PROGRAMA PRINCIPAL
! ==============================================================================
program problema_dois_corpos_simples
    implicit none

    ! Declaração de variáveis simples
    real(8) :: m1, m2, G, dt, t_max
    integer :: n_steps
    real(8) :: R0(4), V0(4)

    ! 1. Configuração dos parâmetros iniciais
    m1 = 1.0d0
    m2 = 1.0d0
    G  = 1.0d0
    
    dt = 1.0d-3
    t_max = 5.0d0
    n_steps = int(t_max / dt)

    ! 2. Condições Iniciais
    ! Posições: R0 = [x1, y1, x2, y2]
    R0(1) = -0.5d0;  R0(2) = 0.0d0
    R0(3) =  0.5d0;  R0(4) = 0.0d0
    
    ! Velocidades: V0 = [vx1, vy1, vx2, vy2]
    V0(1) =  0.0d0;  V0(2) = -0.5d0
    V0(3) =  0.0d0;  V0(4) =  0.5d0

    ! 3. Executar os métodos
    print *, "Executando o metodo de Euler..."
    call metodo_euler(R0, V0, m1, m2, G, dt, n_steps)

    print *, "Executando o metodo Leapfrog..."
    call metodo_leapfrog(R0, V0, m1, m2, G, dt, n_steps)

    print *, "Fim do programa! Arquivos gerados com sucesso."

end program problema_dois_corpos_simples


! ==============================================================================
! SUBROTINA 1: Método de Euler (Simples)
! ==============================================================================
subroutine metodo_euler(R0, V0, m1, m2, G, dt, n_steps)
    implicit none
    
    ! Entradas
    real(8), intent(in) :: R0(4), V0(4)
    real(8), intent(in) :: m1, m2, G, dt
    integer, intent(in) :: n_steps
    
    ! Variáveis locais
    real(8) :: R(4), V(4), A(4)
    real(8) :: E, t
    integer :: i
    
    ! Inicializa o tempo e as variáveis com os valores iniciais
    t = 0.0d0
    R = R0
    V = V0
    
    ! Abre o arquivo para salvar os dados
    open(unit=10, file='euler.dat', status='replace')
    
    ! Loop de integração (Passo a passo)
    do i = 0, n_steps
        ! Calcula a energia atual
        call calcula_energia(R, V, m1, m2, G, E)
        
        ! Salva o instante atual no arquivo
        ! Colunas: t, x1, y1, vx1, vy1, x2, y2, vx2, vy2, E
        write(10, '(10F15.6)') t, R(1), R(2), V(1), V(2), R(3), R(4), V(3), V(4), E
        
        ! --- Fórmulas do Euler ---
        ! 1. Descobre a aceleração atual
        call calcula_aceleracao(R, m1, m2, G, A)
        
        ! 2. Atualiza a posição (r = r + dt * v)
        R = R + dt * V
        
        ! 3. Atualiza a velocidade (v = v + dt * a)
        V = V + dt * A
        
        ! Atualiza o tempo para a próxima rodada
        t = t + dt
    end do
    
    close(10)
end subroutine metodo_euler


! ==============================================================================
! SUBROTINA 2: Método Leapfrog (Simples)
! ==============================================================================
subroutine metodo_leapfrog(R0, V0, m1, m2, G, dt, n_steps)
    implicit none
    
    ! Entradas
    real(8), intent(in) :: R0(4), V0(4)
    real(8), intent(in) :: m1, m2, G, dt
    integer, intent(in) :: n_steps
    
    ! Variáveis locais
    real(8) :: R(4), V(4), A(4), V_meio(4)
    real(8) :: E, t
    integer :: i
    
    ! Inicializa
    t = 0.0d0
    R = R0
    V = V0
    
    ! Abre o arquivo para salvar os dados
    open(unit=20, file='leapfrog.dat', status='replace')
    
    ! Loop de integração
    do i = 0, n_steps
        call calcula_energia(R, V, m1, m2, G, E)
        write(20, '(10F15.6)') t, R(1), R(2), V(1), V(2), R(3), R(4), V(3), V(4), E
        
        ! --- Fórmulas do Leapfrog ---
        ! 1. Calcula a aceleração com a posição atual
        call calcula_aceleracao(R, m1, m2, G, A)
        
        ! 2. Calcula a velocidade no "meio do passo" (n+1/2)
        V_meio = V + 0.5d0 * dt * A
        
        ! 3. Atualiza a posição inteira usando essa velocidade do meio
        R = R + dt * V_meio
        
        ! 4. Recalcula a aceleração usando a NOVA posição que acabamos de achar
        call calcula_aceleracao(R, m1, m2, G, A)
        
        ! 5. Atualiza a velocidade inteira
        V = V_meio + 0.5d0 * dt * A
        
        ! Atualiza o tempo
        t = t + dt
    end do
    
    close(20)
end subroutine metodo_leapfrog


! ==============================================================================
! SUBROTINAS AUXILIARES DE FÍSICA
! ==============================================================================

! Calcula a aceleração da gravidade usando as posições atuais
subroutine calcula_aceleracao(R, m1, m2, G, A)
    implicit none
    real(8), intent(in)  :: R(4), m1, m2, G
    real(8), intent(out) :: A(4)
    real(8) :: dx, dy, dist3
    
    ! Distância no eixo x e y (r2 - r1)
    dx = R(3) - R(1)
    dy = R(4) - R(2)
    
    ! Distância total ao cubo (necessária para a fórmula da força)
    dist3 = (sqrt(dx**2 + dy**2))**3
    
    ! Acelerações (a = F/m)
    ! Corpo 1
    A(1) =  G * m2 * dx / dist3
    A(2) =  G * m2 * dy / dist3
    
    ! Corpo 2
    A(3) = -G * m1 * dx / dist3
    A(4) = -G * m1 * dy / dist3
end subroutine calcula_aceleracao

! Calcula a energia mecânica total (Cinética + Potencial)
subroutine calcula_energia(R, V, m1, m2, G, E)
    implicit none
    real(8), intent(in)  :: R(4), V(4), m1, m2, G
    real(8), intent(out) :: E
    real(8) :: K, U, dist, dx, dy
    
    ! Energia Cinética (K = mv^2 / 2)
    K = 0.5d0 * m1 * (V(1)**2 + V(2)**2) + &
        0.5d0 * m2 * (V(3)**2 + V(4)**2)
        
    ! Energia Potencial Gravitacional (U = -G*m1*m2 / dist)
    dx = R(3) - R(1)
    dy = R(4) - R(2)
    dist = sqrt(dx**2 + dy**2)
    U = -G * m1 * m2 / dist
    
    ! Energia Total
    E = K + U
end subroutine calcula_energia