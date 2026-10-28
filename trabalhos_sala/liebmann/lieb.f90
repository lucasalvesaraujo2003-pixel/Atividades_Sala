! ==============================================================================
! RESUMO DO MÉTODO:
! Este algoritmo implementa o método iterativo de Liebmann (também conhecido 
! como método de Gauss-Seidel aplicado a malhas espaciais) para a resolução 
! numérica da Equação de Laplace bidimensional utilizando o Método das 
! Diferenças Finitas. O código atualiza iterativamente o valor da temperatura 
! T(i,j) em cada nó interno da malha com base nos valores dos nós vizinhos, 
! repetindo o processo até que a maior variação encontrada (erro_max) seja 
! inferior à tolerância estipulada (tol) ou até atingir max_iter.
!
! APLICAÇÃO GERAL:
! Resolução de Equações Diferenciais Parciais (EDPs) elípticas. No contexto 
! físico do código, ele soluciona um problema de transferência de calor 
! bidimensional em regime permanente (estacionário) aplicado a uma aleta. 
! Na engenharia, essa abordagem é amplamente utilizada para mapear o perfil 
! térmico em seções transversais, analisar o desempenho de dissipadores de 
! calor e verificar gradientes de temperatura em componentes estruturais.
!
! INPUTS (Entradas):
! - Geometria do domínio: L e H.
! - Malha de discretização: nx e ny.
! - Condições térmicas: Tb (temperatura da base) e T_inf (temperatura ambiente).
! - Parâmetros numéricos de parada: tol e max_iter.
!
! OUTPUTS (Saídas):
! - Impressão no terminal do resumo da simulação contendo: iter, erro_max 
!   e tempo_total.
! - Geração e exportação de dois arquivos de dados:
!   'saida.dat'  : Mapeamento completo contendo x, y e T de toda a malha.
!   'centro.dat' : Perfil linear de temperatura contendo x e T restritos  
!                  apenas à linha central da aleta.
! ==============================================================================

! Escrever a matriz de temperatura na tela:
! do j = 1, ny
!    write(*,'(*(F12.2))') (T(i,j), i=1, nx)
! end do

program lieb
    implicit none
    integer :: i, j, iter
    ! Parâmetros dados pelo usuário
    real(8) :: L, H, Tb, T_inf, tol
    integer :: max_iter, nx, ny

    ! Construção da malha
    real(8) :: delta_x, delta_y
    real(8), allocatable :: T(:,:), x(:), y(:)
    real(8) :: erro_max, T_old, diferenca, j_centro

    ! Discretização de Laplace
    real(8) :: beta_sq, denominador

    ! Variáveis de tempo
    real(8) :: tempo_inicial, tempo_final, tempo_total
    
    ! print *, 'Transferência de Calor em uma aleta'
    ! print *, 'Digite o comprimento da aleta L: '
    ! read *, L
    ! print *, 'Digite a altura da aleta H: '
    ! read *, H
    ! print *, 'Digite o número de nós na direção x: '
    ! read *, nx
    ! print *, 'Digite o número de nós na direção y: '
    ! read *, ny
    ! print *, 'Digite a temperatura da base (parede esquerda) Tb: '
    ! read *, Tb
    ! print *, 'Digite a temperatura ambiente T_inf: '
    ! read *, T_inf
    ! print *, 'Digite a tolerância de convergência: '
    ! read *, tol
    ! print *, 'Digite o número máximo de iterações: '
    ! read *, max_iter

    L = 40.0d0
    H = 40.0d0
    Tb = 200.0d0
    T_inf = 25.0d0
    tol = 1.0d-6
    max_iter = 100000
    nx = 81
    ny = 81

    ! Criação da malha
    delta_x = L/real(nx - 1)
    delta_y = H/real(ny - 1)

    allocate(T(nx,ny))
    allocate(x(nx))
    allocate(y(ny))
    
    do i = 1, nx
        x(i) = real(i-1)*delta_x
    end do
    do j = 1, ny
        y(j) = real(j-1)*delta_y
    end do

    ! Aplicação das condições de contorno e dos nós internos
    T = T_inf
    do j = 1, ny
        T(1,j) = Tb
    end do

    beta_sq = (delta_x/delta_y)**2
    denominador = 2.0d0*(1.0d0 + beta_sq)
    iter = 0
    erro_max = tol + 1.0d0

    call cpu_time(tempo_inicial)

    do while (erro_max.gt.tol.and.iter.lt.max_iter)
        iter = iter + 1
        erro_max = 0.0d0
        do i = 2, nx-1
            do j = 2, ny-1
                T_old = T(i,j)
                T(i,j) = (T(i,j+1)+T(i,j-1)+beta_sq*(T(i+1,j)+T(i-1,j)))/denominador
                diferenca = abs(T(i,j) - T_old)
                if (diferenca.gt.erro_max) then
                    erro_max = diferenca
                end if
            end do
        end do
    end do

    call cpu_time(tempo_final)
    tempo_total = tempo_final - tempo_inicial

    write(*,*) ''
    write(*,*) '================================================='
    write(*,*) '               RESUMO DA SIMULACAO               '
    write(*,*) '================================================='
    write(*,'(A, I14)')       ' Numero de iteracoes  : ', iter
    write(*,'(A, E14.6)')    ' Erro final iterativo : ', erro_max
    write(*,'(A, F14.4, A)') ' Tempo computacional  : ', tempo_total, ' s'
    write(*,*) '================================================='
    write(*,*) ''

    j_centro = (ny + 1)/2

    open(10, file='saida.dat', status='replace')
    open(11, file='centro.dat', status='replace')
    do i = 1, nx
        do j = 1, ny
            write(10, '(3(F14.6))') x(i), y(j), T(i,j)
            if (j == j_centro) write(11, '(2(F14.6))') x(i), T(i,j)
        end do
        write(10, *)
    end do
    close(10)
    close(11)

end program lieb
