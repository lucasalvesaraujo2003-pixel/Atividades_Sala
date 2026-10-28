! ==============================================================================
! RESUMO DO MÉTODO:
! Este algoritmo implementa os métodos numéricos de Euler (explícito) e 
! Heun (Euler Aprimorado / Preditor-Corretor) para a resolução de Problemas 
! de Valor Inicial (PVI). O código resolve numericamente uma Equação Diferencial 
! Ordinária (EDO) de segunda ordem, tratada no código como um sistema de duas 
! EDOs de primeira ordem. O método de Heun utiliza a inclinação 
! inicial para prever uma solução parcial e, em seguida, corrige esse valor 
! utilizando a média das inclinações nos extremos do intervalo.
!
! APLICAÇÃO GERAL:
! Integração numérica de sistemas dinâmicos no domínio do tempo. Pela função 
! de aceleracao definida internamente, o código simula especificamente 
! um oscilador harmônico amortecido sujeito a uma excitação externa harmônica 
! (o clássico sistema massa-mola-amortecedor). Na engenharia, essa formulação 
! é amplamente empregada para análise de vibrações mecânicas, resposta de 
! suspensões automotivas e análise de transientes estruturais.
!
! INPUTS (Entradas):
! - Parâmetros físicos (constantes): m, k, c, F0 e omega.
! - Parâmetros numéricos: dt e t_max.
!
! OUTPUTS (Saídas):
! - O script gera e exporta dois arquivos de dados em disco ('euler.dat' 
!   e 'heun.dat').
! - Para o método de Euler, as variáveis salvas são: t_euler, x_euler e v_euler.
! - Para o método de Heun, as variáveis salvas são: t_heun, x_heun e v_heun.
! ==============================================================================

program euler_heun
    implicit none
    integer :: i
    real(8), parameter :: m = 1.0d0, k = 4.0d0, c = 0.5d0, F0 = 1.0d0, omega = 2.0d0
    real(8), parameter :: dt = 0.1d0, t_max = 10.0d0
    integer, parameter :: passo = int(t_max/dt) + 1
    real(8) :: t_euler(passo), x_euler(passo), v_euler(passo)
    real(8) :: t_heun(passo), x_heun(passo), v_heun(passo)
    real(8) :: x_prev, v_prev, a_prev, a_atual

    ! Método de Euler
    open(unit=10, file='euler.dat', status='replace')

    x_euler(1) = 0.0d0
    t_euler(1) = 0.0d0
    v_euler(1) = 0.0d0

    do i=2, passo
        x_euler(i) = x_euler(i-1) + dt*v_euler(i-1)
        v_euler(i) = v_euler(i-1) + dt*aceleracao(x_euler(i-1),v_euler(i-1),t_euler(i-1))
        t_euler(i) = t_euler(i-1) + dt
    end do

    do i=1, passo
        write(10,*) t_euler(i), x_euler(i), v_euler(i)
    end do

    close(10)

    ! Método de Heun
    open(unit=20, file='heun.dat', status='replace')

    x_heun(1) = 0.0d0
    t_heun(1) = 0.0d0
    v_heun(1) = 0.0d0
    
    do i=2, passo
        a_atual = aceleracao(x_heun(i-1),v_heun(i-1),t_heun(i-1))
        x_prev = x_heun(i-1) + dt*v_heun(i-1)
        v_prev = v_heun(i-1) + dt*a_atual
        a_prev = aceleracao(x_prev,v_prev,t_heun(i-1)+dt)

        x_heun(i) = x_heun(i-1) + dt*(v_heun(i-1) + v_prev)/2
        v_heun(i) = v_heun(i-1) + dt*(a_prev + a_atual)/2
        t_heun(i) = t_heun(i-1) + dt
    end do

    do i=1, passo
        write(20,*) t_heun(i), x_heun(i), v_heun(i)
    end do

    close(20)
    
    contains
        real(8) function aceleracao(pos, vel, tempo)
            real(8), intent(in) :: pos, vel, tempo
            aceleracao = (- k*pos - c*vel + F0*sin(omega*tempo)) / m
        end function aceleracao
end program euler_heun
