# PPC1 — Método de Runge-Kutta de 4ª Ordem

Código em Fortran para a solução numérica do problema de sedimentação de uma partícula esférica. Os comentários no próprio código identificam as variáveis físicas e as etapas das análises pedidas na PPC1.

```fortran
program ppc1_rk4
    implicit none

    integer :: i, j
    real(8), parameter :: tf = 8.0d0
    real(8), dimension(3) :: St_list
    real(8), dimension(4) :: h_list
    real(8), dimension(6) :: Re_list
    real(8) :: St, Re, h, erro

    ! Valores usados nas analises da atividade
    St_list = (/ 0.2d0, 1.0d0, 5.0d0 /)
    h_list  = (/ 0.4d0, 0.2d0, 0.1d0, 0.05d0 /)
    Re_list = (/ 0.0d0, 0.1d0, 0.25d0, 0.5d0, 0.75d0, 1.0d0 /)

    ! -------------------------------------------------------------
    ! 1) Re -> 0: compara RK4 com a solucao analitica para varios St
    ! -------------------------------------------------------------
    do i = 1, size(St_list)
        St = St_list(i)
        Re = 0.0d0
        h  = 0.01d0
        call simular(St, Re, h, tf, i, 'stokes')
    end do

    ! -------------------------------------------------------------
    ! 2) Estudo do passo de tempo h para St = 1 e Re = 0
    ! -------------------------------------------------------------
    open(20, file='erro_h.dat', status='replace')
    write(20,'(A)') '# h     erro_maximo'

    do i = 1, size(h_list)
        h = h_list(i)
        call erro_maximo(1.0d0, 0.0d0, h, 4.0d0, erro)
        write(20,'(2ES20.10)') h, erro
    end do
    close(20)

    ! -------------------------------------------------------------
    ! 3) e 4) Re != 0: solucao numerica e validacao analitica
    ! Exemplo: St = 1 e Re = 0.5
    ! -------------------------------------------------------------
    call simular(1.0d0, 0.5d0, 0.01d0, tf, 1, 'validacao_re')

    ! -------------------------------------------------------------
    ! 5) Varia Re e observa o desvio em relacao ao regime de Stokes
    ! -------------------------------------------------------------
    do j = 1, size(Re_list)
        Re = Re_list(j)
        call simular(1.0d0, Re, 0.01d0, tf, j, 'reynolds')
    end do

    print *, 'Execucao concluida.'
    print *, 'Arquivos .dat gerados na pasta atual.'

contains

    ! EDO adimensional:
    ! dv/dt = [1 - v - (3/8) Re v^2] / St
    ! v  : velocidade adimensional da particula
    ! t  : tempo adimensional
    ! St : numero de Stokes
    ! Re : numero de Reynolds da particula
    real(8) function f(v, St, Re)
        implicit none
        real(8), intent(in) :: v, St, Re
        f = (1.0d0 - v - (3.0d0/8.0d0)*Re*v*v)/St
    end function f


    ! Realiza um passo do metodo de Runge-Kutta de quarta ordem.
    subroutine passo_rk4(v, St, Re, h, vnovo)
        implicit none
        real(8), intent(in)  :: v, St, Re, h
        real(8), intent(out) :: vnovo
        real(8) :: k1, k2, k3, k4

        k1 = f(v, St, Re)
        k2 = f(v + 0.5d0*h*k1, St, Re)
        k3 = f(v + 0.5d0*h*k2, St, Re)
        k4 = f(v + h*k3, St, Re)

        vnovo = v + (h/6.0d0)*(k1 + 2.0d0*k2 + 2.0d0*k3 + k4)
    end subroutine passo_rk4


    ! Solucao analitica usada para validar o RK4.
    ! Para Re = 0 recupera v = 1 - exp(-t/St).
    real(8) function v_analitica(t, St, Re)
        implicit none
        real(8), intent(in) :: t, St, Re
        real(8) :: S, vinf, q, e

        if (abs(Re) < 1.0d-12) then
            v_analitica = 1.0d0 - exp(-t/St)
        else
            S    = sqrt(1.0d0 + 1.5d0*Re)
            vinf = 2.0d0/(1.0d0 + S)
            q    = (S - 1.0d0)/(S + 1.0d0)
            e    = exp(-S*t/St)
            v_analitica = vinf*(1.0d0 - e)/(1.0d0 + q*e)
        end if
    end function v_analitica


    ! Executa uma simulacao e salva:
    ! tempo, solucao numerica, solucao analitica e erro absoluto.
    subroutine simular(St, Re, h, tfinal, indice, prefixo)
        implicit none
        real(8), intent(in) :: St, Re, h, tfinal
        integer, intent(in) :: indice
        character(len=*), intent(in) :: prefixo

        integer :: unidade
        real(8) :: t, v, vnovo, va, erro_abs, hi
        character(len=100) :: nome

        write(nome,'(A,"_",I0,".dat")') trim(prefixo), indice
        unidade = 30 + indice

        open(unidade, file=trim(nome), status='replace')
        write(unidade,'(A)') '# t   v_numerica   v_analitica   erro_abs'

        t = 0.0d0
        v = 0.0d0

        do while (t < tfinal - 1.0d-14)
            va = v_analitica(t, St, Re)
            erro_abs = abs(v - va)
            write(unidade,'(4ES20.10)') t, v, va, erro_abs

            hi = min(h, tfinal - t)
            call passo_rk4(v, St, Re, hi, vnovo)

            v = vnovo
            t = t + hi
        end do

        va = v_analitica(t, St, Re)
        erro_abs = abs(v - va)
        write(unidade,'(4ES20.10)') t, v, va, erro_abs

        close(unidade)
    end subroutine simular


    ! Calcula o maior erro absoluto da simulacao.
    subroutine erro_maximo(St, Re, h, tfinal, erro)
        implicit none
        real(8), intent(in)  :: St, Re, h, tfinal
        real(8), intent(out) :: erro

        real(8) :: t, v, vnovo, va, hi

        t = 0.0d0
        v = 0.0d0
        erro = 0.0d0

        do while (t < tfinal - 1.0d-14)
            va = v_analitica(t, St, Re)
            erro = max(erro, abs(v - va))

            hi = min(h, tfinal - t)
            call passo_rk4(v, St, Re, hi, vnovo)

            v = vnovo
            t = t + hi
        end do

        va = v_analitica(t, St, Re)
        erro = max(erro, abs(v - va))
    end subroutine erro_maximo

end program ppc1_rk4

```
