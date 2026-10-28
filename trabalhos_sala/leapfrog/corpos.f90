! Resolver por Euler e dps resolver pelo Leapfrog
! arrumar pras normas dos vetores nas fórmulas de energia
program leap
    implicit none
    real(8) :: G, m1, m2, dt, t_max, E, tempo
    real(8), dimension(2) :: F12, r1_old, r2_old, v1_old, v2_old
    real(8), dimension(2) :: r1_new, r2_new, v1_new, v2_new
    integer :: i, n_passos

    G = 1.0d0
    m1 = 1.0d0
    m2 = 1.0d0
    dt = 1d-3
    t_max = 100.0d0
    n_passos = int(t_max/dt)
    tempo = 0.0d0

    r1_old = [-0.5d0, 0.0d0]
    r2_old = [0.5d0, 0.0d0]
    v1_old = [0.0d0, -0.5d0]
    v2_old = [0.0d0, 0.5d0]

    F12(1) = forca(r1_old(1), r2_old(1))
    F12(2) = forca(r1_old(2), r2_old(2))

    ! Método de Euler

    contains
        real(8) function forca(pos1, pos2)
            real(8), intent(in) :: pos1, pos2
            forca = G*m1*m2*(pos2-pos1)/(abs(pos2-pos1)**3)
        end function forca

        real(8) function dvdt(massa, pos1, pos2)
            real(8), intent(in) :: massa, pos1, pos2
            ! Para encontrar dv1, deve-se utilizar m2
            ! Para encontrar dv2, deve-se utilizar m1
            if (massa.eq.m2) then
                dvdt = G*massa*(pos2-pos1)/(abs(pos2-pos1)**3)
            else
                dvdt = - G*massa*(pos2-pos1)/(abs(pos2-pos1)**3)
            end if
        end function dvdt

        real(8) function en_cin(vel1, vel2)
            real(8), intent(in) :: vel1, vel2
            en_cin = m1*(abs(vel1)**2)/2 + m2*(abs(vel2)**2)/2
        end function en_cin

        real(8) function en_pot(pos1, pos2)
            real(8), intent(in) :: pos1, pos2
            en_pot = - G*m1*m2/(abs(pos2-pos1))
        end function en_pot

        subroutine euler
            open(unit=10, file='euler.dat', status='replace')
            write(10, *) tempo, r1_old(1), r1_old(2), v1_old(1), v1_old(2), r2_old(1), r2_old(2), v2_old(1), v2_old(2), E
            do i=1, n_passos
                tempo = tempo + dt
                
                r1_new(1) = r1_old(1) + dt*v1_old(1)
                r1_new(2) = r1_old(2) + dt*v1_old(2)

                r2_new(1) = r2_old(1) + dt*v2_old(1)
                r2_new(2) = r2_old(2) + dt*v2_old(2)

                v1_new(1) = v1_old(1) + dt*dvdt(m2, r1_old(1), r2_old(1))
                v1_new(2) = v1_old(2) + dt*dvdt(m2, r1_old(2), r2_old(2))

                v2_new(1) = v2_old(1) + dt*dvdt(m1, r1_old(1), r2_old(1))
                v2_new(2) = v2_old(2) + dt*dvdt(m1, r1_old(2), r2_old(2))
            end do
        end subroutine euler
end program leap