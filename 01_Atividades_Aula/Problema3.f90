program cna_toolkit_mestre
  ! =========================================================================
  ! TOOLKIT MESTRE DE CÁLCULO NUMÉRICO APLICADO (CNA) - FORTRAN 90
  ! Estrutura Genérica e Modular para a Prova Prática
  ! -> Não altere os módulos matemáticos!
  ! -> Altere apenas as funções no final do arquivo (BLOCO 5).
  ! =========================================================================
  implicit none

  ! Variáveis para testes dos módulos
  real*8 :: raiz_bissecao, raiz_newton
  complex*16 :: raiz_muller
  
  integer, parameter :: N_sys = 3
  real*8, dimension(N_sys, N_sys) :: A_teste, L_mat, U_mat
  real*8, dimension(N_sys) :: b_teste, x_gauss, x_lu
  
  integer, parameter :: N_tri = 4
  real*8, dimension(N_tri) :: e_diag, f_diag, g_diag, r_tri, x_thomas
  
  real*8 :: y_final_rk4
  
  real*8 :: x_otimo_aurea, custo_minimo_aurea
  
  integer :: i

  print *, '================================================================='
  print *, '   BIBLIOTECA GENÉRICA DE MÉTODOS NUMÉRICOS - CNA (FORTRAN)      '
  print *, '================================================================='
  print *, ''

  ! -------------------------------------------------------------------------
  ! TESTE 1: ZEROS DE FUNÇÕES
  ! -------------------------------------------------------------------------
  print *, '--- MÓDULO 1: ZEROS DE FUNÇÕES ---'
  call bissecao(1.0d0, 3.0d0, 1.0d-6, 100, raiz_bissecao)
  print '(A, F10.6)', 'Bissecao       -> Raiz x = ', raiz_bissecao
  
  call newton_raphson(2.0d0, 1.0d-6, 100, raiz_newton)
  print '(A, F10.6)', 'Newton-Raphson -> Raiz x = ', raiz_newton
  
  call muller(0.0d0, 1.0d0, 2.0d0, 1.0d-6, 100, raiz_muller)
  print '(A, F10.6, A, F10.6, A)', 'Muller         -> Raiz x = ', &
        real(raiz_muller), ' + (', aimag(raiz_muller), ')i'
  print *, ''

  ! -------------------------------------------------------------------------
  ! TESTE 2: SISTEMAS LINEARES & DECOMPOSIÇÃO LU
  ! -------------------------------------------------------------------------
  print *, '--- MÓDULO 2: SISTEMAS LINEARES E LU ---'
  A_teste(1,:) = (/  3.0d0, -0.1d0, -0.2d0 /)
  A_teste(2,:) = (/  0.1d0,  7.0d0, -0.3d0 /)
  A_teste(3,:) = (/  0.3d0, -0.2d0, 10.0d0 /)
  b_teste = (/ 7.85d0, -19.30d0, 71.40d0 /)
  
  call gauss_pivoteamento(A_teste, b_teste, N_sys, x_gauss)
  print '(A, 3F10.4)', 'Gauss c/ Pivoteamento -> x = ', x_gauss
  
  call decomposicao_lu(A_teste, b_teste, N_sys, x_lu, L_mat, U_mat)
  print '(A, 3F10.4)', 'Decomposicao LU       -> x = ', x_lu
  print *, ''

  ! -------------------------------------------------------------------------
  ! TESTE 3: ALGORITMO DE THOMAS (MATRIZ TRIDIAGONAL)
  ! -------------------------------------------------------------------------
  print *, '--- MÓDULO 3: ALGORITMO DE THOMAS (TDMA) ---'
  f_diag = (/ 2.04d0, 2.04d0, 2.04d0, 2.04d0 /)
  e_diag = (/ 0.00d0,-1.00d0,-1.00d0,-1.00d0 /)
  g_diag = (/-1.00d0,-1.00d0,-1.00d0, 0.00d0 /)
  r_tri  = (/ 40.8d0,  0.8d0,  0.8d0,200.8d0 /)
  
  call thomas_tdma(e_diag, f_diag, g_diag, r_tri, N_tri, x_thomas)
  print '(A, 4F10.4)', 'Thomas (TDMA) -> x = ', x_thomas
  print *, ''

  ! -------------------------------------------------------------------------
  ! TESTE 4: EDO (RUNGE-KUTTA 4a ORDEM)
  ! -------------------------------------------------------------------------
  print *, '--- MÓDULO 4: INTEGRAÇÃO DE EDO (RK4) ---'
  call rk4_edo(0.0d0, 2.0d0, 1.0d0, 0.1d0, y_final_rk4)
  print '(A, F10.6)', 'Runge-Kutta 4 (t=2.0) -> y = ', y_final_rk4
  print *, ''

  ! -------------------------------------------------------------------------
  ! TESTE 5: OTIMIZAÇÃO (RAZÃO ÁUREA)
  ! -------------------------------------------------------------------------
  print *, '--- MÓDULO 5: OTIMIZAÇÃO POR RAZÃO ÁUREA ---'
  call razao_aurea(0.0d0, 4.0d0, 1.0d-6, 100, x_otimo_aurea, custo_minimo_aurea)
  print '(A, F10.6, A, F10.6)', 'Razao Aurea -> x* = ', x_otimo_aurea, ' | f(x*) = ', custo_minimo_aurea
  print *, '================================================================='

contains

  ! =========================================================================
  ! BLOCO 1: MÓDULOS DE ZEROS DE FUNÇÕES
  ! =========================================================================
  subroutine bissecao(a_in, b_in, tol, max_iter, x_raiz)
    implicit none
    real*8, intent(in) :: a_in, b_in, tol
    integer, intent(in) :: max_iter
    real*8, intent(out) :: x_raiz
    real*8 :: a, b, c, fa, fc
    integer :: iter

    a = a_in; b = b_in
    do iter = 1, max_iter
       c = (a + b) / 2.0d0
       fc = f_raiz(c)
       fa = f_raiz(a)
       
       if (abs(fc) < tol .or. (b - a) / 2.0d0 < tol) exit
       
       if (fa * fc < 0.0d0) then
          b = c
       else
          a = c
       end if
    end do
    x_raiz = c
  end subroutine bissecao

  subroutine newton_raphson(x0, tol, max_iter, x_raiz)
    implicit none
    real*8, intent(in) :: x0, tol
    integer, intent(in) :: max_iter
    real*8, intent(out) :: x_raiz
    real*8 :: x, dx, fx, dfx
    integer :: iter

    x = x0
    do iter = 1, max_iter
       fx = f_raiz(x)
       dfx = f_derivada_raiz(x)
       
       if (abs(dfx) < 1.0d-12) exit
       
       dx = fx / dfx
       x = x - dx
       if (abs(dx) < tol) exit
    end do
    x_raiz = x
  end subroutine newton_raphson

  subroutine muller(x0_in, x1_in, x2_in, tol, max_iter, x_raiz)
    implicit none
    real*8, intent(in) :: x0_in, x1_in, x2_in, tol
    integer, intent(in) :: max_iter
    complex*16, intent(out) :: x_raiz
    complex*16 :: x0, x1, x2, h0, h1, d0, d1, a, b, c, rad, den, dx, f0, f1, f2
    integer :: iter

    x0 = cmplx(x0_in, 0.0d0, 8); x1 = cmplx(x1_in, 0.0d0, 8); x2 = cmplx(x2_in, 0.0d0, 8)

    do iter = 1, max_iter
       f0 = cmplx(f_raiz(real(x0)), 0.0d0, 8)
       f1 = cmplx(f_raiz(real(x1)), 0.0d0, 8)
       f2 = cmplx(f_raiz(real(x2)), 0.0d0, 8)

       h0 = x1 - x0; h1 = x2 - x1
       d0 = (f1 - f0) / h0; d1 = (f2 - f1) / h1
       
       a = (d1 - d0) / (h1 + h0); b = a * h1 + d1; c = f2
       rad = sqrt(b**2 - 4.0d0 * a * c)
       
       if (abs(b + rad) > abs(b - rad)) then
          den = b + rad
       else
          den = b - rad
       end if

       dx = -2.0d0 * c / den
       x2 = x2 + dx
       if (abs(dx) < tol) exit
       x0 = x1; x1 = x2 - dx
    end do
    x_raiz = x2
  end subroutine muller

  ! =========================================================================
  ! BLOCO 2: MÓDULOS DE SISTEMAS LINEARES E MATRIZES
  ! =========================================================================
  subroutine gauss_pivoteamento(A_in, b_in, n, x_sol)
    implicit none
    integer, intent(in) :: n
    real*8, dimension(n,n), intent(in) :: A_in
    real*8, dimension(n), intent(in) :: b_in
    real*8, dimension(n), intent(out) :: x_sol
    real*8, dimension(n,n) :: A
    real*8, dimension(n) :: b
    real*8 :: fator, temp, soma
    integer :: i, j, k, piv_lin

    A = A_in; b = b_in

    do k = 1, n - 1
       piv_lin = k
       do i = k + 1, n
          if (abs(A(i,k)) > abs(A(piv_lin,k))) piv_lin = i
       end do

       if (piv_lin /= k) then
          do j = 1, n
             temp = A(k,j); A(k,j) = A(piv_lin,j); A(piv_lin,j) = temp
          end do
          temp = b(k); b(k) = b(piv_lin); b(piv_lin) = temp
       end if

       do i = k + 1, n
          fator = A(i,k) / A(k,k)
          do j = k + 1, n
             A(i,j) = A(i,j) - fator * A(k,j)
          end do
          b(i) = b(i) - fator * b(k)
          A(i,k) = 0.0d0
       end do
    end do

    x_sol(n) = b(n) / A(n,n)
    do i = n - 1, 1, -1
       soma = b(i)
       do j = i + 1, n
          soma = soma - A(i,j) * x_sol(j)
       end do
       x_sol(i) = soma / A(i,i)
    end do
  end subroutine gauss_pivoteamento

  subroutine decomposicao_lu(A_in, b_in, n, x_sol, L_out, U_out)
    implicit none
    integer, intent(in) :: n
    real*8, dimension(n,n), intent(in) :: A_in
    real*8, dimension(n), intent(in) :: b_in
    real*8, dimension(n), intent(out) :: x_sol
    real*8, dimension(n,n), intent(out) :: L_out, U_out
    real*8, dimension(n) :: y
    real*8 :: soma
    integer :: i, j, k

    L_out = 0.0d0; U_out = 0.0d0
    
    do k = 1, n
       L_out(k,k) = 1.0d0
       do j = k, n
          soma = 0.0d0
          do i = 1, k - 1
             soma = soma + L_out(k,i) * U_out(i,j)
          end do
          U_out(k,j) = A_in(k,j) - soma
       end do
       
       do i = k + 1, n
          soma = 0.0d0
          do j = 1, k - 1
             soma = soma + L_out(i,j) * U_out(j,k)
          end do
          L_out(i,k) = (A_in(i,k) - soma) / U_out(k,k)
       end do
    end do

    do i = 1, n
       soma = b_in(i)
       do j = 1, i - 1
          soma = soma - L_out(i,j) * y(j)
       end do
       y(i) = soma / L_out(i,i)
    end do

    x_sol(n) = y(n) / U_out(n,n)
    do i = n - 1, 1, -1
       soma = y(i)
       do j = i + 1, n
          soma = soma - U_out(i,j) * x_sol(j)
       end do
       x_sol(i) = soma / U_out(i,i)
    end do
  end subroutine decomposicao_lu

  subroutine thomas_tdma(e_in, f_in, g_in, r_in, n, x_sol)
    implicit none
    integer, intent(in) :: n
    real*8, dimension(n), intent(in) :: e_in, f_in, g_in, r_in
    real*8, dimension(n), intent(out) :: x_sol
    real*8, dimension(n) :: e, f, g, r
    integer :: k

    e = e_in; f = f_in; g = g_in; r = r_in

    do k = 2, n
       e(k) = e(k) / f(k-1)
       f(k) = f(k) - e(k) * g(k-1)
       r(k) = r(k) - e(k) * r(k-1)
    end do

    x_sol(n) = r(n) / f(n)
    do k = n - 1, 1, -1
       x_sol(k) = (r(k) - g(k) * x_sol(k+1)) / f(k)
    end do
  end subroutine thomas_tdma

  ! =========================================================================
  ! BLOCO 3: MÓDULO DE EQUAÇÃO DIFERENCIAL ORDINÁRIA (EDO)
  ! =========================================================================
  subroutine rk4_edo(t0, tf, y0, h, y_final)
    implicit none
    real*8, intent(in) :: t0, tf, y0, h
    real*8, intent(out) :: y_final
    real*8 :: t, y, k1, k2, k3, k4
    integer :: passos, i

    passos = int((tf - t0) / h)
    t = t0; y = y0

    do i = 1, passos
       k1 = func_edo(t, y)
       k2 = func_edo(t + 0.5d0 * h, y + 0.5d0 * k1 * h)
       k3 = func_edo(t + 0.5d0 * h, y + 0.5d0 * k2 * h)
       k4 = func_edo(t + h, y + k3 * h)

       y = y + (h / 6.0d0) * (k1 + 2.0d0 * k2 + 2.0d0 * k3 + k4)
       t = t + h
    end do
    y_final = y
  end subroutine rk4_edo

  ! =========================================================================
  ! BLOCO 4: MÓDULO DE OTIMIZAÇÃO UNIDIMENSIONAL
  ! =========================================================================
  subroutine razao_aurea(a_in, b_in, tol, max_iter, x_opt, f_opt)
    implicit none
    real*8, intent(in) :: a_in, b_in, tol
    integer, intent(in) :: max_iter
    real*8, intent(out) :: x_opt, f_opt
    real*8 :: a, b, tau, x1, x2, f1, f2
    integer :: iter

    a = a_in; b = b_in
    tau = (sqrt(5.0d0) - 1.0d0) / 2.0d0

    x1 = b - tau * (b - a)
    x2 = a + tau * (b - a)
    f1 = func_objetivo_otimizacao(x1)
    f2 = func_objetivo_otimizacao(x2)

    do iter = 1, max_iter
       if (abs(b - a) < tol) exit
       
       if (f1 > f2) then
          a = x1; x1 = x2; f1 = f2
          x2 = a + tau * (b - a)
          f2 = func_objetivo_otimizacao(x2)
       else
          b = x2; x2 = x1; f2 = f1
          x1 = b - tau * (b - a)
          f1 = func_objetivo_otimizacao(x1)
       end if
    end do

    x_opt = (a + b) / 2.0d0
    f_opt = func_objetivo_otimizacao(x_opt)
  end subroutine razao_aurea

  ! =========================================================================
  ! BLOCO 5: FUNÇÕES MATEMÁTICAS DA PROVA (ALTERE APENAS AQUI!)
  ! =========================================================================
  
  ! 5.1 -> Função para Busca de Raízes: f(x) = 0
  real*8 function f_raiz(x)
    implicit none
    real*8, intent(in) :: x
    ! Exemplo: x^2 - 4 = 0 (Raiz em x = 2)
    f_raiz = x**2 - 4.0d0
  end function f_raiz

  ! 5.2 -> Derivada da Função para Newton-Raphson: f'(x)
  real*8 function f_derivada_raiz(x)
    implicit none
    real*8, intent(in) :: x
    ! Derivada de x^2 - 4
    f_derivada_raiz = 2.0d0 * x
  end function f_derivada_raiz

  ! 5.3 -> Equação Diferencial Ordinária: dy/dt = f(t, y)
  real*8 function func_edo(t, y)
    implicit none
    real*8, intent(in) :: t, y
    ! Exemplo: dy/dt = -2y + t
    func_edo = -2.0d0 * y + t
  end function func_edo

  ! 5.4 -> Função Objetivo para Otimização (Minimizar Custo/Área/etc.)
  real*8 function func_objetivo_otimizacao(x)
    implicit none
    real*8, intent(in) :: x
    ! Exemplo: Parábola com mínimo em x = 2.5
    func_objetivo_otimizacao = (x - 2.5d0)**2 + 3.0d0
  end function func_objetivo_otimizacao

end program cna_toolkit_mestre
