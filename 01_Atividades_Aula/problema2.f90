program master_cna_toolkit
  ! =========================================================================
  ! TOOLKIT MESTRE DE CÁLCULO NUMÉRICO APLICADO (CNA)
  ! Estrutura Generica e Modular para Consulta na Prova Pratica
  ! =========================================================================
  implicit none

  ! =========================================================================
  ! VARIÁVEIS DE TESTE E DEMONSTRAÇÃO DOS MÓDULOS
  ! =========================================================================
  ! Variaveis para Zeros de Funcao
  real*8 :: raiz_bissecao, raiz_newton
  complex*16 :: raiz_muller
  
  ! Variaveis para Sistemas Lineares (N x N)
  integer, parameter :: N_sys = 3
  real*8, dimension(N_sys, N_sys) :: A_mat, A_lu
  real*8, dimension(N_sys) :: b_vec, x_sol, x_lu
  
  ! Variaveis para Sistema Tridiagonal (Algoritmo de Thomas)
  integer, parameter :: N_tri = 4
  real*8, dimension(N_tri) :: e_diag, f_diag, g_diag, r_tri, x_tri
  
  ! Variaveis para EDO (Runge-Kutta 4a Ordem)
  real*8 :: t0, tf, y0, h_passo, y_final
  
  ! Variaveis para Otimizacao 1D
  real*8 :: x_min, f_min
  
  integer :: i, j

  print *, '================================================================='
  print *, '   BIBLIOTECA GENÉRICA DE MÉTODOS NUMÉRICOS - CNA (FORTRAN 90)   '
  print *, '================================================================='
  print *, ''

  ! -------------------------------------------------------------------------
  ! TESTE DO MÓDULO 1: ZEROS DE FUNÇÕES
  ! -------------------------------------------------------------------------
  print *, '--- MÓDULO 1: ZEROS DE FUNÇÕES ---'
  call bissecao(1.0d0, 3.0d0, 1.0d-6, 100, raiz_bissecao)
  print '(A, F10.6)', 'Raiz por Bissecao     : x = ', raiz_bissecao
  
  call newton_raphson(2.0d0, 1.0d-6, 100, raiz_newton)
  print '(A, F10.6)', 'Raiz por Newton-Raphson: x = ', raiz_newton
  
  call muller(0.0d0, 1.0d0, 2.0d0, 1.0d-6, 100, raiz_muller)
  print '(A, F10.6, A, F10.6, A)', 'Raiz por Muller         : x = ', &
        real(raiz_muller), ' + (', aimag(raiz_muller), ')i'
  print *, ''

  ! -------------------------------------------------------------------------
  ! TESTE DO MÓDULO 2: SISTEMAS LINEARES & DECOMPOSIÇÃO LU
  ! -------------------------------------------------------------------------
  print *, '--- MÓDULO 2: SISTEMAS LINEARES E LU ---'
  ! Matriz de Teste 3x3:
  ! [ 3.0  -0.1  -0.2 ] [ x1 ]   [  7.85 ]
  ! [ 0.1   7.0  -0.3 ] [ x2 ] = [ -19.30 ]
  ! [ 0.3  -0.2  10.0 ] [ x3 ]   [  71.40 ]
  A_mat(1,:) = (/  3.0d0, -0.1d0, -0.2d0 /)
  A_mat(2,:) = (/  0.1d0,  7.0d0, -0.3d0 /)
  A_mat(3,:) = (/  0.3d0, -0.2d0, 10.0d0 /)
  b_vec = (/ 7.85d0, -19.30d0, 71.40d0 /)
  
  call gauss_pivoteamento(A_mat, b_vec, N_sys, x_sol)
  print '(A, 3F10.4)', 'Solucao por Gauss c/ Pivoteamento: x = ', x_sol
  
  call decomposicao_lu_solve(A_mat, b_vec, N_sys, x_lu)
  print '(A, 3F10.4)', 'Solucao por Decomposicao LU        : x = ', x_lu
  print *, ''

  ! -------------------------------------------------------------------------
  ! TESTE DO MÓDULO 3: ALGORITMO DE THOMAS (MATRIZ TRIDIAGONAL)
  ! -------------------------------------------------------------------------
  print *, '--- MÓDULO 3: ALGORITMO DE THOMAS (TDMA) ---'
  f_diag = (/ 2.04d0, 2.04d0, 2.04d0, 2.04d0 /) ! Diagonal Principal
  e_diag = (/ 0.00d0,-1.00d0,-1.00d0,-1.00d0 /) ! Subdiagonal Inferior
  g_diag = (/-1.00d0,-1.00d0,-1.00d0, 0.00d0 /) ! Superdiagonal Superior
  r_tri  = (/ 40.8d0,  0.8d0,  0.8d0,200.8d0 /) ! Termos Independentes
  
  call thomas_tdma(e_diag, f_diag, g_diag, r_tri, N_tri, x_tri)
  print '(A, 4F10.4)', 'Solucao Thomas (TDMA): x = ', x_tri
  print *, ''

  ! -------------------------------------------------------------------------
  ! TESTE DO MÓDULO 4: RESOLUÇÃO DE EDO (RUNGE-KUTTA 4a ORDEM)
  ! -------------------------------------------------------------------------
  print *, '--- MÓDULO 4: INTEGRAÇÃO DE EDO (RK4) ---'
  t0 = 0.0d0
  tf = 2.0d0
  y0 = 1.0d0
  h_passo = 0.1d0
  
  call rk4_edo(t0, tf, y0, h_passo, y_final)
  print '(A, F10.6, A, F10.6)', 'Solucao RK4 em t = ', tf, ' : y = ', y_final
  print *, ''

  ! -------------------------------------------------------------------------
  ! TESTE DO MÓDULO 5: OTIMIZAÇÃO 1D (RAZÃO ÁUREA)
  ! -------------------------------------------------------------------------
  print *, '--- MÓDULO 5: OTIMIZAÇÃO POR RAZÃO ÁUREA ---'
  call razao_aurea(0.0d0, 4.0d0, 1.0d-6, 100, x_min, f_min)
  print '(A, F10.6, A, F10.6)', 'Ponto de Minimo: x* = ', x_min, ' | f(x*) = ', f_min
  print *, '================================================================='

contains

  ! =========================================================================
  ! MÓDULO 1: ZEROS DE FUNÇÕES (SUBROTINAS)
  ! =========================================================================

  ! --- 1.1: METODO DA BISSEÇÃO ---
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
       fc = f_zero(c)
       fa = f_zero(a)
       
       if (abs(fc) < tol .or. (b - a) / 2.0d0 < tol) exit
       
       if (fa * fc < 0.0d0) then
          b = c
       else
          a = c
       end if
    end do
    x_raiz = c
  end subroutine bissecao

  ! --- 1.2: METODO DE NEWTON-RAPHSON ---
  subroutine newton_raphson(x0, tol, max_iter, x_raiz)
    implicit none
    real*8, intent(in) :: x0, tol
    integer, intent(in) :: max_iter
    real*8, intent(out) :: x_raiz
    real*8 :: x, dx, fx, dfx
    integer :: iter

    x = x0
    do iter = 1, max_iter
       fx = f_zero(x)
       dfx = f_derivada(x)
       
       if (abs(dfx) < 1.0d-12) then
          print *, 'Erro: Derivada proxima de zero em Newton-Raphson.'
          exit
       end if
       
       dx = fx / dfx
       x = x - dx
       if (abs(dx) < tol) exit
    end do
    x_raiz = x
  end subroutine newton_raphson

  ! --- 1.3: METODO DE MÜLLER (SUPORTA RAÍZES COMPLEXAS) ---
  subroutine muller(x0_in, x1_in, x2_in, tol, max_iter, x_raiz)
    implicit none
    real*8, intent(in) :: x0_in, x1_in, x2_in, tol
    integer, intent(in) :: max_iter
    complex*16, intent(out) :: x_raiz
    
    complex*16 :: x0, x1, x2, h0, h1, d0, d1, a, b, c, rad, den, dx
    complex*16 :: f0, f1, f2
    integer :: iter

    x0 = cmplx(x0_in, 0.0d0, 8)
    x1 = cmplx(x1_in, 0.0d0, 8)
    x2 = cmplx(x2_in, 0.0d0, 8)

    do iter = 1, max_iter
       f0 = cmplx(f_zero(real(x0)), 0.0d0, 8) ! Para funcoes reais simples
       f1 = cmplx(f_zero(real(x1)), 0.0d0, 8)
       f2 = cmplx(f_zero(real(x2)), 0.0d0, 8)

       h0 = x1 - x0; h1 = x2 - x1
       d0 = (f1 - f0) / h0; d1 = (f2 - f1) / h1
       
       a = (d1 - d0) / (h1 + h0)
       b = a * h1 + d1
       c = f2

       rad = sqrt(b**2 - 4.0d0 * a * c)
       if (abs(b + rad) > abs(b - rad)) then
          den = b + rad
       else
          den = b - rad
       end if

       dx = -2.0d0 * c / den
       x2 = x2 + dx

       if (abs(dx) < tol) exit
       x0 = x1; x1 = x2 - dx; x2 = x2
    end do
    x_raiz = x2
  end subroutine muller

  ! =========================================================================
  ! MÓDULO 2: SISTEMAS LINEARES & MATRIZES
  ! =========================================================================

  ! --- 2.1: ELIMINAÇÃO GAUSSIANA COM PIVOTEAMENTO PARCIAL ---
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

    ! Eliminação Progressiva com Pivoteamento
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

    ! Substituição Regressiva
    x_sol(n) = b(n) / A(n,n)
    do i = n - 1, 1, -1
       soma = b(i)
       do j = i + 1, n
          soma = soma - A(i,j) * x_sol(j)
       end do
       x_sol(i) = soma / A(i,i)
    end do
  end subroutine gauss_pivoteamento

  ! --- 2.2: DECOMPOSIÇÃO LU (DOOLITTLE) E SOLUÇÃO DE SISTEMA ---
  subroutine decomposicao_lu_solve(A_in, b_in, n, x_sol)
    implicit none
    integer, intent(in) :: n
    real*8, dimension(n,n), intent(in) :: A_in
    real*8, dimension(n), intent(in) :: b_in
    real*8, dimension(n), intent(out) :: x_sol
    
    real*8, dimension(n,n) :: L, U
    real*8, dimension(n) :: y
    real*8 :: soma
    integer :: i, j, k

    L = 0.0d0; U = 0.0d0
    
    ! Fatoracao LU
    do k = 1, n
       L(k,k) = 1.0d0
       do j = k, n
          soma = 0.0d0
          do i = 1, k - 1
             soma = soma + L(k,i) * U(i,j)
          end do
          U(k,j) = A_in(k,j) - soma
       end do
       
       do i = k + 1, n
          soma = 0.0d0
          do j = 1, k - 1
             soma = soma + L(i,j) * U(j,k)
          end do
          L(i,k) = (A_in(i,k) - soma) / U(k,k)
       end do
    end do

    ! L * y = b (Substituicao Progressiva)
    do i = 1, n
       soma = b_in(i)
       do j = 1, i - 1
          soma = soma - L(i,j) * y(j)
       end do
       y(i) = soma / L(i,i)
    end do

    ! U * x = y (Substituicao Regressiva)
    x_sol(n) = y(n) / U(n,n)
    do i = n - 1, 1, -1
       soma = y(i)
       do j = i + 1, n
          soma = soma - U(i,j) * x_sol(j)
       end do
       x_sol(i) = soma / U(i,i)
    end do
  end subroutine decomposicao_lu_solve

  ! --- 2.3: ALGORITMO DE THOMAS (TDMA PARA MATRIZ TRIDIAGONAL) ---
  subroutine thomas_tdma(e_in, f_in, g_in, r_in, n, x_sol)
    implicit none
    integer, intent(in) :: n
    real*8, dimension(n), intent(in) :: e_in, f_in, g_in, r_in
    real*8, dimension(n), intent(out) :: x_sol
    
    real*8, dimension(n) :: e, f, g, r
    integer :: k

    e = e_in; f = f_in; g = g_in; r = r_in

    ! Etapa de Decomposição e Substituição Progressiva
    do k = 2, n
       e(k) = e(k) / f(k-1)
       f(k) = f(k) - e(k) * g(k-1)
       r(k) = r(k) - e(k) * r(k-1)
    end do

    ! Substituição Regressiva
    x_sol(n) = r(n) / f(n)
    do k = n - 1, 1, -1
       x_sol(k) = (r(k) - g(k) * x_sol(k+1)) / f(k)
    end do
  end subroutine thomas_tdma

  ! =========================================================================
  ! MÓDULO 3: INTEGRADOR DE EDO (RUNGE-KUTTA 4a ORDEM)
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
       k1 = f_edo(t, y)
       k2 = f_edo(t + 0.5d0 * h, y + 0.5d0 * k1 * h)
       k3 = f_edo(t + 0.5d0 * h, y + 0.5d0 * k2 * h)
       k4 = f_edo(t + h, y + k3 * h)

       y = y + (h / 6.0d0) * (k1 + 2.0d0 * k2 + 2.0d0 * k3 + k4)
       t = t + h
    end do
    y_final = y
  end subroutine rk4_edo

  ! =========================================================================
  ! MÓDULO 4: OTIMIZAÇÃO UNIDIMENSIONAL (RAZÃO ÁUREA)
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
    f1 = f_objetivo(x1)
    f2 = f_objetivo(x2)

    do iter = 1, max_iter
       if (abs(b - a) < tol) exit
       
       if (f1 > f2) then
          a = x1; x1 = x2; f1 = f2
          x2 = a + tau * (b - a)
          f2 = f_objetivo(x2)
       else
          b = x2; x2 = x1; f2 = f1
          x1 = b - tau * (b - a)
          f1 = f_objetivo(x1)
       end if
    end do

    x_opt = (a + b) / 2.0d0
    f_opt = f_objetivo(x_opt)
  end subroutine razao_aurea

  ! =========================================================================
  ! EQUAÇÕES E FUNÇÕES DO PROBLEMA (ALTERAR APENAS AQUI NA HORA DA PROVA!)
  ! =========================================================================

  ! 1. Função para busca de Raízes: f(x) = 0
  real*8 function f_zero(x)
    implicit none
    real*8, intent(in) :: x
    f_zero = x**2 - 4.0d0 ! Exemplo: x^2 - 4 = 0 -> Raiz x = 2
  end function f_zero

  ! 2. Derivada da Função para Newton-Raphson: f'(x)
  real*8 function f_derivada(x)
    implicit none
    real*8, intent(in) :: x
    f_derivada = 2.0d0 * x ! Derivada de x^2 - 4
  end function f_derivada

  ! 3. Função Derivada para EDO: dy/dt = f(t, y)
  real*8 function f_edo(t, y)
    implicit none
    real*8, intent(in) :: t, y
    f_edo = -2.0d0 * y + t ! Exemplo: dy/dt = -2y + t
  end function f_edo

  ! 4. Função Objetivo para Otimização: Minimizar f(x)
  real*8 function f_objetivo(x)
    implicit none
    real*8, intent(in) :: x
    f_objetivo = (x - 2.5d0)**2 + 3.0d0 ! Exemplo: Minimo em x = 2.5
  end function f_objetivo

end program master_cna_toolkit
