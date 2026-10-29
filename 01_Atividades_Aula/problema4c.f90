! ======================================================================
! TOOLKIT MESTRE CNA - FORTRAN 95
! Metodos: bissecao, secante, Newton, Gauss com pivoteamento, LU com
! pivoteamento (PA=LU), Thomas, RK4 para EDO escalar y'=f(t,y) e razao aurea.
!
! ATENCAO: os exemplos matematicos sao demonstrativos. Na prova, adapte
! as funcoes no BLOCO 6 e confirme o modelo fisico da aleta.
! A equacao classica da aleta e de segunda ordem; a rotina RK4 daqui
! resolve uma EDO escalar de primeira ordem e NAO resolve a aleta diretamente.
!
! Compilar:
! gfortran -std=f95 -Wall -Wextra -pedantic cna_toolkit_mestre.f90 -o cna_toolkit_mestre
! ======================================================================
program cna_toolkit_mestre
  implicit none
  integer, parameter :: dp = selected_real_kind(15,300)
  integer, parameter :: ns=3, nt=4
  real(dp) :: rb, rs, rn, amat(ns,ns), bvec(ns), xg(ns), xl(ns)
  real(dp) :: lmat(ns,ns), umat(ns,ns), elo(nt), diag(nt), ghi(nt)
  real(dp) :: rhs(nt), xt(nt), yfinal, xopt, fopt
  integer :: status

  print *, '=== TOOLKIT MESTRE CNA - TESTES DEMONSTRATIVOS ==='
  print *, 'As funcoes de teste devem ser trocadas pelo modelo da prova.'
  print *

  ! TESTE 1: f(x)=x**2-4; raiz positiva conhecida x=2.
  call bissecao(1.0_dp,3.0_dp,1.0e-8_dp,200,rb,status)
  call mostra_status('Bissecao',status)
  print '(A,F12.8)', 'Raiz = ',rb
  call secante(1.0_dp,3.0_dp,1.0e-8_dp,200,rs,status)
  call mostra_status('Secante',status)
  print '(A,F12.8)', 'Raiz = ',rs
  call newton(1.5_dp,1.0e-8_dp,200,rn,status)
  call mostra_status('Newton-Raphson',status)
  print '(A,F12.8)', 'Raiz = ',rn
  print *

  ! TESTE 2: sistema A*x=b. Gauss e LU preservam A e b de entrada.
  amat(1,:)=(/3.0_dp,-0.1_dp,-0.2_dp/)
  amat(2,:)=(/0.1_dp,7.0_dp,-0.3_dp/)
  amat(3,:)=(/0.3_dp,-0.2_dp,10.0_dp/)
  bvec=(/7.85_dp,-19.30_dp,71.40_dp/)
  call gauss_pp(amat,bvec,ns,xg,status)
  call mostra_status('Gauss com pivoteamento',status)
  print '(A,3F12.6)', 'x = ',xg
  print '(A,ES12.4)', 'Residuo infinito = ',residuo_inf(amat,bvec,xg)
  call lu_pp(amat,bvec,ns,xl,lmat,umat,status)
  call mostra_status('LU com pivoteamento',status)
  print '(A,3F12.6)', 'x = ',xl
  print '(A,ES12.4)', 'Residuo infinito = ',residuo_inf(amat,bvec,xl)
  print *

  ! TESTE 3: Thomas exige matriz tridiagonal e pivos adequados.
  diag=(/2.04_dp,2.04_dp,2.04_dp,2.04_dp/)
  elo=(/0.0_dp,-1.0_dp,-1.0_dp,-1.0_dp/)
  ghi=(/-1.0_dp,-1.0_dp,-1.0_dp,0.0_dp/)
  rhs=(/40.8_dp,0.8_dp,0.8_dp,200.8_dp/)
  call thomas(elo,diag,ghi,rhs,nt,xt,status)
  call mostra_status('Thomas',status)
  print '(A,4F12.6)','x = ',xt
  print *

  ! TESTE 4: RK4 para y'=-2y+t, y(0)=1; exemplo de primeira ordem.
  call rk4(0.0_dp,2.0_dp,1.0_dp,0.01_dp,yfinal,status)
  call mostra_status('RK4',status)
  print '(A,F12.8)','y(tf) = ',yfinal
  print *

  ! TESTE 5: minimizacao de (x-2.5)^2+3 em [0,4].
  call razao_aurea(0.0_dp,4.0_dp,1.0e-8_dp,200,xopt,fopt,status)
  call mostra_status('Razao aurea',status)
  print '(A,F12.8)','x minimo = ',xopt
  print '(A,F12.8)','f minimo = ',fopt
  print *
  print *, '=== FIM DOS TESTES ==='

contains

  ! STATUS: 0=sucesso, 2=intervalo invalido, 3=max_iter, 4=derivada/
  ! denominador inadequado, 5=matriz singular/pivo pequeno, 1=entrada invalida.

  ! ---------------- BLOCO 1: RAIZES DE FUNCOES ----------------
  subroutine bissecao(ain,bin,tol,maxit,x,status)
    implicit none
    real(dp),intent(in)::ain,bin,tol
    integer,intent(in)::maxit
    real(dp),intent(out)::x
    integer,intent(out)::status
    real(dp)::a,b,c,fa,fb,fc
    integer::it
    x=ain; status=1
    if(tol<=0.0_dp.or.maxit<1.or.ain>=bin)return
    a=ain; b=bin; fa=f_raiz(a); fb=f_raiz(b)
    if(abs(fa)<=tol)then
      x=a; status=0; return
    endif
    if(abs(fb)<=tol)then
      x=b; status=0; return
    endif
    if(fa*fb>0.0_dp)then
      status=2; return
    endif
    do it=1,maxit
      c=a+0.5_dp*(b-a)
      fc=f_raiz(c); x=c
      if(abs(fc)<=tol.or.0.5_dp*abs(b-a)<=tol)then
        status=0; return
      endif
      if(fa*fc<0.0_dp)then
        b=c; fb=fc
      else
        a=c; fa=fc
      endif
    enddo
    status=3
  end subroutine bissecao

  subroutine secante(x0,x1,tol,maxit,x,status)
    implicit none
    real(dp),intent(in)::x0,x1,tol
    integer,intent(in)::maxit
    real(dp),intent(out)::x
    integer,intent(out)::status
    real(dp)::a,b,fa,fb,den,xn
    integer::it
    a=x0; b=x1; x=x1; status=1
    if(tol<=0.0_dp.or.maxit<1)return
    fa=f_raiz(a); fb=f_raiz(b)
    do it=1,maxit
      if(abs(fb)<=tol)then
        x=b; status=0; return
      endif
      den=fb-fa
      if(abs(den)<=tiny(den))then
        x=b; status=4; return
      endif
      xn=b-fb*(b-a)/den
      if(abs(xn-b)<=tol.and.abs(f_raiz(xn))<=sqrt(tol))then
        x=xn; status=0; return
      endif
      a=b; fa=fb; b=xn; fb=f_raiz(b)
    enddo
    x=b; status=3
  end subroutine secante

  subroutine newton(x0,tol,maxit,x,status)
    implicit none
    real(dp),intent(in)::x0,tol
    integer,intent(in)::maxit
    real(dp),intent(out)::x
    integer,intent(out)::status
    real(dp)::fx,dfx,dx
    integer::it
    x=x0; status=1
    if(tol<=0.0_dp.or.maxit<1)return
    do it=1,maxit
      fx=f_raiz(x)
      if(abs(fx)<=tol)then
        status=0; return
      endif
      dfx=df_raiz(x)
      if(abs(dfx)<=sqrt(tiny(dfx)))then
        status=4; return
      endif
      dx=fx/dfx
      x=x-dx
      if(abs(dx)<=tol.and.abs(f_raiz(x))<=sqrt(tol))then
        status=0; return
      endif
    enddo
    status=3
  end subroutine newton

  ! ---------------- BLOCO 2: GAUSS COM PIVOTEAMENTO ----------------
  subroutine gauss_pp(ain,bin,n,x,status)
    implicit none
    integer,intent(in)::n
    real(dp),intent(in)::ain(n,n),bin(n)
    real(dp),intent(out)::x(n)
    integer,intent(out)::status
    real(dp)::a(n,n),b(n),factor,temp,s,scale
    integer::i,j,k,piv
    x=0.0_dp; status=1
    if(n<1)return
    a=ain; b=bin; scale=max(1.0_dp,maxval(abs(a)))
    do k=1,n-1
      piv=k
      do i=k+1,n
        if(abs(a(i,k))>abs(a(piv,k)))piv=i
      enddo
      if(abs(a(piv,k))<=epsilon(1.0_dp)*scale)then
        status=5; return
      endif
      if(piv/=k)then
        do j=1,n
          temp=a(k,j); a(k,j)=a(piv,j); a(piv,j)=temp
        enddo
        temp=b(k); b(k)=b(piv); b(piv)=temp
      endif
      do i=k+1,n
        factor=a(i,k)/a(k,k); a(i,k)=0.0_dp
        do j=k+1,n
          a(i,j)=a(i,j)-factor*a(k,j)
        enddo
        b(i)=b(i)-factor*b(k)
      enddo
    enddo
    if(abs(a(n,n))<=epsilon(1.0_dp)*scale)then
      status=5; return
    endif
    x(n)=b(n)/a(n,n)
    do i=n-1,1,-1
      s=b(i)
      do j=i+1,n
        s=s-a(i,j)*x(j)
      enddo
      if(abs(a(i,i))<=epsilon(1.0_dp)*scale)then
        status=5; return
      endif
      x(i)=s/a(i,i)
    enddo
    status=0
  end subroutine gauss_pp

  ! ---------------- BLOCO 3: LU COM PIVOTEAMENTO (P*A=L*U) --------
  subroutine lu_pp(ain,bin,n,x,l,u,status)
    implicit none
    integer,intent(in)::n
    real(dp),intent(in)::ain(n,n),bin(n)
    real(dp),intent(out)::x(n),l(n,n),u(n,n)
    integer,intent(out)::status
    real(dp)::a(n,n),b(n),temp,factor,s,y(n),scale
    integer::i,j,k,piv
    x=0.0_dp; l=0.0_dp; u=0.0_dp; status=1
    if(n<1)return
    a=ain; b=bin; scale=max(1.0_dp,maxval(abs(a)))
    do i=1,n
      l(i,i)=1.0_dp
    enddo
    do k=1,n-1
      piv=k
      do i=k+1,n
        if(abs(a(i,k))>abs(a(piv,k)))piv=i
      enddo
      if(abs(a(piv,k))<=epsilon(1.0_dp)*scale)then
        status=5; return
      endif
      if(piv/=k)then
        do j=1,n
          temp=a(k,j); a(k,j)=a(piv,j); a(piv,j)=temp
        enddo
        temp=b(k); b(k)=b(piv); b(piv)=temp
        do j=1,k-1
          temp=l(k,j); l(k,j)=l(piv,j); l(piv,j)=temp
        enddo
      endif
      do i=k+1,n
        factor=a(i,k)/a(k,k); l(i,k)=factor; a(i,k)=0.0_dp
        do j=k+1,n
          a(i,j)=a(i,j)-factor*a(k,j)
        enddo
      enddo
    enddo
    if(abs(a(n,n))<=epsilon(1.0_dp)*scale)then
      status=5; return
    endif
    u=a
    ! b ja foi permutado; resolve L*y=P*b.
    do i=1,n
      s=b(i)
      do j=1,i-1
        s=s-l(i,j)*y(j)
      enddo
      y(i)=s/l(i,i)
    enddo
    x(n)=y(n)/u(n,n)
    do i=n-1,1,-1
      s=y(i)
      do j=i+1,n
        s=s-u(i,j)*x(j)
      enddo
      if(abs(u(i,i))<=epsilon(1.0_dp)*scale)then
        status=5; return
      endif
      x(i)=s/u(i,i)
    enddo
    status=0
  end subroutine lu_pp

  ! ---------------- BLOCO 4: THOMAS PARA MATRIZ TRIDIAGONAL -------
  subroutine thomas(ein,fin,gin,rin,n,x,status)
    implicit none
    integer,intent(in)::n
    real(dp),intent(in)::ein(n),fin(n),gin(n),rin(n)
    real(dp),intent(out)::x(n)
    integer,intent(out)::status
    real(dp)::e(n),f(n),g(n),r(n),fac,scale
    integer::k
    x=0.0_dp; status=1
    if(n<1)return
    e=ein; f=fin; g=gin; r=rin
    scale=max(1.0_dp,maxval(abs(f)))
    do k=2,n
      if(abs(f(k-1))<=epsilon(1.0_dp)*scale)then
        status=5; return
      endif
      fac=e(k)/f(k-1)
      f(k)=f(k)-fac*g(k-1)
      r(k)=r(k)-fac*r(k-1)
    enddo
    if(abs(f(n))<=epsilon(1.0_dp)*scale)then
      status=5; return
    endif
    x(n)=r(n)/f(n)
    do k=n-1,1,-1
      if(abs(f(k))<=epsilon(1.0_dp)*scale)then
        status=5; return
      endif
      x(k)=(r(k)-g(k)*x(k+1))/f(k)
    enddo
    status=0
  end subroutine thomas

  ! ---------------- BLOCO 5: RK4 PARA y'=f(t,y) -------------------
  subroutine rk4(t0,tf,y0,h0,yf,status)
    implicit none
    real(dp),intent(in)::t0,tf,y0,h0
    real(dp),intent(out)::yf
    integer,intent(out)::status
    real(dp)::t,y,h,k1,k2,k3,k4
    integer::steps
    yf=y0; status=1
    if(h0<=0.0_dp.or.tf<t0)return
    t=t0; y=y0; steps=0
    do while(t<tf)
      h=min(h0,tf-t)
      k1=func_edo(t,y)
      k2=func_edo(t+0.5_dp*h,y+0.5_dp*h*k1)
      k3=func_edo(t+0.5_dp*h,y+0.5_dp*h*k2)
      k4=func_edo(t+h,y+h*k3)
      y=y+(h/6.0_dp)*(k1+2.0_dp*k2+2.0_dp*k3+k4)
      t=t+h; steps=steps+1
      if(steps>10000000)then
        yf=y; status=6; return
      endif
    enddo
    yf=y; status=0
  end subroutine rk4

  ! ---------------- BLOCO 6: RAZAO AUREA PARA MINIMIZACAO ---------
  subroutine razao_aurea(ain,bin,tol,maxit,xmin,fmin,status)
    implicit none
    real(dp),intent(in)::ain,bin,tol
    integer,intent(in)::maxit
    real(dp),intent(out)::xmin,fmin
    integer,intent(out)::status
    real(dp)::a,b,tau,x1,x2,y1,y2
    integer::it
    xmin=ain; fmin=func_objetivo(ain); status=1
    if(tol<=0.0_dp.or.maxit<1.or.ain>=bin)return
    a=ain; b=bin; tau=(sqrt(5.0_dp)-1.0_dp)/2.0_dp
    x1=b-tau*(b-a); x2=a+tau*(b-a)
    y1=func_objetivo(x1); y2=func_objetivo(x2)
    do it=1,maxit
      if(abs(b-a)<=tol)then
        xmin=0.5_dp*(a+b); fmin=func_objetivo(xmin); status=0; return
      endif
      if(y1<=y2)then
        b=x2; x2=x1; y2=y1
        x1=b-tau*(b-a); y1=func_objetivo(x1)
      else
        a=x1; x1=x2; y1=y2
        x2=a+tau*(b-a); y2=func_objetivo(x2)
      endif
    enddo
    xmin=0.5_dp*(a+b); fmin=func_objetivo(xmin); status=3
  end subroutine razao_aurea

  ! ==================================================================
  ! BLOCO 7: EDITE ESTAS FUNCOES PARA O ENUNCIADO DA PROVA.
  ! 1) f_raiz: equacao f(x)=0.
  ! 2) df_raiz: derivada de f_raiz, usada por Newton.
  ! 3) func_edo: EDO escalar y'=f(t,y), se aplicavel.
  ! 4) func_objetivo: funcao a minimizar no intervalo.
  ! Sempre confira unidades, condicoes de contorno e significado fisico.
  ! ==================================================================
  real(dp) function f_raiz(x)
    implicit none
    real(dp),intent(in)::x
    f_raiz=x**2-4.0_dp
  end function f_raiz

  real(dp) function df_raiz(x)
    implicit none
    real(dp),intent(in)::x
    df_raiz=2.0_dp*x
  end function df_raiz

  real(dp) function func_edo(t,y)
    implicit none
    real(dp),intent(in)::t,y
    func_edo=-2.0_dp*y+t
  end function func_edo

  real(dp) function func_objetivo(x)
    implicit none
    real(dp),intent(in)::x
    func_objetivo=(x-2.5_dp)**2+3.0_dp
  end function func_objetivo

  ! ---------------- ROTINAS DE DIAGNOSTICO ------------------------
  subroutine mostra_status(nome,st)
    implicit none
    character(len=*),intent(in)::nome
    integer,intent(in)::st
    select case(st)
    case(0)
      print '(A,A)',trim(nome),': OK/convergiu'
    case(2)
      print '(A,A)',trim(nome),': ERRO - intervalo sem troca de sinal'
    case(3)
      print '(A,A)',trim(nome),': AVISO - maximo de iteracoes atingido'
    case(4)
      print '(A,A)',trim(nome),': ERRO - derivada ou denominador inadequado'
    case(5)
      print '(A,A)',trim(nome),': ERRO - matriz singular ou pivo pequeno'
    case(6)
      print '(A,A)',trim(nome),': ERRO - limite de passos excedido'
    case default
      print '(A,A)',trim(nome),': ERRO - argumentos invalidos'
    end select
  end subroutine mostra_status

  real(dp) function residuo_inf(a,b,x)
    implicit none
    real(dp),intent(in)::a(:,:),b(:),x(:)
    real(dp)::r(size(b))
    r=matmul(a,x)-b
    residuo_inf=maxval(abs(r))
  end function residuo_inf

end program cna_toolkit_mestre
