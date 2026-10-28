! ================================================================
! EXERCICIO 2 - METODOS DE NEWTON-RAPHSON E DA SECANTE
!
! OBJETIVO:
! Encontrar a raiz de uma equacao nao linear utilizando dois
! metodos numericos e comparar os resultados.
!
! EQUACAO:
! f(x) = x**3 - x - 2
!
! DIFERENCA EM RELACAO AOS METODOS DA BISECCAO E FALSA POSICAO:
! Bisseccao e falsa posicao utilizam um intervalo [a,b].
! Newton-Raphson utiliza uma estimativa inicial e a derivada.
! Secante utiliza duas estimativas iniciais e dispensa a derivada.
!
! Para adaptar o programa a outra equacao:
! 1. Altere a funcao F.
! 2. Altere a derivada DF para Newton-Raphson.
! 3. Escolha estimativas iniciais adequadas.
! ================================================================

program exercicio_2

    implicit none

    ! IMPLICIT NONE:
    ! Obriga a declarar todas as variaveis.
    ! Ajuda a evitar erros de digitacao nos nomes.

    ! REAL(8):
    ! Utiliza precisao dupla para os calculos numericos.
    !
    ! INTEGER:
    ! Armazena numeros inteiros, como o numero de iteracoes.
    !
    ! LOGICAL:
    ! Armazena .TRUE. ou .FALSE. para indicar convergencia.

    real(8) :: x0, x1
    real(8) :: raiz_newton, raiz_secante
    real(8) :: tolerancia
    integer :: max_iter

    logical :: convergiu_newton, convergiu_secante

    ! Tolerancia define a precisao desejada.
    ! Max_iter limita o numero maximo de repeticoes.

    tolerancia = 1.0d-8
    max_iter = 100

    print *, "=============================================="
    print *, " EXERCICIO 2 - CALCULO NUMERICO EM FORTRAN"
    print *, "=============================================="
    print *, "Equacao: f(x) = x^3 - x - 2"
    print *, "Tolerancia:", tolerancia
    print *, "Maximo de iteracoes:", max_iter
    print *

    ! ============================================================
    ! PARTE 1 - NEWTON-RAPHSON
    ! ============================================================

    print *, "METODO DE NEWTON-RAPHSON"

    ! Formula:
    ! x(n+1) = x(n) - f(x(n))/f'(x(n))
    !
    ! Para esta funcao:
    ! f'(x) = 3*x**2 - 1
    !
    ! O metodo precisa de uma estimativa inicial.

    x0 = 1.5d0

    call newton_raphson(x0, tolerancia, max_iter, &
                        raiz_newton, convergiu_newton)

    if (convergiu_newton) then
        print *, "Raiz encontrada:", raiz_newton
        print *, "Valor de f(raiz):", f(raiz_newton)
    else
        print *, "Newton-Raphson nao convergiu."
    end if

    print *

    ! ============================================================
    ! PARTE 2 - METODO DA SECANTE
    ! ============================================================

    print *, "METODO DA SECANTE"

    ! Formula:
    ! x(n+1) = x(n) - f(x(n))*(x(n)-x(n-1)) /
    !          (f(x(n))-f(x(n-1)))
    !
    ! Diferentemente de Newton-Raphson, nao precisamos calcular
    ! a derivada analitica.
    !
    ! Sao necessarias duas estimativas iniciais.

    x0 = 1.0d0
    x1 = 2.0d0

    call metodo_secante(x0, x1, tolerancia, max_iter, &
                        raiz_secante, convergiu_secante)

    if (convergiu_secante) then
        print *, "Raiz encontrada:", raiz_secante
        print *, "Valor de f(raiz):", f(raiz_secante)
    else
        print *, "Metodo da secante nao convergiu."
    end if

    print *

    ! ============================================================
    ! PARTE 3 - COMPARACAO DOS RESULTADOS
    ! ============================================================

    print *, "COMPARACAO DOS METODOS"

    ! IF / ELSE:
    ! Executa comandos diferentes conforme uma condicao.
    !
    ! .AND. significa E logico.
    ! ABS calcula o valor absoluto.

    if (convergiu_newton .and. convergiu_secante) then

        print *, "Os dois metodos convergiram."
        print *, "Newton-Raphson:", raiz_newton
        print *, "Secante:", raiz_secante
        print *, "Diferenca entre as raizes:"
        print *, abs(raiz_newton - raiz_secante)

    else

        print *, "Pelo menos um metodo nao convergiu."

    end if

    print *
    print *, "Fim do programa."

contains

    ! ============================================================
    ! FUNCAO: f(x) = x**3 - x - 2
    ! ============================================================

    real(8) function f(x)

        implicit none
        real(8), intent(in) :: x

        ! FUNCTION calcula e retorna um valor.
        ! INTENT(IN) indica que x e somente uma entrada.
        ! ** representa potenciacao.

        f = x**3 - x - 2.0d0

    end function f

    ! ============================================================
    ! DERIVADA: f'(x) = 3*x**2 - 1
    ! ============================================================

    real(8) function df(x)

        implicit none
        real(8), intent(in) :: x

        ! Newton-Raphson utiliza a derivada.
        ! Se F mudar, DF tambem deve ser alterada.
        ! A secante nao precisa chamar esta funcao.

        df = 3.0d0*x**2 - 1.0d0

    end function df

    ! ============================================================
    ! SUBROTINA - NEWTON-RAPHSON
    ! ============================================================

    subroutine newton_raphson(x_inicial, tol, nmax, raiz, sucesso)

        implicit none

        real(8), intent(in) :: x_inicial, tol
        integer, intent(in) :: nmax

        real(8), intent(out) :: raiz
        logical, intent(out) :: sucesso

        real(8) :: x_atual, x_novo
        real(8) :: derivada
        integer :: i

        ! SUBROUTINE organiza uma tarefa separada.
        ! CALL executa a subrotina.
        ! DO repete o calculo ate atingir o limite de iteracoes.

        x_atual = x_inicial
        sucesso = .false.

        do i = 1, nmax

            derivada = df(x_atual)

            ! Evita divisao por uma derivada proxima de zero.
            if (abs(derivada) < 1.0d-14) then
                raiz = x_atual
                return
            end if

            ! Atualiza a estimativa usando Newton-Raphson.
            x_novo = x_atual - f(x_atual)/derivada

            ! Criterio de parada:
            ! A diferenca entre duas aproximacoes deve ser menor
            ! que a tolerancia estabelecida.
            !
            ! EXIT encerra o laco imediatamente.

            if (abs(x_novo - x_atual) < tol) then
                raiz = x_novo
                sucesso = .true.
                exit
            end if

            ! Atualiza a estimativa para a proxima iteracao.
            x_atual = x_novo

        end do

        if (.not. sucesso) then
            raiz = x_atual
        end if

    end subroutine newton_raphson

    ! ============================================================
    ! SUBROTINA - METODO DA SECANTE
    ! ============================================================

    subroutine metodo_secante(x_inicial0, x_inicial1, tol, &
                              nmax, raiz, sucesso)

        implicit none

        real(8), intent(in) :: x_inicial0, x_inicial1, tol
        integer, intent(in) :: nmax

        real(8), intent(out) :: raiz
        logical, intent(out) :: sucesso

        real(8) :: x_anterior, x_atual, x_novo
        real(8) :: f_anterior, f_atual, denominador
        integer :: i

        ! A secante trabalha com duas aproximacoes.
        ! A cada iteracao, os pontos sao atualizados.

        x_anterior = x_inicial0
        x_atual = x_inicial1

        sucesso = .false.

        do i = 1, nmax

            f_anterior = f(x_anterior)
            f_atual = f(x_atual)

            denominador = f_atual - f_anterior

            ! Evita divisao por zero ou por valor muito pequeno.
            if (abs(denominador) < 1.0d-14) then
                raiz = x_atual
                return
            end if

            ! Formula do metodo da secante.
            x_novo = x_atual - f_atual * &
                     (x_atual - x_anterior) / denominador

            ! Verifica se a aproximacao atingiu a tolerancia.
            if (abs(x_novo - x_atual) < tol) then
                raiz = x_novo
                sucesso = .true.
                exit
            end if

            ! Move os pontos para a proxima iteracao.
            x_anterior = x_atual
            x_atual = x_novo

        end do

        if (.not. sucesso) then
            raiz = x_atual
        end if

    end subroutine metodo_secante

end program exercicio_2