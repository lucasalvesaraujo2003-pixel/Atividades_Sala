! ================================================================
! ATIVIDADE 2 - METODOS DE NEWTON-RAPHSON E DA SECANTE
! Disciplina: Calculo Numerico Aplicado
!
! Objetivo:
! Encontrar a temperatura de equilibrio de um sistema termico
! utilizando os metodos de Newton-Raphson e da Secante.
!
! Equacao:
! f(T) = h*A*(T - Tinf) + epsilon*sigma*A*(T**4 - Tinf**4) - P
!
! Newton-Raphson:
! Tnovo = T - f(T)/f'(T)
!
! Secante:
! Tnovo = T2 - f(T2)*(T2-T1)/(f(T2)-f(T1))
!
! Criterio de parada:
! |f(T)| < TOL
!
! Comparacao com a Atividade 1:
! Na bissecao e na falsa posicao, utilizamos um intervalo [a,b].
! No metodo de Newton-Raphson, precisamos da derivada da funcao.
! Na secante, utilizamos duas aproximacoes iniciais, sem derivada.
! ================================================================

program atividade_02_newton_secante

    ! Impede que variaveis sejam criadas implicitamente.
    ! Assim, erros de digitacao nos nomes podem ser identificados.
    ! Use implicit none em praticamente todos os programas Fortran.
    implicit none

    ! Define a precisao numerica utilizada nas variaveis reais.
    ! kind(1.0d0) seleciona a precisao de um numero real de dupla
    ! precisao, adequada para reduzir erros de arredondamento.
    ! Na Atividade 1, os valores reais podiam estar declarados
    ! diretamente como REAL; aqui deixamos a precisao explicita.
    integer, parameter :: dp = kind(1.0d0)

    ! ------------------------------------------------------------
    ! Declaracao das variaveis
    ! ------------------------------------------------------------

    ! Dados fisicos do problema.
    ! PARAMETER indica que os valores sao constantes durante a
    ! execucao. Diferentemente de uma variavel comum, nao mudam.
    real(dp), parameter :: A = 0.01_dp
    real(dp), parameter :: h = 10.0_dp
    real(dp), parameter :: epsilon = 0.8_dp
    real(dp), parameter :: Tinf = 300.0_dp
    real(dp), parameter :: P = 20.0_dp
    real(dp), parameter :: sigma = 5.67e-8_dp

    ! TOL e MAX_ITER controlam a precisao e o custo computacional.
    ! A tolerancia e semelhante a da Atividade 1, mas o criterio
    ! aqui e aplicado ao residuo da equacao, |f(T)|.
    real(dp), parameter :: TOL = 1.0e-6_dp
    integer, parameter :: MAX_ITER = 1000

    ! Variaveis do metodo de Newton-Raphson.
    real(dp) :: T_newton
    real(dp) :: T_prox
    real(dp) :: f_atual
    real(dp) :: df_atual
    integer :: iter_newton
    logical :: convergiu_newton

    ! Variaveis do metodo da Secante.
    real(dp) :: T1
    real(dp) :: T2
    real(dp) :: T3
    real(dp) :: f1
    real(dp) :: f2
    real(dp) :: denominador
    integer :: iter_secante
    logical :: convergiu_secante

    ! ------------------------------------------------------------
    ! Metodo de Newton-Raphson
    ! ------------------------------------------------------------

    print *, "=============================================="
    print *, " METODO DE NEWTON-RAPHSON"
    print *, "=============================================="

    ! Estimativa inicial da temperatura.
    ! Newton depende da escolha de uma aproximacao inicial.
    ! Diferentemente da bissecao, nao exige intervalo com mudanca
    ! de sinal, mas uma escolha ruim pode impedir a convergencia.
    T_newton = 350.0_dp

    iter_newton = 0
    convergiu_newton = .false.

    ! DO WHILE repete enquanto a condicao for verdadeira.
    ! Na Atividade 1, a mesma estrutura foi usada para repetir
    ! iteracoes ate atingir a tolerancia ou o limite de iteracoes.
    do while (iter_newton < MAX_ITER .and. .not. convergiu_newton)

        ! Calcula o valor da funcao na estimativa atual.
        f_atual = f(T_newton)

        ! Se o residuo for pequeno o suficiente, encontramos
        ! uma raiz aproximada e podemos encerrar o processo.
        if (abs(f_atual) < TOL) then
            convergiu_newton = .true.
            exit
        end if

        ! Calcula a derivada da funcao no ponto atual.
        ! Newton utiliza a derivada para estimar a raiz.
        df_atual = df(T_newton)

        ! Evita divisao por um valor praticamente nulo.
        ! Se a derivada for muito pequena, o passo de Newton pode
        ! ficar enorme e causar instabilidade numerica.
        if (abs(df_atual) < 1.0e-14_dp) then
            print *, "ERRO: derivada proxima de zero em Newton."
            exit
        end if

        ! Formula de Newton-Raphson:
        ! T_prox = T_atual - f(T_atual)/f'(T_atual)
        !
        ! Diferentemente da bissecao, nao calculamos o ponto medio.
        ! Diferentemente da secante, nao estimamos a derivada.
        T_prox = T_newton - f_atual / df_atual

        ! Atualiza a estimativa para a proxima iteracao.
        T_newton = T_prox

        iter_newton = iter_newton + 1

    end do

    ! Mostra os resultados de Newton-Raphson.
    print *, "Iteracoes:", iter_newton
    print *, "Temperatura encontrada (K):", T_newton
    print *, "Residuo final:", abs(f(T_newton))

    if (convergiu_newton) then
        print *, "Status: convergiu."
    else
        print *, "Status: nao foi confirmada a convergencia."
    end if

    ! ------------------------------------------------------------
    ! Metodo da Secante
    ! ------------------------------------------------------------

    print *
    print *, "=============================================="
    print *, " METODO DA SECANTE"
    print *, "=============================================="

    ! Duas estimativas iniciais.
    ! A secante nao precisa da derivada analitica, mas precisa
    ! de dois valores iniciais diferentes.
    ! Ao contrario da bissecao, eles nao precisam formar um
    ! intervalo garantido com mudanca de sinal.
    T1 = 300.0_dp
    T2 = 400.0_dp

    iter_secante = 0
    convergiu_secante = .false.

    do while (iter_secante < MAX_ITER .and. .not. convergiu_secante)

        f1 = f(T1)
        f2 = f(T2)

        ! Verifica se o residuo da estimativa atual ja e pequeno.
        if (abs(f2) < TOL) then
            convergiu_secante = .true.
            exit
        end if

        ! O denominador e a diferenca entre os valores da funcao.
        denominador = f2 - f1

        ! Se f(T2) e f(T1) forem praticamente iguais, a formula
        ! da secante pode dividir por zero ou gerar um passo enorme.
        if (abs(denominador) < 1.0e-14_dp) then
            print *, "ERRO: denominador proximo de zero na secante."
            exit
        end if

        ! Formula da secante:
        ! T3 = T2 - f(T2)*(T2-T1)/(f(T2)-f(T1))
        !
        ! A inclinacao e estimada usando dois pontos.
        ! Isso evita programar a derivada, mas pode ser menos
        ! robusto dependendo das estimativas iniciais.
        T3 = T2 - f2 * (T2 - T1) / denominador

        ! Desloca os pontos:
        ! o ponto antigo T2 passa a ser T1, e T3 passa a ser T2.
        ! Esse deslocamento permite repetir a formula sem guardar
        ! todos os valores anteriores.
        T1 = T2
        T2 = T3

        iter_secante = iter_secante + 1

    end do

    ! Mostra os resultados da secante.
    print *, "Iteracoes:", iter_secante
    print *, "Temperatura encontrada (K):", T2
    print *, "Residuo final:", abs(f(T2))

    if (convergiu_secante) then
        print *, "Status: convergiu."
    else
        print *, "Status: nao foi confirmada a convergencia."
    end if

    print *
    print *, "=============================================="
    print *, " FIM DO PROGRAMA"
    print *, "=============================================="

contains

    ! ------------------------------------------------------------
    ! Funcao que calcula o residuo do balanco de energia.
    !
    ! Em Fortran, CONTAINS permite declarar funcoes internas.
    ! Isso evita repetir a equacao em varios trechos do programa.
    ! A funcao recebe uma temperatura e retorna f(T).
    ! ------------------------------------------------------------
    real(dp) function f(T)

        implicit none

        ! Temperatura recebida pela funcao.
        real(dp), intent(in) :: T

        ! Balanco de energia:
        ! conveccao + radiacao - potencia fornecida.
        f = h*A*(T - Tinf) + epsilon*sigma*A*(T**4 - Tinf**4) - P

    end function f

    ! ------------------------------------------------------------
    ! Funcao que calcula a derivada analitica de f(T).
    !
    ! Derivada:
    ! f'(T) = h*A + 4*epsilon*sigma*A*T**3
    !
    ! Essa funcao e necessaria para Newton-Raphson.
    ! Para a secante, ela nao e utilizada.
    ! ------------------------------------------------------------
    real(dp) function df(T)

        implicit none

        real(dp), intent(in) :: T

        df = h*A + 4.0_dp*epsilon*sigma*A*T**3

    end function df

end program atividade_02_newton_secante