! ================================================================
! ATIVIDADE 3 - NEWTON-RAPHSON PARA SISTEMAS NAO LINEARES
!
! APLICACAO: REDE DE TUBULACOES COM DOIS RAMOS EM PARALELO
!
! OBJETIVO:
! Calcular as vazoes Q1 e Q2 em dois ramos de uma tubulacao,
! utilizando o metodo de Newton-Raphson para sistemas nao lineares.
!
! O metodo utiliza uma matriz Jacobiana para calcular as
! correcoes das variaveis em cada iteracao.
!
! MODELO UTILIZADO:
!
! Equacao 1 - Conservacao de massa:
! Q1 + Q2 - Qt = 0
!
! Equacao 2 - Igualdade das perdas de pressao:
! R1*Q1**2 - R2*Q2**2 = 0
!
! Dados ilustrativos:
! Qt = 0.03 m3/s
! R1 = 200000
! R2 = 800000
!
! ATENCAO:
! Os coeficientes R1 e R2 sao parametros ilustrativos do modelo.
! Em uma aplicacao real, devem ser obtidos a partir das
! caracteristicas hidraulicas de cada ramo.
!
! DIFERENCA EM RELACAO AOS EXERCICIOS ANTERIORES:
!
! Bisseccao e falsa posicao:
! Calculam uma raiz de uma equacao usando um intervalo.
!
! Newton-Raphson e secante:
! Calculam uma raiz de uma equacao utilizando estimativas.
!
! Newton-Raphson para sistemas:
! Calcula varias incognitas simultaneamente, utilizando
! um sistema de equacoes e sua matriz Jacobiana.
!
! ARQUIVO: atividade_03_newton_sistemas.md
! ================================================================

program atividade_03

    implicit none

    ! IMPLICIT NONE:
    ! Obriga a declarar todas as variaveis.
    ! Evita erros provocados por nomes digitados incorretamente.

    ! REAL(8):
    ! Utiliza precisao dupla nos calculos.
    ! E adequado para reduzir erros de arredondamento.

    ! INTEGER:
    ! Armazena numeros inteiros, como o contador de iteracoes.

    ! LOGICAL:
    ! Armazena valores .TRUE. ou .FALSE.
    ! Sera utilizado para indicar se o metodo convergiu.

    real(8) :: q1, q2
    real(8) :: q1_novo, q2_novo
    real(8) :: qt, r1, r2
    real(8) :: f1, f2
    real(8) :: j11, j12, j21, j22
    real(8) :: determinante
    real(8) :: delta1, delta2
    real(8) :: tol, erro
    integer :: iteracao, max_iter

    logical :: convergiu

    ! ============================================================
    ! 1. DEFINICAO DOS DADOS DO PROBLEMA
    ! ============================================================

    ! Vazao total que entra na tubulacao.
    qt = 0.03d0

    ! Coeficientes de resistencia hidraulica dos dois ramos.
    r1 = 200000.0d0
    r2 = 800000.0d0

    ! Estimativas iniciais das vazoes.
    !
    ! Newton-Raphson para sistemas precisa de um chute inicial
    ! para cada incognita.
    !
    ! Diferentemente da bisseccao, nao e necessario fornecer
    ! um intervalo para cada variavel.

    q1 = 0.015d0
    q2 = 0.015d0

    ! Tolerancia utilizada no criterio de convergencia.
    tol = 1.0d-10

    ! Numero maximo de iteracoes.
    max_iter = 100

    convergiu = .false.

    ! ============================================================
    ! 2. APRESENTACAO DOS DADOS
    ! ============================================================

    print *, "=============================================="
    print *, " ATIVIDADE 3 - NEWTON-RAPHSON PARA SISTEMAS"
    print *, "=============================================="

    print *, "Vazao total:", qt
    print *, "Coeficiente R1:", r1
    print *, "Coeficiente R2:", r2
    print *, "Estimativa inicial Q1:", q1
    print *, "Estimativa inicial Q2:", q2
    print *

    ! ============================================================
    ! 3. PROCESSO ITERATIVO
    ! ============================================================

    ! DO:
    ! Repete os calculos enquanto houver iteracoes disponiveis.
    !
    ! Neste exercicio, cada iteracao atualiza Q1 e Q2 ao mesmo
    ! tempo, utilizando a matriz Jacobiana.

    do iteracao = 1, max_iter

        ! --------------------------------------------------------
        ! 3.1. CALCULO DAS EQUACOES DO SISTEMA
        ! --------------------------------------------------------

        ! Equacao 1:
        ! A soma das vazoes dos ramos deve ser igual a vazao total.

        f1 = q1 + q2 - qt

        ! Equacao 2:
        ! Para ramos em paralelo, as perdas de pressao devem ser
        ! iguais, conforme o modelo quadratico adotado.

        f2 = r1*q1**2 - r2*q2**2

        ! --------------------------------------------------------
        ! 3.2. CONSTRUCAO DA MATRIZ JACOBIANA
        ! --------------------------------------------------------

        ! A matriz Jacobiana contem as derivadas parciais
        ! das equacoes em relacao as incognitas.
        !
        ! Sistema:
        !
        ! F1(Q1,Q2) = Q1 + Q2 - Qt
        ! F2(Q1,Q2) = R1*Q1**2 - R2*Q2**2
        !
        ! Jacobiana:
        !
        !        | dF1/dQ1   dF1/dQ2 |
        ! J =    |                    |
        !        | dF2/dQ1   dF2/dQ2 |
        !
        ! Portanto:
        !
        ! J = | 1          1           |
        !     | 2*R1*Q1   -2*R2*Q2     |

        j11 = 1.0d0
        j12 = 1.0d0
        j21 = 2.0d0*r1*q1
        j22 = -2.0d0*r2*q2

        ! --------------------------------------------------------
        ! 3.3. CALCULO DO DETERMINANTE
        ! --------------------------------------------------------

        ! Para uma matriz 2x2:
        !
        !       | a  b |
        ! J  =  |      |
        !       | c  d |
        !
        ! det(J) = a*d - b*c
        !
        ! Se o determinante for zero ou muito proximo de zero,
        ! nao sera possivel resolver o sistema linear desta etapa
        ! de maneira numericamente segura.

        determinante = j11*j22 - j12*j21

        if (abs(determinante) < 1.0d-14) then

            print *, "Erro: matriz Jacobiana singular."
            print *, "Tente outras estimativas iniciais."

            exit

        end if

        ! --------------------------------------------------------
        ! 3.4. CALCULO DAS CORRECOES
        ! --------------------------------------------------------

        ! O metodo de Newton-Raphson para sistemas utiliza:
        !
        ! J * delta = -F
        !
        ! Onde:
        ! J     = matriz Jacobiana
        ! delta = vetor de correcoes
        ! F     = vetor com os valores das equacoes
        !
        ! Para este sistema 2x2, podemos resolver diretamente
        ! utilizando a regra de Cramer.
        !
        ! delta1 = (-F1*j22 + j12*f2)/determinante
        !
        ! delta2 = (-j11*f2 + j21*f1)/determinante
        !
        ! As expressoes abaixo correspondem a solucao de:
        ! J*delta = -F.

        delta1 = (-f1*j22 + j12*f2)/determinante

        delta2 = (-j11*f2 + j21*f1)/determinante

        ! --------------------------------------------------------
        ! 3.5. ATUALIZACAO DAS INCOGNITAS
        ! --------------------------------------------------------

        ! As novas estimativas sao obtidas somando as correcoes
        ! calculadas as estimativas anteriores.

        q1_novo = q1 + delta1
        q2_novo = q2 + delta2

        ! --------------------------------------------------------
        ! 3.6. CALCULO DO ERRO
        ! --------------------------------------------------------

        ! MAX:
        ! Retorna o maior valor entre os argumentos.
        !
        ! ABS:
        ! Calcula o valor absoluto.
        !
        ! Utilizamos a maior correcao absoluta para verificar
        ! se ambas as incognitas atingiram a precisao desejada.

        erro = max(abs(delta1), abs(delta2))

        ! Apresenta os resultados de cada iteracao.
        print *, "Iteracao:", iteracao
        print *, "Q1:", q1_novo
        print *, "Q2:", q2_novo
        print *, "Erro:", erro
        print *

        ! --------------------------------------------------------
        ! 3.7. CRITERIO DE PARADA
        ! --------------------------------------------------------

        ! IF:
        ! Executa o bloco quando a condicao for verdadeira.
        !
        ! Se a maior correcao for menor que a tolerancia,
        ! consideramos que o metodo convergiu.

        if (erro < tol) then

            q1 = q1_novo
            q2 = q2_novo

            convergiu = .true.

            exit

        end if

        ! Atualiza as estimativas para a proxima iteracao.
        q1 = q1_novo
        q2 = q2_novo

    end do

    ! ============================================================
    ! 4. APRESENTACAO DOS RESULTADOS FINAIS
    ! ============================================================

    print *, "=============================================="
    print *, " RESULTADOS FINAIS"
    print *, "=============================================="

    if (convergiu) then

        print *, "Metodo convergiu."
        print *, "Iteracoes realizadas:", iteracao
        print *, "Vazao no ramo 1 (m3/s):", q1
        print *, "Vazao no ramo 2 (m3/s):", q2

        ! Verifica novamente as equacoes para avaliar os residuos.
        !
        ! Residuos proximos de zero indicam que as equacoes
        ! estao sendo satisfeitas pela solucao encontrada.

        f1 = q1 + q2 - qt
        f2 = r1*q1**2 - r2*q2**2

        print *, "Residuo da equacao 1:", f1
        print *, "Residuo da equacao 2:", f2

        print *
        print *, "Verificacao da conservacao de massa:"
        print *, q1 + q2

    else

        print *, "O metodo nao convergiu dentro do limite."
        print *, "Verifique os dados e as estimativas iniciais."

    end if

    print *
    print *, "Fim do programa."

end program atividade_03
