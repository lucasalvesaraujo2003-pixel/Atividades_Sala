! ================================================================
! EXERCICIO 5 - ELIMINACAO GAUSSIANA: SISTEMA DE BARRAS
! ================================================================
!
! OBJETIVO:
! Resolver um sistema de equacoes lineares usando o metodo da
! eliminacao gaussiana com pivotamento parcial.
!
! APLICACAO:
! Em problemas de engenharia estrutural, sistemas lineares podem
! surgir ao calcular deslocamentos, esforcos ou outras incognitas
! de uma estrutura formada por barras.
!
! MODELO MATEMATICO:
!
!                   A * x = b
!
! A = matriz dos coeficientes
! x = vetor das incognitas
! b = vetor dos termos independentes
!
! DIFERENCA EM RELACAO AO EXERCICIO 4:
! No exercicio 4, resolvemos um sistema associado a uma rede de
! reatores sem trocar linhas.
!
! Neste exemplo, implementamos o PIVOTAMENTO PARCIAL:
! antes de eliminar uma coluna, procuramos o maior valor absoluto
! entre os elementos disponiveis daquela coluna e trocamos as
! linhas, se necessario.
!
! Isso ajuda a evitar divisoes por valores muito pequenos e pode
! melhorar a estabilidade numerica.
!
! ATENCAO:
! Os coeficientes abaixo sao ilustrativos, nao foram confirmados
! como os valores do enunciado original da sua lista.
!
! ================================================================

program eliminacao_gaussiana_barras

    implicit none

    ! ------------------------------------------------------------
    ! 1. DECLARACAO DAS VARIAVEIS
    ! ------------------------------------------------------------

    ! Numero de equacoes e de incognitas.
    integer, parameter :: n = 5

    ! Matrizes e vetores do sistema.
    real(8) :: A(n,n), A_original(n,n)
    real(8) :: b(n), b_original(n)
    real(8) :: x(n)

    ! Variaveis auxiliares para os lacos.
    integer :: i, j, k

    ! Variaveis utilizadas no pivotamento parcial.
    integer :: linha_pivo, linha_temporaria
    real(8) :: maior_valor, valor_temporario

    ! Fator utilizado para eliminar os elementos abaixo do pivo.
    real(8) :: fator

    ! Variavel para a substituicao regressiva.
    real(8) :: soma

    ! Residuo de cada equacao e maior residuo encontrado.
    real(8) :: residuo(n), maior_residuo

    ! Tolerancia para identificar um pivo numericamente pequeno.
    real(8), parameter :: tolerancia = 1.0d-12

    ! Contador aproximado de operacoes aritmeticas.
    integer :: operacoes

    ! ------------------------------------------------------------
    ! 2. MONTAGEM DO SISTEMA LINEAR
    ! ------------------------------------------------------------

    ! Sistema ilustrativo:
    !
    !  4*x1 - x2              =  2
    ! -x1 + 4*x2 - x3         =  4
    !      - x2 + 4*x3 - x4  =  6
    !             - x3 + 4*x4 - x5 = 8
    !                    - x4 + 3*x5 = 11
    !
    ! Solucao esperada:
    ! x = [1, 2, 3, 4, 5]
    !
    ! Em um problema real, A e b devem ser preenchidos com os
    ! coeficientes obtidos a partir das equacoes das barras.

    A = 0.0d0
    b = 0.0d0

    A(1,:) = (/  4.0d0, -1.0d0,  0.0d0,  0.0d0,  0.0d0 /)
    A(2,:) = (/ -1.0d0,  4.0d0, -1.0d0,  0.0d0,  0.0d0 /)
    A(3,:) = (/  0.0d0, -1.0d0,  4.0d0, -1.0d0,  0.0d0 /)
    A(4,:) = (/  0.0d0,  0.0d0, -1.0d0,  4.0d0, -1.0d0 /)
    A(5,:) = (/  0.0d0,  0.0d0,  0.0d0, -1.0d0,  3.0d0 /)

    b = (/ 2.0d0, 4.0d0, 6.0d0, 8.0d0, 11.0d0 /)

    ! Guardamos o sistema original para verificar a resposta
    ! depois da eliminacao.
    !
    ! IMPORTANTE:
    ! A e b serao modificados durante a eliminacao. Por isso,
    ! precisamos preservar suas copias originais.

    A_original = A
    b_original = b

    ! Inicializa o contador de operacoes.
    operacoes = 0

    print *, "=============================================="
    print *, " ELIMINACAO GAUSSIANA - SISTEMA DE BARRAS"
    print *, "=============================================="

    print *
    print *, "Matriz A original:"
    do i = 1, n
        write(*,'(5(F10.3,1X))') A(i,:)
    end do

    print *
    print *, "Vetor b original:"
    write(*,'(5(F10.3,1X))') b

    ! ------------------------------------------------------------
    ! 3. ELIMINACAO PROGRESSIVA COM PIVOTAMENTO PARCIAL
    ! ------------------------------------------------------------

    ! O laco percorre as colunas de 1 ate n-1.
    !
    ! Em cada etapa:
    ! 1. Encontramos o melhor pivo disponivel.
    ! 2. Trocamos as linhas, se necessario.
    ! 3. Zeramos os elementos abaixo do pivo.
    !
    ! O pivo e o elemento da diagonal utilizado para eliminar
    ! os elementos abaixo dele.

    do k = 1, n-1

        ! --------------------------------------------------------
        ! 3.1. ENCONTRAR A LINHA DO PIVO
        ! --------------------------------------------------------

        ! Inicialmente, consideramos o elemento da diagonal
        ! como o maior candidato a pivo.
        linha_pivo = k
        maior_valor = abs(A(k,k))

        ! Procuramos o maior valor absoluto na coluna k,
        ! considerando somente as linhas ainda nao eliminadas.
        do i = k+1, n

            if (abs(A(i,k)) > maior_valor) then
                maior_valor = abs(A(i,k))
                linha_pivo = i
            end if

        end do

        ! Se o maior candidato for muito pequeno, o sistema pode
        ! ser singular ou numericamente inadequado para esta
        ! eliminacao.
        if (maior_valor < tolerancia) then
            print *
            print *, "ERRO: pivo muito pequeno."
            print *, "Nao foi possivel continuar a eliminacao."
            stop
        end if

        ! --------------------------------------------------------
        ! 3.2. TROCAR AS LINHAS, SE NECESSARIO
        ! --------------------------------------------------------

        ! Se a melhor linha nao for a linha atual, trocamos
        ! as duas linhas da matriz e os termos independentes.
        !
        ! Nao basta trocar apenas A: o vetor b tambem precisa
        ! acompanhar a troca das equacoes.

        if (linha_pivo /= k) then

            do j = 1, n

                valor_temporario = A(k,j)
                A(k,j) = A(linha_pivo,j)
                A(linha_pivo,j) = valor_temporario

            end do

            valor_temporario = b(k)
            b(k) = b(linha_pivo)
            b(linha_pivo) = valor_temporario

        end if

        ! Verificacao adicional do pivo usado na divisao.
        if (abs(A(k,k)) < tolerancia) then
            print *
            print *, "ERRO: divisao por pivo muito pequeno."
            stop
        end if

        ! --------------------------------------------------------
        ! 3.3. ELIMINAR OS ELEMENTOS ABAIXO DO PIVO
        ! --------------------------------------------------------

        do i = k+1, n

            ! O fator determina quanto da linha do pivo deve
            ! ser subtraido da linha atual.
            !
            ! Formula:
            ! fator = A(i,k) / A(k,k)

            fator = A(i,k) / A(k,k)
            operacoes = operacoes + 1

            ! Atualizamos os elementos restantes da linha.
            do j = k+1, n

                A(i,j) = A(i,j) - fator*A(k,j)

                ! Uma multiplicacao e uma subtracao.
                operacoes = operacoes + 2

            end do

            ! Atualizamos o termo independente da mesma linha.
            b(i) = b(i) - fator*b(k)
            operacoes = operacoes + 2

            ! O elemento abaixo do pivo deve ser zero.
            A(i,k) = 0.0d0

        end do

        ! Mostra o resultado parcial de cada etapa.
        print *
        print *, "Matriz apos a etapa ", k

        do i = 1, n
            write(*,'(5(F10.3,1X))') A(i,:)
        end do

    end do

    ! Verificamos o ultimo elemento diagonal antes de iniciar
    ! a substituicao regressiva.
    if (abs(A(n,n)) < tolerancia) then
        print *
        print *, "ERRO: ultimo pivo muito pequeno."
        print *, "Sistema singular ou numericamente inadequado."
        stop
    end if

    ! ------------------------------------------------------------
    ! 4. SUBSTITUICAO REGRESSIVA
    ! ------------------------------------------------------------

    ! Depois da eliminacao, temos um sistema triangular superior.
    !
    ! Calculamos primeiro a ultima incognita:
    !
    ! x(n) = b(n) / A(n,n)
    !
    ! Em seguida, calculamos as demais incognitas de baixo para
    ! cima, utilizando os valores ja encontrados.

    x(n) = b(n) / A(n,n)
    operacoes = operacoes + 1

    do i = n-1, 1, -1

        soma = 0.0d0

        ! Soma os produtos dos coeficientes pelas incognitas
        ! que ja foram calculadas.
        do j = i+1, n

            soma = soma + A(i,j)*x(j)
            operacoes = operacoes + 2

        end do

        ! Isolamos a incognita da linha atual.
        x(i) = (b(i) - soma) / A(i,i)
        operacoes = operacoes + 2

    end do

    ! ------------------------------------------------------------
    ! 5. APRESENTAR A SOLUCAO
    ! ------------------------------------------------------------

    print *
    print *, "=============================================="
    print *, " SOLUCAO DO SISTEMA"
    print *, "=============================================="

    do i = 1, n
        write(*,'("x(",I1,") = ",F12.6)') i, x(i)
    end do

    ! ------------------------------------------------------------
    ! 6. CALCULAR OS RESIDUOS
    ! ------------------------------------------------------------

    ! Para verificar a solucao, usamos o sistema ORIGINAL:
    !
    ! residuo(i) = b_original(i) - A_original(i,:)*x
    !
    ! Se a solucao estiver correta, os residuos devem ser
    ! proximos de zero, considerando a precisao numerica.

    maior_residuo = 0.0d0

    print *
    print *, "=============================================="
    print *, " VERIFICACAO DOS RESIDUOS"
    print *, "=============================================="

    do i = 1, n

        soma = 0.0d0

        do j = 1, n
            soma = soma + A_original(i,j)*x(j)
        end do

        residuo(i) = b_original(i) - soma

        write(*,'("Residuo(",I1,") = ",ES12.4)') i, residuo(i)

        if (abs(residuo(i)) > maior_residuo) then
            maior_residuo = abs(residuo(i))
        end if

    end do

    print *
    print *, "Maior residuo absoluto: ", maior_residuo

    ! ------------------------------------------------------------
    ! 7. CONTAGEM APROXIMADA DE OPERACOES
    ! ------------------------------------------------------------

    ! Este contador inclui as operacoes aritmeticas principais
    ! da eliminacao e da substituicao regressiva.
    !
    ! Nao inclui todas as operacoes da busca pelo pivo, das
    ! trocas de linhas, da verificacao dos residuos ou da impressao.
    !
    ! Portanto, e uma contagem aproximada, nao uma medicao
    ! completa do custo computacional.

    print *
    print *, "Operacoes aritmeticas contabilizadas: ", operacoes

    print *
    print *, "Fim do programa."

end program eliminacao_gaussiana_barras
