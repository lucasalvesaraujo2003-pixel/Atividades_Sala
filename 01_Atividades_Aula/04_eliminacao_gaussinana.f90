! =====================================================================
! ATIVIDADE 4 - ELIMINACAO GAUSSIANA
! APLICACAO: REDE DE REATORES
!
! OBJETIVO:
! Resolver um sistema linear de 5 equacoes e 5 incognitas utilizando
! o metodo da eliminacao gaussiana.
!
! MODELO MATEMATICO:
!
! A * x = b
!
! A = matriz dos coeficientes
! x = vetor das incognitas
! b = vetor dos termos independentes
!
! ETAPAS DO METODO:
!
! 1. Montagem da matriz A e do vetor b.
! 2. Eliminacao progressiva para obter uma matriz triangular superior.
! 3. Substituicao retroativa para calcular as incognitas.
! 4. Calculo dos residuos para verificar a solucao.
! 5. Contagem aproximada das operacoes aritmeticas (FLOPs).
!
! IMPORTANTE:
! Os coeficientes abaixo sao ilustrativos, pois os valores numericos
! originais do enunciado nao estao disponiveis neste momento.
! Substitua A e b pelos dados fornecidos na sua atividade.
!
! Este exemplo possui solucao conhecida:
! x = [1, 2, 3, 4, 5]
!
! COMPARACAO COM OS EXERCICIOS ANTERIORES:
!
! Bisseccao e falsa posicao:
! Encontram a raiz de uma equacao usando um intervalo.
!
! Newton-Raphson e secante:
! Encontram a raiz de uma equacao por aproximacoes sucessivas.
!
! Newton-Raphson para sistemas:
! Resolve sistemas nao lineares por meio de uma matriz Jacobiana.
!
! Eliminacao gaussiana:
! Resolve sistemas lineares por operacoes nas linhas da matriz.
!
! DIFERENCA IMPORTANTE:
! Este codigo nao utiliza pivoteamento. Essa tecnica sera estudada
! na atividade seguinte.
!
! ARQUIVO: 04_eliminacao_gaussiana_reatores.md
! =====================================================================

program atividade_04

    implicit none

    ! IMPLICIT NONE:
    ! Obriga a declarar todas as variaveis.
    ! Evita erros de digitacao nos nomes das variaveis.

    ! PARAMETER:
    ! Define uma constante que nao pode ser alterada durante a execucao.
    ! Use quando o tamanho do sistema for conhecido antecipadamente.

    integer, parameter :: n = 5

    ! REAL(8):
    ! Utiliza numeros reais em precisao dupla.
    ! E util para reduzir erros de arredondamento nos calculos.

    real(8) :: A(n,n), A_original(n,n)
    real(8) :: b(n), b_original(n)
    real(8) :: x(n)
    real(8) :: fator, soma, residuo
    real(8) :: maior_residuo
    real(8) :: tolerancia_pivo

    integer :: i, j, k
    integer :: flops

    ! ================================================================
    ! 1. MONTAGEM DO SISTEMA LINEAR
    ! ================================================================

    ! A matriz A representa os coeficientes das cinco equacoes.
    !
    ! O vetor b contem os termos independentes.
    !
    ! Para adaptar o programa a outro problema, altere os valores
    ! de A e b sem modificar o algoritmo de eliminacao.
    !
    ! Os dados deste exemplo foram escolhidos para que a solucao
    ! seja conhecida e os resultados possam ser conferidos.

    A(1,:) = (/ 10.0d0, -1.0d0,  2.0d0,  0.0d0,  0.0d0 /)
    A(2,:) = (/ -1.0d0, 11.0d0, -1.0d0,  3.0d0,  0.0d0 /)
    A(3,:) = (/  2.0d0, -1.0d0, 10.0d0, -1.0d0,  2.0d0 /)
    A(4,:) = (/  0.0d0,  3.0d0, -1.0d0,  8.0d0, -1.0d0 /)
    A(5,:) = (/  0.0d0,  0.0d0,  2.0d0, -1.0d0,  5.0d0 /)

    b = (/ 14.0d0, 30.0d0, 36.0d0, 30.0d0, 27.0d0 /)

    ! ================================================================
    ! 2. ARMAZENAMENTO DOS DADOS ORIGINAIS
    ! ================================================================

    ! A eliminacao gaussiana modifica os coeficientes da matriz A
    ! e os elementos do vetor b.
    !
    ! Por isso, guardamos copias dos dados originais para calcular
    ! posteriormente os residuos da solucao.
    !
    ! O simbolo = em atribuicoes de vetores e matrizes copia os dados.

    A_original = A
    b_original = b

    ! Inicializa o contador de operacoes aritmeticas.
    ! A contagem abaixo considera operacoes escalares de soma,
    ! subtracao, multiplicacao e divisao.

    flops = 0

    ! Tolerancia para identificar pivôs nulos ou numericamente
    ! muito pequenos. O teste nao substitui o pivoteamento.
    tolerancia_pivo = 1.0d-14

    print *, "=============================================="
    print *, " ATIVIDADE 4 - ELIMINACAO GAUSSIANA"
    print *, "=============================================="
    print *, "Numero de equacoes:", n
    print *

    print *, "MATRIZ ORIGINAL A:"
    do i = 1, n
        print *, A(i,:)
    end do

    print *
    print *, "VETOR ORIGINAL b:"
    print *, b
    print *

    ! ================================================================
    ! 3. ELIMINACAO PROGRESSIVA
    ! ================================================================

    ! OBJETIVO:
    ! Transformar a matriz A em uma matriz triangular superior.
    !
    ! A ideia e eliminar os elementos abaixo da diagonal principal.
    !
    ! DO k = 1, n-1:
    ! Escolhe a coluna do pivo.
    !
    ! DO i = k+1, n:
    ! Percorre as linhas abaixo do pivo.
    !
    ! DO j = k+1, n:
    ! Atualiza os coeficientes restantes da linha.
    !
    ! Esta estrutura de tres lacos e comum em algoritmos matriciais.
    ! A ordem dos lacos importa: cada um controla uma dimensao diferente.

    do k = 1, n-1

        ! O elemento A(k,k) e o pivo da etapa atual.
        !
        ! Sem pivoteamento, um pivo nulo impede a divisao.
        ! Mesmo um pivo nao nulo pode causar erros numericos se
        ! tiver magnitude muito pequena.

        if (abs(A(k,k)) < tolerancia_pivo) then

            print *, "ERRO: pivo nulo ou muito pequeno."
            print *, "A matriz exige pivoteamento."
            print *, "A eliminacao foi interrompida."

            stop

        end if

        ! Percorre todas as linhas abaixo da linha do pivo.

        do i = k+1, n

            ! FATOR DE ELIMINACAO:
            !
            ! fator = A(i,k) / A(k,k)
            !
            ! O fator indica quanto da linha do pivo deve ser
            ! subtraido da linha atual para zerar A(i,k).
            !
            ! Contabilizamos uma divisao.

            fator = A(i,k) / A(k,k)

            flops = flops + 1

            ! Atualiza os coeficientes a direita do pivo.
            !
            ! A(i,j) = A(i,j) - fator*A(k,j)
            !
            ! Para cada j:
            ! 1 multiplicacao;
            ! 1 subtracao.
            !
            ! O elemento A(i,k) e zerado explicitamente depois.

            do j = k+1, n

                A(i,j) = A(i,j) - fator*A(k,j)

                flops = flops + 2

            end do

            ! Atualiza o termo independente da mesma equacao.
            !
            ! b(i) = b(i) - fator*b(k)
            !
            ! Sem atualizar b, o sistema transformado seria
            ! matematicamente inconsistente.

            b(i) = b(i) - fator*b(k)

            flops = flops + 2

            ! O valor abaixo da diagonal deve ser zero.
            ! Essa atribuicao representa a estrutura triangular
            ! do sistema transformado.

            A(i,k) = 0.0d0

        end do

    end do

    ! Verifica o ultimo pivo antes da substituicao retroativa.
    if (abs(A(n,n)) < tolerancia_pivo) then

        print *, "ERRO: ultimo pivo nulo ou muito pequeno."
        print *, "A eliminacao foi interrompida."

        stop

    end if

    ! ================================================================
    ! 4. APRESENTACAO DA MATRIZ TRIANGULAR SUPERIOR
    ! ================================================================

    ! A matriz agora possui zeros abaixo da diagonal principal.
    ! O vetor b tambem foi transformado durante a eliminacao.

    print *, "MATRIZ TRIANGULAR SUPERIOR:"

    do i = 1, n
        print *, A(i,:)
    end do

    print *
    print *, "VETOR b APOS A ELIMINACAO:"
    print *, b
    print *

    ! ================================================================
    ! 5. SUBSTITUICAO RETROATIVA
    ! ================================================================

    ! A matriz triangular superior permite calcular as incognitas
    ! da ultima equacao para a primeira.
    !
    ! Formula:
    !
    ! x(i) = [b(i) - SOMA(A(i,j)*x(j))] / A(i,i)
    !
    ! O laco externo percorre as linhas de baixo para cima.
    ! O laco interno utiliza as incognitas ja calculadas.
    !
    ! Em Fortran, DO i = n, 1, -1 percorre:
    ! n, n-1, n-2, ..., 1.

    do i = n, 1, -1

        soma = 0.0d0

        do j = i+1, n

            soma = soma + A(i,j)*x(j)

            ! Uma multiplicacao e uma soma por termo.
            flops = flops + 2

        end do

        ! Subtrai a soma dos termos conhecidos do lado direito.
        soma = b(i) - soma

        flops = flops + 1

        ! Calcula a incognita da linha atual.
        x(i) = soma / A(i,i)

        flops = flops + 1

    end do

    ! ================================================================
    ! 6. APRESENTACAO DA SOLUCAO
    ! ================================================================

    print *, "=============================================="
    print *, " SOLUCAO DO SISTEMA LINEAR"
    print *, "=============================================="

    do i = 1, n

        ! Escreve cada incognita e seu valor calculado.
        ! A numeracao indica a posicao no vetor de solucao.

        print *, "x(", i, ") = ", x(i)

    end do

    print *

    ! ================================================================
    ! 7. CALCULO DOS RESIDUOS
    ! ================================================================

    ! O residuo de cada equacao e definido por:
    !
    ! r(i) = SOMA(A_original(i,j)*x(j)) - b_original(i)
    !
    ! Se a solucao estiver correta, os residuos devem ser
    ! proximos de zero, considerando a precisao numerica.
    !
    ! E necessario utilizar A_original e b_original porque A e b
    ! ja foram modificados pela eliminacao gaussiana.

    maior_residuo = 0.0d0

    print *, "RESIDUOS DAS EQUACOES:"
    print *

    do i = 1, n

        residuo = 0.0d0

        do j = 1, n

            residuo = residuo + A_original(i,j)*x(j)

        end do

        residuo = residuo - b_original(i)

        print *, "Residuo da equacao ", i, ": ", residuo

        ! MAX identifica o maior valor entre dois numeros.
        ! ABS calcula o modulo, evitando que residuos negativos
        ! sejam interpretados como erros pequenos.

        maior_residuo = max(maior_residuo, abs(residuo))

    end do

    print *
    print *, "Maior residuo absoluto:", maior_residuo

    ! ================================================================
    ! 8. CONTAGEM DAS OPERACOES ARITMETICAS
    ! ================================================================

    ! FLOPs:
    ! Estimativa do numero de operacoes aritmeticas realizadas
    ! durante a eliminacao e a substituicao retroativa.
    !
    ! Esta contagem considera somas, subtracoes, multiplicacoes
    ! e divisoes explicitas no algoritmo.
    !
    ! Nao contabiliza comparacoes, atribuicoes ou operacoes usadas
    ! exclusivamente para calcular os residuos.

    print *
    print *, "Operacoes aritmeticas contabilizadas:", flops

    ! ================================================================
    ! 9. VERIFICACAO FINAL
    ! ================================================================

    if (maior_residuo < 1.0d-10) then

        print *, "Verificacao: solucao consistente."

    else

        print *, "Atencao: residuos significativos."
        print *, "Verifique os dados e a estabilidade numerica."

    end if

    print *
    print *, "Fim do programa."

end program atividade_04
