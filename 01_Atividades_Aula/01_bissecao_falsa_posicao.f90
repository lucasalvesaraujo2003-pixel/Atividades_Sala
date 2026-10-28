! =====================================================================
! ATIVIDADE 1 - METODOS DA BISSECAO E DA FALSA POSICAO
! Disciplina: Calculo Numerico Aplicado
!
! OBJETIVO:
! Encontrar a profundidade H da agua em um canal retangular aberto
! usando dois metodos numericos:
!   1) Metodo da bissecao;
!   2) Metodo da falsa posicao.
!
! Depois, calcular a velocidade media U = Q/(B*H), comparar os metodos
! e imprimir o historico das iteracoes de cada um.
!
! COMO ESTUDAR ESTE CODIGO:
! - Leia os comentarios antes de cada bloco de comandos.
! - Para adaptar o programa a outro problema, normalmente sera necessario
!   mudar os dados de entrada, a funcao f(H), o intervalo inicial e,
!   quando aplicavel, a formula usada para calcular alguma grandeza final.
! - Os dois metodos compartilham a mesma funcao, mas possuem formulas
!   diferentes para escolher a nova aproximacao da raiz.
! =====================================================================

program atividade_01_bissecao_falsa_posicao

    ! IMPLICIT NONE obriga a declarar todas as variaveis.
    ! Isso ajuda a evitar erros causados por nomes digitados incorretamente.
    implicit none

    ! KIND=8 solicita variaveis reais com precisao maior que a real padrao
    ! na maioria dos compiladores Fortran comuns.
    ! Essa precisao e importante em calculos numericos iterativos.
    integer, parameter :: dp = kind(1.0d0)

    ! -----------------------------
    ! Dados fisicos do problema
    ! -----------------------------

    ! S: inclinacao do canal, adimensional.
    ! B: largura do canal, em metros.
    ! n_manning: coeficiente de rugosidade de Manning.
    ! Q: vazao volumetrica desejada, em m^3/s.
    real(dp), parameter :: S = 0.0002_dp
    real(dp), parameter :: B = 20.0_dp
    real(dp), parameter :: n_manning = 0.03_dp
    real(dp), parameter :: Q = 5.0_dp

    ! Hl e Hu sao os limites inferior e superior do intervalo inicial.
    ! A atividade pede que seja usado o intervalo de 0,5 m a 1,0 m.
    real(dp), parameter :: H_inferior_inicial = 0.5_dp
    real(dp), parameter :: H_superior_inicial = 1.0_dp

    ! TOL e a tolerancia usada para decidir se a raiz foi encontrada.
    ! O criterio solicitado e |f(Hr)| < 10^(-6).
    real(dp), parameter :: TOL = 1.0e-6_dp

    ! MAX_ITER impede que o programa fique preso em um laco infinito
    ! caso haja um erro no algoritmo ou o metodo nao consiga convergir.
    integer, parameter :: MAX_ITER = 1000

    ! -----------------------------
    ! Variaveis da bissecao
    ! -----------------------------

    ! Cada metodo tem seus proprios limites e resultados. Isso e importante:
    ! a bissecao e a falsa posicao devem comecar do mesmo intervalo inicial,
    ! mas nao devem compartilhar limites que sao atualizados separadamente.
    real(dp) :: Hl_bis, Hu_bis
    real(dp) :: fHl_bis, fHu_bis
    real(dp) :: Hr_bis, fHr_bis
    real(dp) :: velocidade_bis
    integer :: iter_bis
    logical :: convergiu_bis

    ! -----------------------------
    ! Variaveis da falsa posicao
    ! -----------------------------

    real(dp) :: Hl_fp, Hu_fp
    real(dp) :: fHl_fp, fHu_fp
    real(dp) :: Hr_fp, fHr_fp
    real(dp) :: velocidade_fp
    real(dp) :: denominador
    integer :: iter_fp
    logical :: convergiu_fp

    ! Variaveis auxiliares para comparar os resultados.
    real(dp) :: diferenca_H, diferenca_U

    ! -----------------------------------------------------------------
    ! PARTE 1 - MOSTRAR OS DADOS E VERIFICAR O INTERVALO INICIAL
    ! -----------------------------------------------------------------

    ! A equacao de Manning para um canal retangular pode ser rearranjada
    ! para encontrar a profundidade H resolvendo f(H) = 0:
    !
    ! f(H) = (sqrt(S)/n) * (B*H)^(5/3) / (B + 2*H)^(2/3) - Q
    !
    ! Quando f(H) = 0, a vazao calculada pela equacao de Manning coincide
    ! com a vazao Q fornecida no enunciado.
    ! A funcao sera implementada uma unica vez, no final deste arquivo,
    ! e chamada pelos dois metodos. Isso evita duplicar a mesma formula.

    print *, '=========================================================='
    print *, 'ATIVIDADE 1: BISSECAO E FALSA POSICAO'
    print *, '=========================================================='
    print *, 'Dados do canal:'
    print '(A, F10.6)', ' Inclinacao S = ', S
    print '(A, F10.4)', ' Largura B (m) = ', B
    print '(A, F10.4)', ' Rugosidade n = ', n_manning
    print '(A, F10.4)', ' Vazao Q (m^3/s) = ', Q
    print *, ' '

    ! Avaliamos a funcao nos dois extremos antes de iterar.
    ! Para os metodos intervalares, precisamos de uma mudanca de sinal:
    ! f(Hl)*f(Hu) < 0.
    ! Se isso ocorrer e a funcao for continua no intervalo, existe pelo
    ! menos uma raiz entre Hl e Hu.
    print *, 'VERIFICACAO DO INTERVALO INICIAL'
    print '(A, ES16.8)', ' f(0,5) = ', funcao(H_inferior_inicial)
    print '(A, ES16.8)', ' f(1,0) = ', funcao(H_superior_inicial)

    if (funcao(H_inferior_inicial) * funcao(H_superior_inicial) >= 0.0_dp) then
        ! Se o produto for zero ou positivo, o intervalo nao satisfaz
        ! a condicao de mudanca estrita de sinal. Paramos para evitar
        ! iniciar um metodo intervalar sem a verificacao necessaria.
        print *, 'ERRO: o intervalo inicial nao apresenta mudanca de sinal.'
        print *, 'Escolha outro intervalo antes de continuar.'
        stop
    else
        print *, 'Mudanca de sinal confirmada: o intervalo e adequado.'
    end if

    print *, ' '

    ! =================================================================
    ! PARTE 2 - METODO DA BISSECAO
    ! =================================================================

    ! O metodo da bissecao sempre escolhe o ponto medio:
    !
    !              Hr = (Hl + Hu)/2
    !
    ! Depois calcula f(Hr) e substitui o extremo que possui o mesmo sinal
    ! de f(Hr). Dessa maneira, o intervalo continua contendo a raiz.

    ! Inicializamos os limites da bissecao com o intervalo do enunciado.
    Hl_bis = H_inferior_inicial
    Hu_bis = H_superior_inicial

    ! Calculamos f nos extremos. Guardar esses valores evita recalcular
    ! a funcao em cada iteracao quando apenas um extremo foi atualizado.
    fHl_bis = funcao(Hl_bis)
    fHu_bis = funcao(Hu_bis)

    ! O contador comeca em zero; a variavel logica indica se ja convergiu.
    iter_bis = 0
    convergiu_bis = .false.

    print *, '=========================================================='
    print *, 'METODO DA BISSECAO'
    print *, '=========================================================='
    print '(A)', 'Iteracao        Hl (m)          Hu (m)          Hr (m)          f(Hr)'
    print '(A)', '------------------------------------------------------------------------'

    ! DO WHILE repete os comandos enquanto ambas as condicoes forem validas:
    ! ainda nao atingimos MAX_ITER e ainda nao convergimos.
    ! .and. significa E; .not. significa NAO.
    do while (iter_bis < MAX_ITER .and. .not. convergiu_bis)

        ! Conta a iteracao atual. A primeira iteracao sera a de numero 1.
        iter_bis = iter_bis + 1

        ! Formula central da bissecao: media aritmetica dos extremos.
        Hr_bis = (Hl_bis + Hu_bis) / 2.0_dp

        ! Calcula o valor da funcao no ponto medio.
        fHr_bis = funcao(Hr_bis)

        ! Imprime o intervalo usado nesta iteracao, a aproximacao e o
        ! valor da funcao. ES16.8 mostra notacao cientifica com oito casas.
        print '(I6, 4(2X, ES16.8))', iter_bis, Hl_bis, Hu_bis, Hr_bis, fHr_bis

        ! O criterio solicitado e o modulo do residuo da funcao menor que
        ! a tolerancia. ABS(x) devolve o valor absoluto de x.
        if (abs(fHr_bis) < TOL) then
            convergiu_bis = .true.
        else
            ! Se f(Hl) e f(Hr) possuem sinais opostos, a raiz esta entre
            ! Hl e Hr; portanto, o limite superior passa a ser Hr.
            if (fHl_bis * fHr_bis < 0.0_dp) then
                Hu_bis = Hr_bis
                fHu_bis = fHr_bis
            else
                ! Caso contrario, a mudanca de sinal esta entre Hr e Hu.
                ! Assim, o limite inferior passa a ser Hr.
                Hl_bis = Hr_bis
                fHl_bis = fHr_bis
            end if
        end if

    end do

    ! Depois do laco, Hr_bis guarda a ultima aproximacao calculada.
    ! Se convergiu_bis for falso, o limite de iteracoes foi atingido.
    if (convergiu_bis) then
        print *, 'Bissecao convergiu com sucesso.'
    else
        print *, 'AVISO: a bissecao atingiu o limite de iteracoes.'
    end if

    print '(A, F16.10, A)', ' Profundidade H = ', Hr_bis, ' m'
    print '(A, I8)', ' Numero de iteracoes = ', iter_bis
    print '(A, ES16.8)', ' Residuo |f(H)| = ', abs(fHr_bis)

    ! Para um canal retangular, a area molhada e A = B*H.
    ! Como Q = U*A, isolamos a velocidade media: U = Q/(B*H).
    velocidade_bis = Q / (B * Hr_bis)
    print '(A, F16.10, A)', ' Velocidade media U = ', velocidade_bis, ' m/s'
    print *, ' '

    ! =================================================================
    ! PARTE 3 - METODO DA FALSA POSICAO
    ! =================================================================

    ! A falsa posicao tambem conserva um intervalo com mudanca de sinal,
    ! mas escolhe a nova aproximacao pela intersecao com o eixo horizontal
    ! da reta que liga (Hl, f(Hl)) e (Hu, f(Hu)):
    !
    ! Hr = Hu - f(Hu)*(Hu - Hl)/(f(Hu) - f(Hl))
    !
    ! Forma equivalente:
    ! Hr = [Hl*f(Hu) - Hu*f(Hl)]/[f(Hu) - f(Hl)]
    !
    ! Esta atividade usa novamente o intervalo original, e nao o intervalo
    ! final que sobrou depois da bissecao.

    Hl_fp = H_inferior_inicial
    Hu_fp = H_superior_inicial
    fHl_fp = funcao(Hl_fp)
    fHu_fp = funcao(Hu_fp)

    iter_fp = 0
    convergiu_fp = .false.

    print *, '=========================================================='
    print *, 'METODO DA FALSA POSICAO'
    print *, '=========================================================='
    print '(A)', 'Iteracao        Hl (m)          Hu (m)          Hr (m)          f(Hr)'
    print '(A)', '------------------------------------------------------------------------'

    do while (iter_fp < MAX_ITER .and. .not. convergiu_fp)

        iter_fp = iter_fp + 1

        ! O denominador e f(Hu)-f(Hl). Se for zero, a formula ficaria
        ! indefinida. TINY(denominador) fornece um numero positivo muito
        ! pequeno do mesmo tipo numerico, para proteger contra divisao por
        ! zero ou por um valor praticamente nulo.
        denominador = fHu_fp - fHl_fp

        if (abs(denominador) <= tiny(denominador)) then
            print *, 'ERRO: denominador quase nulo na falsa posicao.'
            exit
        end if

        ! Formula da falsa posicao: intersecao da reta com o eixo H.
        Hr_fp = Hu_fp - fHu_fp * (Hu_fp - Hl_fp) / denominador

        ! Avalia a funcao na nova estimativa.
        fHr_fp = funcao(Hr_fp)

        ! Mostra o historico para permitir acompanhar a convergencia.
        print '(I6, 4(2X, ES16.8))', iter_fp, Hl_fp, Hu_fp, Hr_fp, fHr_fp

        ! Usa o mesmo criterio de convergencia da bissecao:
        ! o modulo de f(Hr) deve ser menor que 10^(-6).
        if (abs(fHr_fp) < TOL) then
            convergiu_fp = .true.
        else
            ! Atualiza um extremo, mantendo a mudanca de sinal no intervalo.
            if (fHl_fp * fHr_fp < 0.0_dp) then
                Hu_fp = Hr_fp
                fHu_fp = fHr_fp
            else
                Hl_fp = Hr_fp
                fHl_fp = fHr_fp
            end if
        end if

    end do

    if (convergiu_fp) then
        print *, 'Falsa posicao convergiu com sucesso.'
    else
        print *, 'AVISO: a falsa posicao nao confirmou convergencia.'
    end if

    print '(A, F16.10, A)', ' Profundidade H = ', Hr_fp, ' m'
    print '(A, I8)', ' Numero de iteracoes = ', iter_fp
    print '(A, ES16.8)', ' Residuo |f(H)| = ', abs(fHr_fp)

    velocidade_fp = Q / (B * Hr_fp)
    print '(A, F16.10, A)', ' Velocidade media U = ', velocidade_fp, ' m/s'
    print *, ' '

    ! =================================================================
    ! PARTE 4 - COMPARACAO DOS DOIS METODOS
    ! =================================================================

    print *, '=========================================================='
    print *, 'COMPARACAO FINAL'
    print *, '=========================================================='

    print '(A, F16.10, A)', ' Profundidade pela bissecao = ', Hr_bis, ' m'
    print '(A, F16.10, A)', ' Profundidade pela falsa posicao = ', Hr_fp, ' m'
    print '(A, I8)', ' Iteracoes da bissecao = ', iter_bis
    print '(A, I8)', ' Iteracoes da falsa posicao = ', iter_fp
    print '(A, F16.10, A)', ' Velocidade pela bissecao = ', velocidade_bis, ' m/s'
    print '(A, F16.10, A)', ' Velocidade pela falsa posicao = ', velocidade_fp, ' m/s'

    ! Diferencas absolutas ajudam a verificar se os resultados dos metodos
    ! estao proximos. ABS(a-b) calcula a distancia positiva entre os valores.
    diferenca_H = abs(Hr_bis - Hr_fp)
    diferenca_U = abs(velocidade_bis - velocidade_fp)

    print '(A, ES16.8, A)', ' Diferenca entre profundidades = ', diferenca_H, ' m'
    print '(A, ES16.8, A)', ' Diferenca entre velocidades = ', diferenca_U, ' m/s'

    ! -----------------------------------------------------------------
    ! RESPOSTAS CONCEITUAIS PARA REVISAR ANTES DA PROVA
    ! -----------------------------------------------------------------
    ! a) Os dois metodos devem encontrar profundidades muito proximas,
    !    pois resolvem a mesma equacao e usam o mesmo intervalo inicial.
    !    Pequenas diferencas aparecem porque cada metodo gera aproximacoes
    !    diferentes e para quando |f(Hr)| < TOL.
    !
    ! b) Compare os contadores impressos ao final para identificar qual
    !    metodo precisou de menos iteracoes neste problema. Nao se deve
    !    assumir que a falsa posicao sempre sera mais rapida: isso depende
    !    da funcao e do intervalo.
    !
    ! c) A bissecao escolhe sempre o ponto medio do intervalo. A falsa
    !    posicao usa a intersecao com o eixo horizontal da reta que liga
    !    os pontos (Hl,f(Hl)) e (Hu,f(Hu)).
    !
    ! d) Ambos sao metodos intervalares porque mantem um intervalo que
    !    contem a raiz e reduzem esse intervalo a cada iteracao.
    !
    ! e) A condicao f(Hl)*f(Hu) < 0 indica que os extremos tem sinais
    !    opostos. Para uma funcao continua, isso garante pelo menos uma
    !    raiz no intervalo. A atualizacao deve preservar essa mudanca.
    !
    ! f) As velocidades sao calculadas por U=Q/(B*H). Como Q e B sao iguais
    !    nos dois calculos, profundidades proximas produzem velocidades
    !    proximas.
    !
    ! -----------------------------------------------------------------
    ! COMO REUTILIZAR PARTES DESTE CODIGO EM OUTRAS QUESTOES
    ! -----------------------------------------------------------------
    !
    ! ESQUELETO DA BISSECAO:
    ! 1. Defina Hl e Hu.
    ! 2. Calcule f(Hl) e f(Hu).
    ! 3. Verifique se f(Hl)*f(Hu) < 0.
    ! 4. Inicialize o contador e o indicador de convergencia.
    ! 5. Dentro de um DO WHILE, incremente o contador.
    ! 6. Calcule Hr=(Hl+Hu)/2.
    ! 7. Calcule f(Hr) e verifique o criterio de parada.
    ! 8. Atualize Hl ou Hu para preservar a mudanca de sinal.
    ! 9. Imprima a aproximacao e o numero de iteracoes.
    !
    ! ESQUELETO DA FALSA POSICAO:
    ! Use os mesmos passos de verificacao, laco, criterio de parada e
    ! atualizacao de intervalo da bissecao. Troque somente a formula de Hr
    ! pela formula da falsa posicao e confira o denominador antes de dividir.
    !
    ! PARA TROCAR O PROBLEMA FISICO:
    ! - altere os parametros declarados no inicio do programa;
    ! - altere a funcao funcao(H), no final do arquivo;
    ! - escolha um intervalo que apresente mudanca de sinal;
    ! - altere o calculo final se a grandeza desejada for outra.
    !
    ! SOBRE CRITERIOS DE PARADA:
    ! abs(f(Hr)) < TOL verifica o residuo da funcao, como pede esta atividade.
    ! abs(Hr_novo-Hr_antigo) < TOL verifica a mudanca entre aproximacoes.
    ! Esses criterios nao sao identicos; use aquele pedido no enunciado.
    !
    ! DO WHILE, IF, contador, variaveis logicas e funcao separada sao
    ! estruturas que voltarao a aparecer nas proximas atividades.
    ! =====================================================================

contains

    ! -----------------------------------------------------------------
    ! FUNCAO DO PROBLEMA: f(H)
    ! -----------------------------------------------------------------
    ! CONTAINS permite declarar uma funcao interna ao programa principal.
    ! Assim, tanto a bissecao quanto a falsa posicao podem chamar funcao(H)
    ! sem copiar e colar a formula em dois lugares diferentes.
    !
    ! Entrada:
    !   H - profundidade da agua, em metros.
    !
    ! Saida:
    !   valor de f(H). A raiz procurada e o H para o qual f(H)=0.
    !
    ! A formula geral vem da equacao de Manning para canal retangular:
    ! Q_calc = (sqrt(S)/n) * (B*H)^(5/3) / (B+2H)^(2/3)
    ! f(H) = Q_calc - Q
    !
    ! Se outro exercicio tiver uma equacao diferente, esta e a parte
    ! principal a substituir, mantendo os metodos numericos se forem
    ! aplicaveis a nova funcao.
    function funcao(H) result(valor)

        ! INTENT(IN) informa que H e apenas entrada: a funcao nao altera H.
        real(dp), intent(in) :: H

        ! valor e o numero devolvido pela funcao ao programa que a chamou.
        real(dp) :: valor

        ! Calcula a vazao prevista por Manning e subtrai a vazao especificada.
        ! Os expoentes fracionarios sao escritos com _dp para manter a
        ! precisao real usada no restante do programa.
        valor = (sqrt(S) / n_manning) * (B * H)**(5.0_dp / 3.0_dp) &
                / (B + 2.0_dp * H)**(2.0_dp / 3.0_dp) - Q

    end function funcao

end program atividade_01_bissecao_falsa_posicao