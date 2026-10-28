# PPC2 — Método de Bairstow para determinação de raízes de polinômios

**Cálculo Numérico Aplicado — ENM0227**  
**Aluno:** Lucas Alves Araujo  
**Matrícula:** 221021705  
**Professor:** Rafael Gabler Gontijo  
**Data:** 24 de setembro de 2026

---

## 1. Objetivo

Este projeto implementa computacionalmente o **método de Bairstow** para a determinação das raízes de polinômios, incluindo raízes reais e complexas conjugadas.

Além da implementação do algoritmo, o programa realiza:

1. validação com um polinômio de sétimo grau cujas raízes são conhecidas;
2. análise da convergência para diferentes valores iniciais de `r0` e `s0`;
3. aplicação do método ao polinômio característico obtido na APC2;
4. interpretação das raízes encontradas para o sistema dinâmico;
5. varredura do plano `(r0,s0)` para construção do fractal de Bairstow.

O desenvolvimento segue a proposta do Programa para Casa #2, que solicita a implementação do método sem o uso de bibliotecas que resolvam diretamente o problema de raízes.

---

## 2. Fundamentação matemática

Considera-se um polinômio escrito na forma

\[
P(x)=a_0+a_1x+a_2x^2+\cdots+a_nx^n.
\]

O método de Bairstow procura um fator quadrático da forma

\[
D(x)=x^2-rx-s.
\]

A divisão sintética é realizada por meio das relações

\[
b_n=a_n,
\]

\[
b_{n-1}=a_{n-1}+rb_n,
\]

\[
b_i=a_i+rb_{i+1}+sb_{i+2}.
\]

O restante da divisão é utilizado para determinar as correções de Newton-Raphson dos parâmetros `r` e `s`. As expressões empregadas no código seguem as relações apresentadas no material da disciplina.

Quando o resto tende a zero, o fator quadrático é considerado convergido. Suas raízes são obtidas por

\[
x_{1,2}=\frac{r\pm\sqrt{r^2+4s}}{2}.
\]

O processo é repetido sobre o polinômio deflacionado até a obtenção de todas as raízes.

---

## 3. Arquivos do projeto

```text
PPC2_Bairstow/
│
├── main.f90
├── plot_fractal.py
├── fractal_bairstow.png
├── validacao.txt
├── convergencia.txt
├── raizes_apc2.txt
├── fractal.dat
└── resumo.txt
```

O arquivo `main.f90` contém o algoritmo de Bairstow e as rotinas necessárias para a execução das quatro etapas do trabalho. O código está comentado em português para facilitar a compreensão de cada parte da implementação.

---

## 4. Implementação

A implementação principal está no arquivo `main.f90`.

A organização do código foi feita em um módulo denominado `bairstow_module`, contendo as seguintes rotinas:

- `multiply_polynomials`: realiza multiplicação de polinômios;
- `bairstow_factor`: executa uma fatoração quadrática pelo método de Bairstow;
- `bairstow_solve`: realiza a deflação sucessiva até encontrar todas as raízes;
- `analyze_initial_pair`: testa diretamente uma escolha de `r0` e `s0`;
- `evaluate_polynomial`: avalia um polinômio pelo método de Horner.

A rotina principal `ppc2_bairstow` executa os testes e gera os arquivos de saída.

---

## 5. Validação com polinômio de sétimo grau

Para validar o método, foi utilizado o polinômio

\[
P(x)=
26.304
-6.272x
-10.168x^2
+2.668x^3
-6.272x^4
-7.46x^5
+0.2x^6
+x^7.
\]

Esse polinômio foi construído a partir dos fatores

\[
(x-1)(x+2)(x-3)
(x^2-0.8x+1.6)
(x^2+3x+2.74).
\]

Portanto, as raízes conhecidas são

\[
1,
\qquad
-2,
\qquad
3,
\]

\[
0.4\pm1.2i,
\]

\[
-1.5\pm0.7i.
\]

O programa recuperou as raízes com erro compatível com a precisão numérica adotada. A execução apresentou convergência e as raízes complexas apareceram corretamente em pares conjugados.

### Resultado obtido

```text
1.0000000000 + 0.0000000000 i
-2.0000000000 + 0.0000000000 i
-1.5000000000 + 0.7000000000 i
-1.5000000000 - 0.7000000000 i
0.4000000000 + 1.1999999999 i
0.4000000000 - 1.1999999999 i
3.0000000000 + 0.0000000000 i
```

Foram necessárias **38 iterações totais** para as sucessivas fatorações do polinômio de sétimo grau.

Esse teste atende à etapa de validação solicitada no PPC2, que pede a construção de um polinômio de sétima ordem com raízes previamente conhecidas e a verificação da recuperação dessas raízes.

---

## 6. Polinômio característico da APC2

Na APC2 foi adotado um sistema com dois graus de liberdade e foram utilizadas as matrizes

\[
M=
\begin{bmatrix}
2 & 0\\
0 & 1.5
\end{bmatrix},
\]

\[
C=
\begin{bmatrix}
0.8 & -0.2\\
-0.2 & 0.5
\end{bmatrix},
\]

\[
K=
\begin{bmatrix}
20 & -5\\
-5 & 12
\end{bmatrix}.
\]

O polinômio característico obtido na APC2 é

\[
\boxed{
P(\lambda)=
3\lambda^4
+2.2\lambda^3
+54.36\lambda^2
+17.6\lambda
+215
}.
\]

Esse polinômio possui os coeficientes

\[
(a_4,a_3,a_2,a_1,a_0)
=
(3,2.2,54.36,17.6,215).
\]

---

## 7. Raízes obtidas para a APC2

Aplicando o método de Bairstow ao polinômio da APC2, foram encontradas as seguintes raízes:

| Raiz | Parte real | Parte imaginária |
|---|---:|---:|
| $\lambda_1$ | $-0.1233234561$ | $+2.4351138561$ |
| $\lambda_2$ | $-0.1233234561$ | $-2.4351138561$ |
| $\lambda_3$ | $-0.2433432106$ | $+3.4634890229$ |
| $\lambda_4$ | $-0.2433432106$ | $-3.4634890229$ |

A execução apresentou convergência em **7 iterações** para a fatoração do polinômio de quarto grau.

As raízes são complexas conjugadas e todas possuem parte real negativa. Escrevendo

\[
\lambda=\alpha+i\beta,
\]

as partes reais negativas indicam decaimento exponencial da resposta, enquanto as partes imaginárias diferentes de zero indicam comportamento oscilatório.

A interpretação física das raízes é justamente uma das etapas solicitadas no PPC2.

---

## 8. Frequências associadas às raízes

Para uma raiz complexa

\[
\lambda=\alpha+i\beta,
\]

a componente imaginária pode ser associada à frequência angular de oscilação.

Para as duas famílias de raízes encontradas:

\[
\omega_1 \approx 2.4351\;\text{rad/s},
\]

\[
\omega_2 \approx 3.4635\;\text{rad/s}.
\]

As frequências em hertz são aproximadamente

\[
f_1=\frac{2.4351}{2\pi}\approx0.3876\;\text{Hz},
\]

\[
f_2=\frac{3.4635}{2\pi}\approx0.5512\;\text{Hz}.
\]

Como as partes reais são negativas, as oscilações apresentam decaimento ao longo do tempo.

---

## 9. Análise da convergência

Foram testados diferentes pares iniciais para `r0` e `s0`. Os resultados obtidos para o primeiro processo de fatoração do polinômio da APC2 foram:

| $r_0$ | $s_0$ | Iterações | Convergiu | Resto |
|---:|---:|---:|:---:|---:|
| 0.0 | 0.0 | 7 | Sim | $2.84\times10^{-14}$ |
| 1.0 | 0.0 | 7 | Sim | $5.88\times10^{-11}$ |
| -1.0 | 0.0 | 7 | Sim | $2.84\times10^{-14}$ |
| 0.0 | 1.0 | 7 | Sim | $1.78\times10^{-15}$ |
| 1.0 | 1.0 | 7 | Sim | $4.53\times10^{-10}$ |
| -1.0 | 1.0 | 7 | Sim | $1.78\times10^{-15}$ |
| 2.0 | -1.0 | 8 | Sim | $2.84\times10^{-14}$ |
| -2.0 | 1.0 | 7 | Sim | $7.35\times10^{-9}$ |

Neste conjunto de testes, todas as escolhas analisadas convergiram. A maior quantidade observada foi de **8 iterações**, enquanto as demais escolhas convergiram em 7 iterações.

O resultado mostra que, para este polinômio específico e para as sementes testadas, o algoritmo apresentou convergência rápida. Isso não significa que todas as escolhas possíveis de `r0` e `s0` necessariamente apresentem o mesmo comportamento, pois a própria formulação do método admite dependência da escolha inicial.

---

## 10. Fractal de Bairstow

Para analisar de forma sistemática a influência dos valores iniciais, foi realizada uma varredura do plano

\[
r_0\in[-3,3],
\qquad
s_0\in[-3,3],
\]

usando uma malha de **201 × 201 pontos**, totalizando **40.401 combinações**.

Para cada ponto foi registrado o número de iterações necessárias para que a primeira fatoração quadrática atingisse a tolerância estabelecida.

Neste caso específico, todos os **40.401 pontos** analisados convergiram. O número de iterações observado variou entre **6 e 9**.

A figura abaixo apresenta o mapa obtido para o número de iterações:

![Fractal de Bairstow](fractal_bairstow.png)

O arquivo `fractal.dat` contém os dados utilizados para a construção da figura. O enunciado do PPC2 solicita justamente uma varredura dos valores iniciais `r0` e `s0` e a construção do fractal associado ao comportamento de convergência do método.

Como todos os pontos analisados convergiram neste experimento, não foram observadas regiões de divergência dentro do intervalo testado. Ainda assim, o mapa permite visualizar a variação do número de iterações conforme a condição inicial é alterada.

---

## 11. Tolerância e critérios de convergência

Foi adotada a tolerância

\[
\varepsilon=10^{-10}.
\]

A convergência da fatoração quadrática é considerada atingida quando o resto da divisão, representado pelos termos `b0` e `b1`, apresenta magnitude suficientemente pequena em relação à escala dos coeficientes do polinômio.

Foi utilizado um limite máximo de **200 iterações por tentativa** para impedir loops indefinidos em casos de não convergência.

Também foi incluída uma verificação do determinante do sistema linear de Newton. Quando esse determinante se aproxima de zero, a tentativa é considerada numericamente inadequada e outra semente inicial pode ser utilizada pelo solucionador geral.

---

## 12. Como compilar e executar

### Fortran

Com o compilador `gfortran` instalado, execute:

```bash
gfortran -std=f2008 -Wall -Wextra -O2 -o ppc2 main.f90
```

Depois:

```bash
./ppc2
```

No Windows, utilizando o `gfortran` pelo terminal:

```powershell
gfortran -std=f2008 -Wall -Wextra -O2 -o ppc2.exe main.f90
```

e então:

```powershell
./ppc2.exe
```

Após a execução serão gerados:

```text
validacao.txt
convergencia.txt
raizes_apc2.txt
fractal.dat
resumo.txt
```

### Figura do fractal

Para gerar a figura a partir do arquivo `fractal.dat`, execute:

```bash
python plot_fractal.py
```

Será criado o arquivo:

```text
fractal_bairstow.png
```

---

## 13. Conclusão

O método de Bairstow foi implementado utilizando apenas operações e estruturas básicas da linguagem Fortran, permitindo a determinação de raízes reais e complexas de polinômios.

A validação com um polinômio de sétimo grau mostrou que o algoritmo foi capaz de recuperar as raízes previamente conhecidas com elevada precisão numérica.

Na aplicação ao polinômio característico da APC2, foram obtidas quatro raízes complexas conjugadas com partes reais negativas, caracterizando respostas oscilatórias amortecidas para o modelo estudado.

A análise de diferentes valores iniciais de `r0` e `s0` mostrou convergência rápida para os casos avaliados. A varredura de 40.401 pontos no plano de condições iniciais também permitiu gerar o mapa de iterações associado ao comportamento do método.

Assim, a implementação contempla as principais etapas solicitadas no PPC2: implementação do método, validação, análise de convergência, aplicação à APC2 e construção do mapa associado às condições iniciais.

---

## 14. Referências

1. CHAPRA, S. C.; CANALE, R. P. **Métodos Numéricos para Engenharia**. 5. ed. McGraw-Hill, 2008.

2. GONTIJO, R. G. **Notas de aula do curso de Cálculo Numérico Aplicado**. 2026.

---

## 11. Código completo em Fortran

O código abaixo é o mesmo utilizado na execução do projeto. Ele pode ser copiado integralmente para um arquivo chamado `main.f90` e compilado com um compilador Fortran compatível com Fortran 2008.

> **Observação:** o bloco abaixo é exibido pelo README para documentação do projeto. O GitHub não executa código Fortran diretamente dentro do arquivo `.md`; para execução, o mesmo conteúdo deve estar salvo como `main.f90` no repositório.

```fortran
! ============================================================================
! PPC2 - Calculo Numerico Aplicado - ENM0227
! Metodo de Bairstow para determinacao de raizes de polinomios
!
! Aluno: Lucas Alves Araujo
! Matricula: 221021705
! Professor: Rafael Gabler Gontijo
!
! Convencao adotada:
!   P(x) = a(0) + a(1)x + ... + a(n)x^n
!
! O codigo implementa o metodo de Bairstow usando somente estruturas basicas
! da linguagem Fortran. Nao sao utilizadas bibliotecas de busca de raizes.
!
! O programa realiza quatro tarefas principais:
!   1) valida o algoritmo com um polinomio de 7o grau de raizes conhecidas;
!   2) analisa a convergencia para diferentes valores iniciais de r0 e s0;
!   3) resolve o polinomio caracteristico obtido na APC2;
!   4) gera uma varredura para o fractal de Bairstow.
! ============================================================================

module bairstow_module
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    real(dp), parameter :: TOL = 1.0d-10
    integer, parameter :: MAX_ITER = 200
    real(dp), parameter :: DET_TOL = 1.0d-14

contains

    ! ------------------------------------------------------------------------
    ! Multiplica dois polinomios escritos em ordem crescente de potencia:
    ! p(x) = p(0) + p(1)x + ...
    ! ------------------------------------------------------------------------
    subroutine multiply_polynomials(p, np, q, nq, result)
        integer, intent(in) :: np, nq
        real(dp), intent(in) :: p(0:np), q(0:nq)
        real(dp), allocatable, intent(out) :: result(:)
        integer :: i, j

        allocate(result(0:np+nq))
        result = 0.0_dp

        do i = 0, np
            do j = 0, nq
                result(i+j) = result(i+j) + p(i)*q(j)
            end do
        end do
    end subroutine multiply_polynomials

    ! ------------------------------------------------------------------------
    ! Fatoracao quadratica do Bairstow.
    ! Para um polinomio de grau n, procura-se um fator
    !       x^2 - r*x - s
    ! que torne o resto aproximadamente nulo.
    !
    ! As relacoes de recorrencia usadas sao:
    !   b(n)   = a(n)
    !   b(n-1) = a(n-1) + r*b(n)
    !   b(i)   = a(i) + r*b(i+1) + s*b(i+2)
    !
    ! Os coeficientes c sao usados para montar o sistema linear de Newton:
    !   c1*dr + c2*ds = -b0
    !   c2*dr + c3*ds = -b1
    ! ------------------------------------------------------------------------
    subroutine bairstow_factor(a, n, r, s, max_used, success, remainder)
        integer, intent(in) :: n
        real(dp), intent(in) :: a(0:n)
        real(dp), intent(inout) :: r, s
        integer, intent(out) :: max_used
        logical, intent(out) :: success
        real(dp), intent(out) :: remainder

        real(dp) :: b(0:n), c(0:n)
        real(dp) :: dr, ds, det
        real(dp) :: scale, err
        integer :: i, it

        success = .false.
        max_used = MAX_ITER
        remainder = huge(1.0_dp)

        do it = 1, MAX_ITER

            ! Calcula os coeficientes b da divisao sintetica quadratica.
            b(n) = a(n)
            b(n-1) = a(n-1) + r*b(n)

            do i = n-2, 0, -1
                b(i) = a(i) + r*b(i+1) + s*b(i+2)
            end do

            remainder = max(abs(b(0)), abs(b(1)))

            ! O vetor c fornece as derivadas necessarias ao passo de Newton.
            c(n) = b(n)
            c(n-1) = b(n-1) + r*c(n)

            do i = n-2, 1, -1
                c(i) = b(i) + r*c(i+1) + s*c(i+2)
            end do

            ! Verificacao relativa do resto.
            scale = 1.0_dp + maxval(abs(a))
            err = remainder/scale

            if (err < TOL) then
                success = .true.
                max_used = it
                return
            end if

            ! Determinante do sistema 2x2 para dr e ds.
            det = c(1)*c(3) - c(2)*c(2)

            ! Quando o sistema linear fica quase singular, esta tentativa
            ! nao e confiavel e deve ser considerada como falha.
            if (abs(det) < DET_TOL) then
                return
            end if

            ! Regra de Cramer aplicada ao sistema de Newton.
            dr = (-b(0)*c(3) + b(1)*c(2))/det
            ds = (-c(1)*b(1) + c(2)*b(0))/det

            ! Atualiza os parametros do fator quadratico.
            r = r + dr
            s = s + ds

            if (.not. (abs(r) < 1.0d12 .and. abs(s) < 1.0d12)) then
                return
            end if
        end do

        remainder = max(abs(b(0)), abs(b(1)))
    end subroutine bairstow_factor

    ! ------------------------------------------------------------------------
    ! Tenta diferentes sementes para aumentar a robustez do metodo.
    ! A rotina e importante porque o proprio enunciado alerta para a
    ! dependencia da convergencia em relacao a r0 e s0.
    ! ------------------------------------------------------------------------
    subroutine bairstow_solve(a, n, roots, nroots, total_iterations, success)
        integer, intent(in) :: n
        real(dp), intent(in) :: a(0:n)
        complex(dp), allocatable, intent(out) :: roots(:)
        integer, intent(out) :: nroots, total_iterations
        logical, intent(out) :: success

        real(dp), allocatable :: work(:), b(:)
        real(dp), parameter :: seeds(2,8) = reshape([ &
            0.0_dp,  0.0_dp, &
            1.0_dp,  0.0_dp, &
           -1.0_dp,  0.0_dp, &
            0.0_dp,  1.0_dp, &
            1.0_dp,  1.0_dp, &
           -1.0_dp,  1.0_dp, &
            2.0_dp, -1.0_dp, &
           -2.0_dp,  1.0_dp], [2,8])

        real(dp) :: r, s, disc
        complex(dp) :: sqrt_disc
        integer :: degree, i, attempt, it, root_pos
        logical :: ok
        real(dp) :: rem

        allocate(work(0:n))
        allocate(roots(1:n))
        work = a
        roots = cmplx(0.0_dp, 0.0_dp, kind=dp)

        degree = n
        root_pos = 0
        total_iterations = 0
        success = .true.

        ! Continua deflacionando o polinomio em fatores quadraticos.
        do while (degree > 2)
            ok = .false.

            do attempt = 1, size(seeds,2)
                r = seeds(1,attempt)
                s = seeds(2,attempt)

                call bairstow_factor(work, degree, r, s, it, ok, rem)

                if (ok) exit
            end do

            if (.not. ok) then
                success = .false.
                return
            end if

            total_iterations = total_iterations + it

            ! Recalcula b para obter explicitamente o polinomio deflacionado.
            allocate(b(0:degree))
            b(degree) = work(degree)
            b(degree-1) = work(degree-1) + r*b(degree)

            do i = degree-2, 0, -1
                b(i) = work(i) + r*b(i+1) + s*b(i+2)
            end do

            ! O novo polinomio e formado pelos coeficientes b(2)...b(n).
            work(0:degree-2) = b(2:degree)
            deallocate(b)

            disc = r*r + 4.0_dp*s
            sqrt_disc = sqrt(cmplx(disc, 0.0_dp, kind=dp))

            root_pos = root_pos + 1
            roots(root_pos) = (r + sqrt_disc)/2.0_dp
            root_pos = root_pos + 1
            roots(root_pos) = (r - sqrt_disc)/2.0_dp

            degree = degree - 2
        end do

        ! Caso final de segundo grau.
        if (degree == 2) then
            disc = work(1)*work(1) - 4.0_dp*work(2)*work(0)
            sqrt_disc = sqrt(cmplx(disc, 0.0_dp, kind=dp))

            root_pos = root_pos + 1
            roots(root_pos) = (-work(1) + sqrt_disc)/(2.0_dp*work(2))
            root_pos = root_pos + 1
            roots(root_pos) = (-work(1) - sqrt_disc)/(2.0_dp*work(2))

        ! Caso final de primeiro grau.
        else if (degree == 1) then
            root_pos = root_pos + 1
            roots(root_pos) = cmplx(-work(0)/work(1), 0.0_dp, kind=dp)
        end if

        nroots = root_pos
    end subroutine bairstow_solve

    ! ------------------------------------------------------------------------
    ! Resolve somente a primeira fatoracao quadratica, permitindo estudar
    ! diretamente a influencia de diferentes pares (r0,s0).
    ! ------------------------------------------------------------------------
    subroutine analyze_initial_pair(a, n, r0, s0, iterations, converged, rem)
        integer, intent(in) :: n
        real(dp), intent(in) :: a(0:n)
        real(dp), intent(in) :: r0, s0
        integer, intent(out) :: iterations
        logical, intent(out) :: converged
        real(dp), intent(out) :: rem

        real(dp) :: r, s

        r = r0
        s = s0
        call bairstow_factor(a, n, r, s, iterations, converged, rem)
    end subroutine analyze_initial_pair

    ! ------------------------------------------------------------------------
    ! Avalia o polinomio em x usando Horner.
    ! ------------------------------------------------------------------------
    pure function evaluate_polynomial(a, n, x) result(value)
        integer, intent(in) :: n
        real(dp), intent(in) :: a(0:n)
        complex(dp), intent(in) :: x
        complex(dp) :: value
        integer :: i

        value = cmplx(a(n), 0.0_dp, kind=dp)
        do i = n-1, 0, -1
            value = value*x + cmplx(a(i), 0.0_dp, kind=dp)
        end do
    end function evaluate_polynomial

end module bairstow_module


program ppc2_bairstow
    use bairstow_module
    implicit none

    integer, parameter :: n_test = 7
    integer, parameter :: n_apc = 4
    integer, parameter :: ngrid = 201

    real(dp) :: p_test(0:n_test)
    real(dp) :: p_apc(0:n_apc)
    complex(dp), allocatable :: roots(:)

    integer :: nroots, iterations, i, j
    integer :: unit_out, unit_conv, unit_fractal
    logical :: success, converged
    real(dp) :: r0, s0, rem
    real(dp) :: rmin, rmax, smin, smax
    integer :: conv_count

    ! ------------------------------------------------------------------------
    ! POLINOMIO DE VALIDACAO DE 7o GRAU
    !
    ! Raizes conhecidas:
    !   1, -2, 3,
    !   0.4 +/- 1.2 i,
    !  -1.5 +/- 0.7 i
    !
    ! Fatores usados:
    !   (x-1)(x+2)(x-3)
    !   (x^2 - 0.8x + 1.6)
    !   (x^2 + 3x + 2.74)
    ! ------------------------------------------------------------------------

    ! Para manter o programa simples e totalmente autonomo, o polinomio de
    ! teste e inserido diretamente na forma expandida abaixo. Esses coeficientes
    ! resultam da multiplicacao dos fatores comentados acima.
    p_test = [ &
        26.304_dp, -6.272_dp, -10.168_dp, 2.668_dp, &
        -6.272_dp, -7.460_dp, 0.2_dp, 1.0_dp ]

    ! ------------------------------------------------------------------------
    ! POLINOMIO CARACTERISTICO DA APC2
    ! P(lambda) = 3 lambda^4 + 2.2 lambda^3 + 54.36 lambda^2
    !             + 17.6 lambda + 215
    ! ------------------------------------------------------------------------

    p_apc = [215.0_dp, 17.6_dp, 54.36_dp, 2.2_dp, 3.0_dp]

    ! ------------------------------------------------------------------------
    ! 1) VALIDACAO DO METODO
    ! ------------------------------------------------------------------------

    call bairstow_solve(p_test, n_test, roots, nroots, iterations, success)

    open(newunit=unit_out, file='validacao.txt', status='replace', action='write')

    write(unit_out,'(A)') 'VALIDACAO DO METODO DE BAIRSTOW - POLINOMIO DE 7o GRAU'
    write(unit_out,'(A)') '-----------------------------------------------------------'
    write(unit_out,'(A)') 'Raizes teoricas:'
    write(unit_out,'(A)') '1.0'
    write(unit_out,'(A)') '-2.0'
    write(unit_out,'(A)') '3.0'
    write(unit_out,'(A)') '0.4 +/- 1.2 i'
    write(unit_out,'(A)') '-1.5 +/- 0.7 i'
    write(unit_out,'(A)') ''
    write(unit_out,'(A,L1)') 'Convergencia: ', success
    write(unit_out,'(A,I0)') 'Iteracoes totais: ', iterations
    write(unit_out,'(A)') 'Raizes calculadas:'

    do i = 1, nroots
        write(unit_out,'(F18.10,1X,F18.10,A)') real(roots(i),dp), aimag(roots(i)), ' i'
    end do

    write(unit_out,'(A)') ''
    write(unit_out,'(A)') 'Observacao: as raizes complexas devem aparecer em pares conjugados.'

    close(unit_out)

    ! ------------------------------------------------------------------------
    ! 2) ANALISE DE DIFERENTES VALORES INICIAIS r0 e s0
    ! ------------------------------------------------------------------------

    open(newunit=unit_conv, file='convergencia.txt', status='replace', action='write')

    write(unit_conv,'(A)') 'ANALISE DA CONVERGENCIA - POLINOMIO DA APC2'
    write(unit_conv,'(A)') 'r0              s0              iteracoes       convergiu       resto'

    block
        real(dp), parameter :: r0_list(8) = [0.0_dp, 1.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, -1.0_dp, 2.0_dp, -2.0_dp]
        real(dp), parameter :: s0_list(8) = [0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, -1.0_dp, 1.0_dp]

        do j = 1, 8
            r0 = r0_list(j)
            s0 = s0_list(j)

        call analyze_initial_pair(p_apc, n_apc, r0, s0, iterations, converged, rem)

            write(unit_conv,'(F12.4,4X,F12.4,8X,I6,10X,L5,6X,ES14.6)') &
                r0, s0, iterations, converged, rem
        end do
    end block

    close(unit_conv)

    ! ------------------------------------------------------------------------
    ! 3) RESOLUCAO DO POLINOMIO DA APC2
    ! ------------------------------------------------------------------------

    call bairstow_solve(p_apc, n_apc, roots, nroots, iterations, success)

    open(newunit=unit_out, file='raizes_apc2.txt', status='replace', action='write')

    write(unit_out,'(A)') 'RAIZES DO POLINOMIO CARACTERISTICO DA APC2'
    write(unit_out,'(A)') 'P(lambda) = 3 lambda^4 + 2.2 lambda^3 + 54.36 lambda^2 + 17.6 lambda + 215'
    write(unit_out,'(A,L1)') 'Convergencia: ', success
    write(unit_out,'(A,I0)') 'Iteracoes totais: ', iterations
    write(unit_out,'(A)') ''
    write(unit_out,'(A)') 'Numero      Parte real              Parte imaginaria'

    do i = 1, nroots
        write(unit_out,'(I4,8X,F18.10,8X,F18.10)') i, real(roots(i),dp), aimag(roots(i))
    end do

    write(unit_out,'(A)') ''
    write(unit_out,'(A)') 'Interpretacao:'
    write(unit_out,'(A)') 'lambda = alpha + i beta'
    write(unit_out,'(A)') 'alpha < 0: decaimento da resposta.'
    write(unit_out,'(A)') 'beta != 0: comportamento oscilatorio.'

    close(unit_out)

    ! ------------------------------------------------------------------------
    ! 4) FRACTAL DE BAIRSTOW
    ! Varredura de r0 e s0 no intervalo proposto no enunciado: [-3,3].
    ! Para cada ponto, registra-se o numero de iteracoes da primeira
    ! fatoracao quadratica. Valor -1 indica divergencia/falha.
    ! ------------------------------------------------------------------------

    rmin = -3.0_dp
    rmax =  3.0_dp
    smin = -3.0_dp
    smax =  3.0_dp

    open(newunit=unit_fractal, file='fractal.dat', status='replace', action='write')

    write(unit_fractal,'(A)') '# r0 s0 iteracoes'

    conv_count = 0

    do j = 0, ngrid-1
        s0 = smin + (smax-smin)*real(j,dp)/real(ngrid-1,dp)

        do i = 0, ngrid-1
            r0 = rmin + (rmax-rmin)*real(i,dp)/real(ngrid-1,dp)

            call analyze_initial_pair(p_apc, n_apc, r0, s0, iterations, converged, rem)

            if (converged) then
                write(unit_fractal,'(F12.6,1X,F12.6,1X,I6)') r0, s0, iterations
                conv_count = conv_count + 1
            else
                write(unit_fractal,'(F12.6,1X,F12.6,1X,I6)') r0, s0, -1
            end if
        end do
    end do

    close(unit_fractal)

    open(newunit=unit_out, file='resumo.txt', status='replace', action='write')
    write(unit_out,'(A)') 'RESUMO DA EXECUCAO DO PPC2'
    write(unit_out,'(A)') '============================='
    write(unit_out,'(A)') '1. Validacao: validacao.txt'
    write(unit_out,'(A)') '2. Convergencia: convergencia.txt'
    write(unit_out,'(A)') '3. Raizes da APC2: raizes_apc2.txt'
    write(unit_out,'(A)') '4. Dados do fractal: fractal.dat'
    write(unit_out,'(A,I0)') '5. Pontos convergentes no fractal: ', conv_count
    write(unit_out,'(A,I0)') '6. Total de pontos analisados: ', ngrid*ngrid
    close(unit_out)

    write(*,'(A)') 'PPC2 executado com sucesso.'
    write(*,'(A)') 'Arquivos gerados:'
    write(*,'(A)') '  validacao.txt'
    write(*,'(A)') '  convergencia.txt'
    write(*,'(A)') '  raizes_apc2.txt'
    write(*,'(A)') '  fractal.dat'
    write(*,'(A)') '  resumo.txt'

end program ppc2_bairstow
```

## 12. Compilação e execução

Com o arquivo salvo como `main.f90`, utilizando o `gfortran`:

```bash
gfortran -O2 -Wall -Wextra -std=f2008 main.f90 -o ppc2_bairstow
./ppc2_bairstow
```

A execução gera os arquivos de resultados utilizados na análise:

```text
validacao.txt
convergencia.txt
raizes_apc2.txt
fractal.dat
resumo.txt
```

---

## 13. Observação sobre a implementação

O algoritmo foi implementado diretamente em Fortran, utilizando vetores, funções, subrotinas, laços e estruturas condicionais. Não foi utilizada uma biblioteca externa para calcular diretamente as raízes do polinômio.

O código implementa a divisão sintética quadrática, o ajuste de `r` e `s` pelo sistema linear de Newton, a obtenção das raízes pela fórmula quadrática e a deflação sucessiva do polinômio.

A rotina também testa diferentes sementes iniciais para melhorar a robustez da fatoração e realiza a varredura sistemática de `r0` e `s0` utilizada na construção do mapa de Bairstow.
