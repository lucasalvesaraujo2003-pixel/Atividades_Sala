! ==============================================================================
! PROGRAMA: atividadeNR (Resolução de Equações Não Lineares)
!
! 1. RESUMO DOS MÉTODOS NUMÉRICOS:
!    O programa implementa algoritmos baseados na família do Método de Newton-Raphson:
!    - Método de Newton-Raphson Padrão (comentado): Utiliza a primeira derivada
!      da função para traçar uma reta tangente ao ponto atual. A interseção dessa
!      reta com o eixo x gera a nova aproximação.
!    - Método de Newton-Raphson Modificado (Múltiplas Raízes): Variante ativada 
!      no loop principal. Ela incorpora a derivada segunda da função para manter
!      a taxa de convergência quadrática e evitar instabilidades quando a equação 
!      possui raízes com multiplicidade maior que um (raízes duplas, triplas, etc.).
!
! 2. PARA QUE SERVE (APLICAÇÃO GERAL DOS MÉTODOS):
!    O Método de Newton-Raphson e as suas variantes servem para aproximar raízes 
!    de equações algébricas e transcendentais de forma extremamente rápida. Por 
!    apresentar convergência quadrática (o número de dígitos corretos quase dobra 
!    a cada iteração), é amplamente empregado na engenharia e computação científica 
!    para resolver sistemas de equações não lineares, problemas de otimização, 
!    simulações de fluxo de fluidos, análise estrutural não linear, cinemática de 
!    mecanismos e calibração instantânea de modelos matemáticos complexos.
!
! 3. INPUTS (ENTRADAS DEFINIDAS NO CÓDIGO):
!    - x1   : Chute/estimativa inicial para o início da busca (Valor: 3.1)
!    - f(x) : A função polinomial avaliada: x^3 - 6x^2 + 9x - 4
!    - f_d1 : Primeira derivada analítica da função: 3x^2 - 12x + 9
!    - f_d2 : Segunda derivada analítica da função: 6x - 12
!
! 4. OUTPUTS (SAÍDAS GERADAS):
!    Os dados são exibidos diretamente no terminal (output padrão) a cada iteração:
!    - Mensagens em tela contendo as seguintes colunas ordenadas:
!      [Número da Iteração (iter)]  [x_atual (x1)]  [x_próximo (x2)]  [f(x_próximo)]
! ==============================================================================

program atividadeNR
implicit none

real x1		! xi
real x2		! xi+1
real fx		! f(x)
integer iter	! iterações

x1 = 3.1	! Chute inicial
fx = f(x1)	! f(x) para o chute inicial
iter = 1.0

! Método de NR normal
 do while (fx.ne.0)
 x2 = x1 - f(x1)/f_d1(x1)
 fx = f(x2)
 write(*,*) iter, x1, x2, fx
 iter = iter + 1
 x1 = x2
 end do

! Método de NR para múltiplas raízes
!do while (fx.ne.0)
!x2 = x1 - (f(x1)*f_d1(x1))/(f_d1(x1)**2 - f(x1)*f_d2(x1))
!fx = f(x2)
!write(*,*) iter, x1, x2, fx
!iter = iter + 1
!x1 = x2
!end do

contains
real function f(x)	! Função f(x)
implicit none
real, intent(in) :: x
f = x**3.0 - 6.0*x**2.0 + 9.0*x - 4.0
end function f

real function f_d1(x)	! Derivada primeira de f(x)
implicit none
real, intent(in) :: x
f_d1 = 3.0*x**2.0 - 12.0*x + 9.0
end function f_d1

real function f_d2(x)	! Derivada segunda de f(x)
implicit none
real, intent(in) :: x
f_d2 = 6.0*x - 12.0
end function f_d2

end program atividadeNR
