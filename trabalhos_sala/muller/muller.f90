! ==============================================================================
! PROGRAMA: muller (Resolução de Equações Não Lineares)
!
! 1. RESUMO DO MÉTODO NUMÉRICO:
!    O programa implementa o Método de Müller. Diferente do método da Secante, 
!    que aproxima a curva da função por uma reta ligando dois pontos, o método 
!    de Müller aproxima a função por uma parábola (polinómio de segundo grau) 
!    que passa por três pontos dados (x0, x1 e x2). A interseção dessa parábola 
!    com o eixo x que estiver mais próxima do último ponto estimado torna-se a 
!    nova aproximação (x3). O uso de termos quadráticos e raízes quadradas permite 
!    trabalhar nativamente com números e raízes complexas.
!
! 2. PARA QUE SERVE (APLICAÇÃO GERAL DO MÉTODO):
!    O Método de Müller serve para encontrar tanto raízes reais quanto raízes 
!    complexas de equações algébricas e transcendentais complexas, sem a 
!    necessidade de calcular as derivadas da função. É amplamente utilizado na 
!    engenharia e computação científica para mapear modos de vibração estrutural 
!    com amortecimento, analisar estabilidade em sistemas de controlo, resolver 
!    equações características em mecânica quântica e eletromagnetismo, e encontrar 
!    raízes de polinómios de alto grau onde se espera obter respostas no plano 
!    complexo.
!
! 3. INPUTS (ENTRADAS DEFINIDAS NO CÓDIGO):
!    - x0, x1, x2 : Três estimativas (chutes) iniciais consecutivas calculadas 
!                   a partir de um ponto inicial (-15.0) e um passo (0.5).
!    - passo      : Distância incremental para definir x1 e x2 (Valor: 0.5)
!    - tol        : Tolerância aceitável para o critério de paragem (Valor: 1.0e-04)
!    - f(x)       : A função cuja raiz se deseja encontrar (definida na secção contains).
!
! 4. OUTPUTS (SAÍDAS GERADAS):
!    Os dados finais são exibidos no terminal (output padrão) após convergência:
!    - Exibição de texto contendo os valores calculados:
!      [Raiz Aproximada (x3)] e [Número de Iterações realizadas]
! ==============================================================================

program muller
implicit none

complex x0, x1, x2, x3		! Chutes iniciais
complex fx0, fx1, fx2, fx3		! Funções dos chutes
complex h0, h1
complex delta0, delta1
complex a, b, c
real passo, tol, erro
real bask1, bask2

passo = 0.5
x0 = -15.0
x1 = x0 + passo
x2 = x1 + passo
erro = 1
tol = 1.0e-04

do while (erro.gt.tol)
fx0 = f(x0)
fx1 = f(x1)
fx2 = f(x2)

h0 = x1 - x0
h1 = x2 - x1

delta0 = (fx1 - fx0)/h0
delta1 = (fx2 - fx1)/h1

a = (delta1 - delta0)/(h1 + h0)
b = a*h1 + delta1
c = fx2

bask1 = abs(b + sqrt(b**2 - 4*a*c))
bask2 = abs(b - sqrt(b**2 - 4*a*c))

if (bask1.gt.bask2) then
x3 = x2 - (2*c)/(b + sqrt(b**2 - 4*a*c))
else
x3 = x2 - (2*c)/(b - sqrt(b**2 - 4*a*c))
end if

erro = abs(x3 - x2)/abs(x3)

x2 = x1
x1 = x0
x0 = x3

write(*,*) x3, erro
end do

contains
complex function f(x)	! Função f(x)
implicit none
complex, intent(in) :: x
f = x**3.0 - 13.0*x - 12.0
end function f

end program muller
