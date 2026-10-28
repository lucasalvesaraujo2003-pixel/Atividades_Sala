! ==============================================================================
! PROGRAMA: bisseccao (Resolução de Equações Não Lineares)
!
! 1. RESUMO DOS MÉTODOS NUMÉRICOS:
!    O programa resolve uma equação não linear f(x) = 0 comparando dois métodos:
!    - Método da Bissecção: Divide o intervalo de busca puramente ao meio a cada 
!      iteração (xm = (xi + xs) / 2). É extremamente robusto e sempre converge,
!      desde que haja uma troca de sinal no intervalo inicial.
!    - Método da Falsa Posição (Regula Falsi): Conecta os pontos f(xi) e f(xs)
!      por uma linha reta. A interseção dessa reta com o eixo x define a nova
!      aproximação (xr). Geralmente converge mais rápido que a bissecção.
!
! 2. PARA QUE SERVE (APLICAÇÃO GERAL DOS MÉTODOS):
!    Estes métodos servem para encontrar raízes (ou zeros) de funções algébricas 
!    ou transcendentais complexas onde não é possível isolar a variável 
!    independente por métodos analíticos tradicionais (álgebra comum). São 
!    amplamente utilizados na engenharia e ciências exatas para resolver 
!    problemas de otimização, determinação de pontos de equilíbrio, cálculo de 
!    fatores de atrito, dimensionamento estrutural e calibração de modelos 
!    físicos complexos através de aproximações sucessivas controladas.
!
! 3. INPUTS (ENTRADAS DEFINIDAS NO CÓDIGO):
!    - xi     : Limite inferior do intervalo inicial de busca (Valor: 12.0)
!    - xs     : Limite superior do intervalo inicial de busca (Valor: 16.0)
!    - tol    : Tolerância aceitável para o critério de parada (Valor: 1.0e-05)
!    - raizv  : Valor analítico/verdadeiro da raiz para calibração (14.7802...)
!
! 4. OUTPUTS (SAÍDAS GERADAS):
!    O programa gera dois arquivos de texto com os resultados das iterações:
!    - 'saida.dat'  : Resultados detalhados do Método da Bissecção.
!    - 'saida2.dat' : Resultados detalhados do Método da Falsa Posição.
!
!    Estrutura das colunas em ambos os arquivos de saída:
!    [Iteração]  [xi]  [xs]  [Raiz Aproximada (xm ou xr)]  [Erro Est.]  [Erro Verd.]
! ==============================================================================

program bisseccao

! Desejamos escrever um programa que resolva a raiz da seguinte equação:

! f(x) = (667.38/x)*(1.0 - exp(-0.146843*x)) - 40.0

! O método a ser utilizado deve ser o método da bissecção e pela falsa posição

! Primeiro, vamos declarar nossas variáveis

real xi, xs     ! limite inferior e superior do intervalo de busca
real xm, fm     ! valor parcial da raiz e da função f(xm) - Método da Bissecção
real fi, fs     ! valor da função f(x) em xi e xs respectivamente
real xr, fr     ! nova aproximação da raiz e da função - Método da Falsa Posição
real tol        ! tolerância para avaliação de convergência
real erro       ! erro relativo para avaliação da convergência
real xmd	! valor provisório da raiz depois de um passo
real xrd
integer iter	! número da iteração

! Definindo o valor da raiz verdadeira
raizv = 14.780208593679468

! Método da bissecção:
! Definindo o intervalo inicial de busca
xi=12.0
xs=16.0

! Definindo a tolerância para fins de convergência
tol = 1.0e-05
erro =1.0
iter =0

open(1, file='saida.dat')
do while(erro.gt.tol)

	! Avaliando os valores da função nos intervalos xi e xs
	fi=funcao(xi)
	fs=funcao(xs)

	! Averiguando se a função muda de sinal no intervalo
	if ((fi*fs).lt.0.0) then
		xm= (xi+xs)/2.0 ! valor parcial da raiz

		! Checando agora se "xi < raiz < xm" ou "xm < raiz < xs"
		fm=funcao(xm)

		if((fi*fm).lt.0.0) then
			xs=xm
		else
			xi=xm
		end if
		xmd = (xi+xs)/2.0
		erro= abs(xm-xmd)	! Erro relativo
		errov = abs(xmd-raizv)	! Erro verdadeiro
		iter =iter + 1
		write(1,*) iter, xi, xs, xm, erro, errov
	
	end if
end do
close(1)

! Método da falsa posição:
! Definindo o intervalo inicial de busca
xi=12.0
xs=16.0

! Definindo a tolerância para fins de convergência
tol = 1.0e-05
erro =1.0
errov = 0.0
iter =0

open(2, file='saida2.dat')
do while (erro.gt.tol)
	! Avaliando os valores da função nos intervalos xi e xs
	fi=funcao(xi)
	fs=funcao(xs)
	
	! Averiguando se a função muda de sinal no intervalo
	if ((fi*fs).lt.0.0) then
		xr = xs - (fs*(xi-xs))/(fi-fs)
		fr=funcao(xr)
		
		if((fi*fr).lt.0.0) then
			xs=xr
		else
			xi=xr
		end if
		fi = funcao(xi)
		fs = funcao(xs)
		xrd = xs - (fs*(xi-xs))/(fi-fs)
		erro= abs(xr-xrd)	! Erro relativo
		errov = abs(xrd-raizv)	! Erro verdadeiro
		iter =iter + 1
		write(2,*) iter, xi, xs, xr, erro, errov
	end if
end do

! Defindo a função f(x)
contains
function funcao(x) result(f)
	real, intent(in) :: x
	real f
	f = (667.38/x)*(1.0 - exp(-0.146843*x)) - 40.0
end function funcao

end program bisseccao
