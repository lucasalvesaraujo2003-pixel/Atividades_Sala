
# Cálculo Numérico Aplicado (CNA) — Engenharia Mecânica

- **Instituição:** Universidade de Brasília (UnB)
- **Faculdade:** Faculdade de Tecnologia (FT)
- **Departamento:** Departamento de Engenharia Mecânica (ENM)
- **Aluno:** Lucas Alves Araujo
- **Matrícula:** 221021705
- **Semestre letivo:** 2026/2

## Sobre o repositório

Este repositório tem como objetivo reunir, organizar e documentar os algoritmos, códigos-fonte, relatórios e atividades desenvolvidos na disciplina de **Cálculo Numérico Aplicado (CNA)**, do curso de Engenharia Mecânica da Universidade de Brasília (UnB), durante o semestre letivo de 2026/1.

Os materiais estão organizados por práticas de programação computacional (PPC), atividades práticas complementares (APC) e trabalhos desenvolvidos em sala de aula. O projeto busca facilitar a consulta aos conteúdos, a revisão dos métodos numéricos e a reprodução dos resultados obtidos.

## Motivação para o estudo dos métodos numéricos

Na Engenharia Mecânica, diversos problemas envolvem modelos matemáticos complexos cuja solução nem sempre pode ser obtida por métodos analíticos. Nesse contexto, os métodos numéricos permitem encontrar soluções aproximadas para problemas de engenharia, sendo aplicados em áreas como transferência de calor, mecânica dos fluidos, análise estrutural, termodinâmica e otimização.

O estudo desses métodos possibilita compreender como os problemas matemáticos podem ser transformados em algoritmos computacionais. Mais do que utilizar ferramentas prontas, é importante conhecer os fundamentos matemáticos, as etapas de implementação, os critérios de convergência e as possíveis fontes de erro envolvidas nos cálculos.

A programação em Fortran é utilizada neste repositório para implementar e analisar diferentes algoritmos numéricos, permitindo relacionar a teoria estudada em sala de aula com aplicações práticas da Engenharia Mecânica. Também são utilizados recursos complementares, como Python e Gnuplot, para o processamento de resultados e a construção de gráficos.

Outro objetivo é adotar boas práticas de organização, documentação e validação dos códigos. A separação dos arquivos por atividade, a descrição das variáveis e o registro dos resultados facilitam a compreensão dos programas, a identificação de erros e a reprodução dos cálculos.

Dessa forma, este repositório funciona como um material de apoio ao aprendizado, reunindo implementações computacionais e registros das atividades desenvolvidas ao longo da disciplina.

## Estrutura do repositório

O projeto está organizado em diretórios correspondentes às práticas computacionais e aos trabalhos desenvolvidos em sala de aula.

```text
cna-unb/
├── README.md
├── livro-cna.pdf
│
├── PPC_01/
│   ├── README.md
│   ├── PPC1 e APC1.pdf
│   ├── APC_1_[identificacao].pdf
│   ├── ppc1-rk4.f90
│   ├── plot1.gnu
│   └── outputs/
│
├── PPC_02/
│   ├── README.md
│   ├── PPC2 e APC2.pdf
│   ├── APC_2_[identificacao].pdf
│   ├── ppc2-bairstow.f90
│   ├── varredura.f90
│   ├── gerar_fractal.f90
│   └── outputs/
│       ├── dados_varredura.csv
│       └── fractal_bairstow.png
│
├── PPC_03/
│   ├── README.md
│   ├── PPC3 e APC3.pdf
│   ├── APC_3_[identificacao].pdf
│   ├── ppc3-thomas.f90
│   ├── campo_com_geracao.gnu
│   ├── evolucao_com_geracao.gnu
│   ├── campo_sem_geracao.gnu
│   ├── evolucao_sem_geracao.gnu
│   └── outputs/
│
├── Atividades_Sala

*Observação: a árvore acima é uma representação resumida da organização proposta. Os nomes dos arquivos e diretórios devem corresponder aos que realmente existem no repositório.*

## Tecnologias utilizadas

- **Fortran:** implementação dos algoritmos numéricos.
- **Python:** processamento de dados e geração de gráficos.
- **Gnuplot:** visualização gráfica dos resultados.
- **Git e GitHub:** versionamento, organização e armazenamento dos arquivos.

## Organização das atividades

Cada diretório de prática contém, conforme as necessidades da atividade:

- **Código-fonte:** implementação do método numérico em Fortran.
- **Documentação:** descrição do problema, fundamentação teórica, variáveis, entradas, saídas e funcionamento do algoritmo.
- **Roteiro da atividade:** enunciado e orientações fornecidas para a prática.
- **Relatório final:** documentação dos procedimentos, resultados e conclusões.
- **Arquivos de saída:** dados numéricos, tabelas, registros e gráficos gerados durante a execução.

O conteúdo de cada diretório pode variar de acordo com os requisitos específicos da atividade.

## Objetivos de aprendizagem

- Compreender os fundamentos matemáticos dos métodos numéricos estudados.
- Implementar algoritmos computacionais aplicados a problemas de engenharia.
- Analisar a precisão, a convergência e as limitações dos métodos utilizados.
- Comparar diferentes abordagens para a resolução de problemas numéricos.
- Desenvolver habilidades de programação científica e interpretação de resultados.
- Adotar práticas de documentação, organização e validação de códigos.

## Considerações finais

Este repositório representa um registro do desenvolvimento acadêmico na disciplina de Cálculo Numérico Aplicado, integrando conceitos matemáticos, programação e aplicações da Engenharia Mecânica.

A organização dos materiais busca tornar o conteúdo mais acessível para consultas futuras, revisões e aprofundamento dos métodos estudados, contribuindo para a formação técnica e computacional em engenharia.
