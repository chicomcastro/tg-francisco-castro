# Proposta de TG

# Informações primárias

I. nome(s) do(s) Orientador(es): Maisa de Oliveira Terra

Aluno: Francisco Matheus Moreira de Castro

II. título do TG:

Estudo de Otimização de Parâmetros em Transferências Terra-Marte com Swing-by Propulsado por Vênus 

## Proposta

Este Trabalho de Graduação objetiva estudar a Otimização de Parâmetros em Projetos de Trajetórias de Transferência Terra-Marte com Swing-by Propulsado por Vênus, buscando minimizar o custo total e/ou tempo de voo a partir da chamada abordagem de Problemas de Três Corpos Colados (Patched-3B approach).

Castro et al. (2019) estudaram o projeto de transferências Terra-Marte a partir da modelagem de cônicas acopladas e aplicação de algoritmo evolutivo para obtenção de soluções ótimas, visando minimizar custo total da missão.

Recentemente, Terra e Prado (2019) investigaram o papel do campo gravitacional de Júpiter em Transferências Terra-Marte (EMT) com a inclusão de uma manobra de Swing-by por Vênus, explorando soluções do Problema Restrito de Três Corpos. Esta análise mostrou que, em geral, após a passagem por Vênus a espaçonave se aproxima da  esfera de influência de Marte com uma velocidade muito mais elevada que a velocidade orbital de Marte, implicando num custo maior desta missão em relação à uma EMT direta. Porém, mostra-se que uma economia de custo de até 6% é possível através da escolha adequada da fase entre Sol e Júpiter no instante inicial da trajetória de transferência. 

Por outro lado, Ferreira (2017) explorou como uma manobra de Swing-by com aplicação de  impulso de módulo e direção variáveis em um ponto arbitrário na passagem pelo primário pode otimizar o ganho ou perda de energia do veículo espacial, mostrando inclusive que resultados melhores do que aqueles obtidos com o impulso aplicado no periapsis podem ser obtidos.

Desta forma, pretende-se neste Trabalho estudar o projeto de uma EMT incluindo a aplicação de um impulso instantâneo de módulo e direção variáveis ao longo da passagem por Vênus, buscando definir e resolver um problema de otimização.  Com esta abordagem pretende-se obter soluções ótimas de baixo custo em projetos preliminares com dupla missão, i.e., passagem por Vênus e Marte. Como parte final, pode-se incluir a análise para o retorno a Terra.

## Plano de Trabalho

Durante o TG-1, as atividades a serem realizadas envolvem:

1. Levantamento bibliográfico dos aspectos relevantes do problema:
(i) Modelos adotados em projeto de EMT da literatura e resultados obtidos (incluindo Terra e Prado (2019) e Conte et al. (2018)).
(ii) Modelagem da manobra de Swing-by propulsada, incluindo Ferreira, (2017, 2015) e referências contidas.
(iii) Estratégias de otimização multiparamétrica, tais como a aplicada em Castro (2019).

2. Estudo do modelo matemático do Swing-by propulsado a partir do Modelo Restrito de Três Corpos e definição de parâmetros associados relevantes para o problema de otimização.

3. Entender modelagem de EMT com Swing-by por Vênus clássico (sem propulsão) baseado na dinâmica do Sistema Sol-Júpiter para posterior extensão ao caso propulsado.

4. Elaboração de algoritmo de solução considerando a inclusão de propulsão em pontos variáveis,
com levantamento de parâmetros a explorar para possíveis abordagens de cálculo: (a) Cola de Sistemas de Três Corpos, (b) Modelagem direta de N-corpos.

5. Implementação do código com a inclusão de otimização aplicada em Castro et al. (2019) para solução do EMT com Swing-by por Vênus sem aplicação de propulsão. Validação de resultados com comparação com resultados de  Terra e Prado (2019).

Durante o TG-2, as atividades a serem realizadas envolvem:

1. Inclusão de aplicação de propulsão e busca de soluções ótimas. Análise de dependência com parâmetros relevantes.

2. Análise de soluções de retorno a Terra.

3. Com base nas análises realizadas, elaborar propostas de trabalhos futuros e aprimoramentos.

## Referências

CASTRO, F.M.M.; MARINOT, G; MEDEIROS, G.; MENDES, R.; SIVIERI, L.E., Projeto de transferência Terra-Marte baseado em patched conics otimizado com algoritmo evolutivo, Trabalho de Final de Disciplina MVO-41 (Mecânica Orbital), ITA, 2019.

CONTE et al., Earth-Mars transfers through Moon Distant Retrograde Orbits, Acta , v. 143, p. 372–379, 2018.

FERREIRA, A.F.S., Manobras orbitais combinadas com uso de propulsão impulsiva
e passagem por um corpo intermediário, Tese de Doutorado, INPE, São José dos Campos, 2017.

FERREIRA, A., PRADO, A.; WINTER, O., A numerical study of powered swing-bys around the moon, Advances in Space Research, v. 56, p. 252–272, 2015.

TERRA, M.O., PRADO, A.B.A., Contributions of Venus Swing-by maneuver in Earth-Mars transfers, 70th International Astronautical Congress, Washington D.C., Paper number IAC-19-C1.4.6, 2019.
