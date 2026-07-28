# Isomorphism

Framework em Python para experimentação com algoritmos de isomorfismo de grafos e problemas correlatos (subgrafo comum máximo, distância estrutural, comparação de proteínas), desenvolvido como parte de um Trabalho de Conclusão de Curso (TCC) sobre algoritmos de isomorfismo em grafos com complexidade polinomial sob classes restritas de grafos.

> **Aviso de escopo:** este projeto é experimental e acadêmico. Nenhum dos algoritmos aqui implementados resolve o problema de isomorfismo de grafos em sua forma geral — cada um deles é válido apenas dentro de uma classe restrita de grafos (árvores, grafos outerplanares, grafos moleculares, etc.) ou fornece garantias parciais (heurísticas, invariantes necessários mas não suficientes). Veja a seção [Limitações por algoritmo](#limitações-por-algoritmo) antes de usar os resultados como prova de corretude.

## Estrutura do repositório

```
Isomorphism/
├── main.py                 # CLI para comparar dois grafos usando qualquer algoritmo
├── structures/
│   └── Grafo.py             # Classe Grafo (wrapper sobre networkx.Graph)
├── algorithms/
│   ├── hausdorff.py                      # Distância de Hausdorff entre árvores (Kelenc, 2021)
│   ├── mcs_bbp_algorithm.py              # MCS via block-and-bridge preserving (Schietgat et al., 2013)
│   ├── mol_graph_iso.py                  # Isomorfismo molecular (Faulon, 1998)
│   ├── polynomial_algorithm_undirected.py # Isomorfismo via invariantes espectrais (He et al., 2019)
│   └── protein_iso.py                    # Distância quantitativa para proteínas (Guo et al., 2022)
└── tests/
    ├── hausdorff_tests.py
    ├── mcs_bbp_algorithm_tests.py
    ├── mol_graph_iso_tests.py
    ├── polynimial_algorithm_undirected_tests.py
    └── protein_iso_tests.py
```

## Instalação

Requer **Python 3.13** (o repositório foi desenvolvido e testado nessa versão; versões 3.10+ provavelmente funcionam, mas não foram testadas).

```bash
git clone https://github.com/Leo20031999/Isomorphism.git
cd Isomorphism
pip install -r requirements.txt
```

### Dependências e versões

O repositório não possui `requirements.txt` no momento — crie um na raiz com o seguinte conteúdo (ajuste as versões para as que você efetivamente usou, via `pip freeze`):

```
networkx>=3.0
numpy>=1.24
scipy>=1.10
psutil>=5.9
```

| Pacote | Uso no projeto |
|---|---|
| `networkx` | Estrutura de grafo subjacente (`structures/Grafo.py`), leitura de `.graphml`/`.gml`, checagem de outerplanaridade |
| `numpy` | Matrizes de adjacência, cálculo de autovalores/potências (He et al., Guo et al.) |
| `scipy` | Atribuição ótima (`linear_sum_assignment`, Hausdorff), SVD e matrizes esparsas (Guo et al.) |
| `psutil` | Apenas nos scripts de teste, para monitorar CPU/memória durante os testes de estresse |

## Uso

### Via linha de comando

```bash
python main.py grafo1.txt grafo2.txt --algoritmo isomorfismo
python main.py proteina1.graphml proteina2.graphml --algoritmo todos
python main.py molecula1.gml molecula2.gml --algoritmo mcs --verbose
```

**Algoritmos disponíveis** (`--algoritmo`):

| Valor | Algoritmo | Referência |
|---|---|---|
| `isomorfismo` | Isomorfismo em grafos não direcionados via invariantes | He et al., 2019 |
| `mcs` | Subgrafo comum máximo (grafos outerplanares) | Schietgat et al., 2013 |
| `isomorfismo_molecular` | Isomorfismo/automorfismo para grafos moleculares | Faulon, 1998 |
| `distancia_proteina` | Distância quantitativa estrutural | Guo et al., 2022 |
| `hausdorff` | Distância de Hausdorff entre árvores | Kelenc, 2021 |
| `automorfismos` | Automorfismos do primeiro grafo | Faulon, 1998 |
| `todos` (padrão) | Executa todos os algoritmos aplicáveis | — |

Outras opções: `--verbose` (exibe estatísticas do grafo e detalhes da execução) e `--formato {auto,graphml,gml,txt}` (padrão detecta pela extensão do arquivo).

### Formato de arquivo `.txt` (formato próprio)

```
vertices:
1: C
2: C
3: O

arestas:
1 2
2 3

rotulos:
1: Carbono

atributos:
1 atom_type:C
```

Seções `rotulos:` e `atributos:` são opcionais. Arquivos `.graphml` e `.gml` também são aceitos diretamente (usa `networkx.read_graphml` / `read_gml`).

### Exemplo de entrada mínima

`arvore1.txt`:
```
vertices:
a
b
c

arestas:
a b
b c
```

```bash
python main.py arvore1.txt arvore1.txt --algoritmo hausdorff --verbose
```

### Uso como biblioteca

```python
from structures.Grafo import Grafo
from algorithms.hausdorff import HausdorffDistanceBetweenTrees

t1 = Grafo()
t1.adicionar_multiplos_vertices(['a', 'b', 'c'])
t1.adicionar_multiplas_arestas([('a', 'b'), ('b', 'c')])

distancia, _ = HausdorffDistanceBetweenTrees(t1, t1)
print(distancia)  # 0.0
```

## Executando os testes

Os arquivos em `tests/` funcionam tanto como scripts standalone quanto com `pytest`:

```bash
# individualmente
python tests/hausdorff_tests.py

# ou, a partir da raiz do repositório
pytest tests/
```

**Atenção:** alguns arquivos de teste (particularmente `mcs_bbp_algorithm_tests.py` e `protein_iso_tests.py`) incluem testes de estresse e de escalabilidade que podem levar vários minutos e consumir bastante memória em máquinas modestas, por trabalharem com grafos de até milhares de vértices.

## Limitações por algoritmo

Resumo do domínio de validade, tipo de retorno e comportamento fora do domínio de cada implementação. Para a justificativa teórica completa, ver o TCC associado a este repositório.

| Algoritmo | Domínio válido | Tipo de retorno | Timeout / comportamento no limite |
|---|---|---|---|
| **Hausdorff** (`hausdorff.py`) | Árvores (inclui casos triviais de 1–2 vértices) | Distância numérica exata | Sem timeout explícito; fallback interno de distância aproximada em caso de erro no cálculo de atribuição ótima |
| **MCS/BBP** (`mcs_bbp_algorithm.py`) | Grafos outerplanares | Subgrafo comum máximo, exato dentro do domínio | Timeout de 30s, limite de 10.000 vértices. **Não valida outerplanaridade da entrada** — grafos fora do domínio (ex: grafo de Petersen) são processados sem erro, mas o resultado não tem garantia de corretude |
| **Isomorfismo molecular** (`mol_graph_iso.py`) | Grafos rotulados (moléculas) | Decisão de isomorfismo/automorfismo | Busca exata apenas para grafos com ≤8 vértices; entre 9–20 vértices usa refinamento heurístico por classes de equivalência; acima de 20, modo conservador (poucas simetrias detectadas). Timeout de 5s |
| **He et al.** (`polynomial_algorithm_undirected.py`) | Grafos simples não direcionados | Decisão via cadeia de invariantes necessários (grau, somas de potências, autovalores, SVD) | Nunca prova isomorfismo — apenas descarta pares que **não** são isomorfos. Grafos não isomorfos podem, em tese, compartilhar todos os invariantes testados. Sem timeout |
| **Guo et al.** (`protein_iso.py`) | Grafos genéricos (usado aqui para proteínas) | Distância quantitativa contínua (não é uma decisão binária de isomorfismo) | Distância 0.0 indica identidade estrutural completa; valores maiores indicam graus crescentes de dissimilaridade. Sem timeout |

## Licença

Este projeto está licenciado sob a licença MIT — veja o arquivo [LICENSE](LICENSE) para o texto completo.

> Ajuste o nome do titular dos direitos e o ano no arquivo `LICENSE` antes de publicar, e confirme que a licença MIT é realmente a que você pretende usar (é a opção mais comum para projetos acadêmicos abertos, mas outras como Apache 2.0 ou BSD também são compatíveis com a afirmação de código aberto).

## Referências

- FAULON, J.-L. Isomorphism, automorphism partitioning, and canonical labeling can be solved in polynomial-time for molecular graphs. *Journal of Chemical Information and Computer Sciences*, v. 38, n. 3, p. 432–444, 1998.
- SCHIETGAT, L.; RAMON, J.; BRUYNOOGHE, M. A polynomial-time maximum common subgraph algorithm for outerplanar graphs and its application to chemoinformatics. *Annals of Mathematics and Artificial Intelligence*, v. 69, n. 4, p. 343–376, 2013.
- GUO, M. et al. A new technique in protein structure quantitative identification. *Procedia Computer Science*, v. 214, p. 1546–1553, 2022.
- KELENC, A. Determining the Hausdorff distance between trees in polynomial time. (2021).
- HE, J. et al. A polynomial-time algorithm for simple undirected graph isomorphism. *Concurrency and Computation: Practice and Experience*, v. 33, n. 7, 2021.