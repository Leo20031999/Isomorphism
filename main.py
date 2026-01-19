import argparse
import sys
import os
from pathlib import Path
import traceback

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

try:
    from structures.Grafo import Grafo
except ImportError:
    print("ERRO: Classe Grafo não encontrada. Certifique-se de que Grafo.py existe no diretório atual.")
    sys.exit(1)

try:
    from algorithms.polynomial_algorithm_undirected import are_isomorphic
    from algorithms.mcs_bbp_algorithm import calcular_mcs_outerplanar
    from algorithms.mol_graph_iso import isomorfismo_molecular, automorfismos_moleculares
    from algorithms.protein_iso import ProteinGraphDistance
    from algorithms.hausdorff import HausdorffDistanceBetweenTrees, is_tree
    from structures.Grafo import Grafo
except ImportError as e:
    print(f"ERRO: Algum algoritmo não pôde ser importado: {e}")
    print("Certifique-se de que todos os arquivos de algoritmo estão no diretório atual.")
    sys.exit(1)


def carregar_grafo(arquivo: str, formato: str = 'auto') -> Grafo:
    """
    Carrega um grafo a partir de arquivo com suporte a múltiplos formatos
    """
    path = Path(arquivo)
    if not path.exists():
        raise FileNotFoundError(f"Arquivo não encontrado: {arquivo}")

    grafo = Grafo()

    if formato == 'auto':
        extensao = path.suffix.lower()
        if extensao in ['.graphml', '.xml']:
            formato = 'graphml'
        elif extensao == '.gml':
            formato = 'gml'
        else:
            formato = 'txt'

    try:
        if formato == 'graphml':
            import networkx as nx
            nx_graph = nx.read_graphml(arquivo)

            for node in nx_graph.nodes():
                atributos = dict(nx_graph.nodes[node])
                rotulo = atributos.pop('label', None) if 'label' in atributos else str(node)
                grafo.adicionar_vertice(node, rotulo, atributos)

            for u, v in nx_graph.edges():
                atributos = dict(nx_graph.edges[u, v])
                rotulo = atributos.pop('label', None) if 'label' in atributos else None
                grafo.adicionar_aresta(u, v, rotulo, atributos)

        elif formato == 'gml':
            import networkx as nx
            nx_graph = nx.read_gml(arquivo)

            for node in nx_graph.nodes():
                atributos = dict(nx_graph.nodes[node])
                rotulo = atributos.pop('label', None) if 'label' in atributos else str(node)
                grafo.adicionar_vertice(node, rotulo, atributos)

            for u, v in nx_graph.edges():
                atributos = dict(nx_graph.edges[u, v])
                rotulo = atributos.pop('label', None) if 'label' in atributos else None
                grafo.adicionar_aresta(u, v, rotulo, atributos)

        else:  
            with open(arquivo, 'r', encoding='utf-8') as f:
                conteudo = f.read().strip()

            linhas = conteudo.split('\n')
            modo = None

            for linha in linhas:
                linha = linha.strip()
                if not linha or linha.startswith('#'):
                    continue

                if linha.lower() == 'vertices:':
                    modo = 'vertices'
                    continue
                elif linha.lower() == 'arestas:':
                    modo = 'arestas'
                    continue
                elif linha.lower() == 'rotulos:':
                    modo = 'rotulos'
                    continue
                elif linha.lower() == 'atributos:':
                    modo = 'atributos'
                    continue

                if modo == 'vertices':
                    if ':' in linha:
                        v, rotulo = linha.split(':', 1)
                        grafo.adicionar_vertice(v.strip(), rotulo.strip())
                    else:
                        grafo.adicionar_vertice(linha, linha)

                elif modo == 'arestas':
                    partes = linha.split()
                    if len(partes) >= 2:
                        u, v = partes[0], partes[1]
                        rotulo = None
                        if len(partes) > 2:
                            if ':' in partes[2]:
                                rotulo = partes[2].split(':', 1)[1]
                            else:
                                rotulo = partes[2]
                        grafo.adicionar_aresta(u, v, rotulo)

                elif modo == 'rotulos':
                    if ':' in linha:
                        v, rotulo = linha.split(':', 1)
                        if grafo.existe_vertice(v.strip()):
                            grafo.set_rotulo_vertice(v.strip(), rotulo.strip())

                elif modo == 'atributos':
                    if ' ' in linha:
                        partes = linha.split()
                        v = partes[0]
                        if grafo.existe_vertice(v):
                            atributos = {}
                            for attr in partes[1:]:
                                if ':' in attr:
                                    chave, valor = attr.split(':', 1)
                                    atributos[chave] = valor
                            current_attrs = grafo.get_atributos_vertice(v)
                            current_attrs.update(atributos)

    except Exception as e:
        raise ValueError(f"Erro ao carregar arquivo {arquivo} (formato: {formato}): {str(e)}")

    if len(grafo.vertices()) == 0:
        raise ValueError(f"Arquivo {arquivo} não contém um grafo válido")

    return grafo

def comparar_grafos(grafo1: Grafo, grafo2: Grafo, algoritmo: str, verbose: bool = True) -> dict:
    """
    Compara dois grafos usando o algoritmo especificado
    """
    resultados = {
        'algoritmo': algoritmo,
        'sucesso': False,
        'resultado': None,
        'mensagem': ''
    }

    try:
        if algoritmo == 'isomorfismo':
            # Algoritmo HE et al., 2019
            resultado = are_isomorphic(grafo1, grafo2)
            resultados.update({
                'sucesso': True,
                'resultado': resultado,
                'mensagem': 'Os grafos são isomorfos' if resultado else 'Os grafos NÃO são isomorfos'
            })

        elif algoritmo == 'mcs':
            # Algoritmo SCHIETGAT et al., 2013
            mcs_grafo, tamanho = calcular_mcs_outerplanar(grafo1, grafo2)
            resultados.update({
                'sucesso': True,
                'resultado': {
                    'grafo_mcs': mcs_grafo,
                    'tamanho': tamanho
                },
                'mensagem': f'MCS encontrado com tamanho {tamanho:.2f}'
            })

        elif algoritmo == 'isomorfismo_molecular':
            # Algoritmo FAULON, 1998
            resultado = isomorfismo_molecular(grafo1, grafo2)
            resultados.update({
                'sucesso': True,
                'resultado': resultado,
                'mensagem': 'Os grafos moleculares são isomorfos' if resultado else 'Os grafos moleculares NÃO são isomorfos'
            })

        elif algoritmo == 'distancia_proteina':
            # Algoritmo GUO et al., 2022
            calculador = ProteinGraphDistance()
            distancia = calculador.quantitative_distance(grafo1, grafo2, verbose=verbose)
            resultados.update({
                'sucesso': True,
                'resultado': distancia,
                'mensagem': f'Distância quantitativa: {distancia:.4f}'
            })

        elif algoritmo == 'hausdorff':
            # Algoritmo KELENC, 2021 (para árvores)
            if not is_tree(grafo1) or not is_tree(grafo2):
                resultados['mensagem'] = 'AVISO: Um ou ambos os grafos não são árvores. Resultado pode ser impreciso.'

            distancia, _ = HausdorffDistanceBetweenTrees(grafo1, grafo2)
            resultados.update({
                'sucesso': True,
                'resultado': distancia,
                'mensagem': f'Distância de Hausdorff: {distancia:.4f}'
            })

        elif algoritmo == 'automorfismos':
            # Automorfismos para grafo1
            automorfismos = automorfismos_moleculares(grafo1)
            resultados.update({
                'sucesso': True,
                'resultado': automorfismos,
                'mensagem': f'Encontrados {len(automorfismos)} automorfismos'
            })

        elif algoritmo == 'todos':
            # Executar todos os algoritmos
            todos_resultados = {}

            # Isomorfismo
            try:
                iso_result = are_isomorphic(grafo1, grafo2)
                todos_resultados['isomorfismo'] = iso_result
            except Exception as e:
                todos_resultados['isomorfismo'] = f'Erro: {str(e)}'

            # MCS
            try:
                mcs_grafo, tamanho_mcs = calcular_mcs_outerplanar(grafo1, grafo2)
                todos_resultados['mcs'] = tamanho_mcs
            except Exception as e:
                todos_resultados['mcs'] = f'Erro: {str(e)}'

            # Isomorfismo Molecular
            try:
                iso_mol = isomorfismo_molecular(grafo1, grafo2)
                todos_resultados['isomorfismo_molecular'] = iso_mol
            except Exception as e:
                todos_resultados['isomorfismo_molecular'] = f'Erro: {str(e)}'

            # Distância Proteica
            try:
                calculador = ProteinGraphDistance()
                distancia_prot = calculador.quantitative_distance(grafo1, grafo2, verbose=False)
                todos_resultados['distancia_proteina'] = distancia_prot
            except Exception as e:
                todos_resultados['distancia_proteina'] = f'Erro: {str(e)}'

            # Hausdorff
            try:
                if is_tree(grafo1) and is_tree(grafo2):
                    dist_hausdorff, _ = HausdorffDistanceBetweenTrees(grafo1, grafo2)
                    todos_resultados['hausdorff'] = dist_hausdorff
                else:
                    todos_resultados['hausdorff'] = 'Não aplicável (não são árvores)'
            except Exception as e:
                todos_resultados['hausdorff'] = f'Erro: {str(e)}'

            resultados.update({
                'sucesso': True,
                'resultado': todos_resultados,
                'mensagem': 'Comparação completa realizada'
            })

        else:
            resultados['mensagem'] = f'Algoritmo desconhecido: {algoritmo}'

    except Exception as e:
        resultados['mensagem'] = f'Erro ao executar algoritmo {algoritmo}: {str(e)}'
        if verbose:
            traceback.print_exc()

    return resultados


def exibir_resultados(resultados: dict, verbose: bool = True):
    """
    Exibe os resultados de forma organizada
    """
    print(f"\n=== RESULTADOS - {resultados['algoritmo'].upper()} ===")

    if not resultados['sucesso']:
        print(f"❌ {resultados['mensagem']}")
        return

    print(f"✅ {resultados['mensagem']}")

    if verbose and resultados['resultado'] is not None:
        if resultados['algoritmo'] == 'todos':
            print("\nDetalhes:")
            for algo, resultado in resultados['resultado'].items():
                print(f"  {algo}: {resultado}")
        elif resultados['algoritmo'] == 'automorfismos':
            automorfismos = resultados['resultado']
            if automorfismos:
                print(f"\nAutomorfismos encontrados ({len(automorfismos)}):")
                for i, automorfismo in enumerate(automorfismos[:5]):  # Mostrar apenas os primeiros 5
                    print(f"  {i + 1}: {automorfismo}")
                if len(automorfismos) > 5:
                    print(f"  ... e mais {len(automorfismos) - 5} automorfismos")
            else:
                print("  Nenhum automorfismo encontrado")
        elif resultados['algoritmo'] == 'mcs':
            mcs_info = resultados['resultado']
            print(
                f"  Grafo MCS: {len(mcs_info['grafo_mcs'].vertices())} vértices, {len(mcs_info['grafo_mcs'].arestas())} arestas")
            print(f"  Tamanho ponderado: {mcs_info['tamanho']:.2f}")
        else:
            print(f"Resultado detalhado: {resultados['resultado']}")


def main():
    """
    Função principal com interface de linha de comando
    """
    parser = argparse.ArgumentParser(
        description='Comparação de Grafos - Ferramenta de Análise',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
    Exemplos de uso:
      python main.py grafo1.txt grafo2.txt --algoritmo isomorfismo
      python main.py proteina1.graphml proteina2.graphml --algoritmo todos
      python main.py molecula1.gml molecula2.gml --algoritmo mcs --verbose
    
    Algoritmos disponíveis:
      • isomorfismo        - Verifica isomorfismo (HE et al., 2019)
      • mcs                - Maximum Common Subgraph (SCHIETGAT et al., 2013)  
      • isomorfismo_molecular - Isomorfismo para grafos moleculares (FAULON, 1998)
      • distancia_proteina - Distância quantitativa para proteínas (GUO et al., 2022)
      • hausdorff          - Distância de Hausdorff para árvores (KELENC, 2021)
      • automorfismos      - Calcula automorfismos do primeiro grafo
      • todos              - Executa todos os algoritmos
            '''
    )

    parser.add_argument('arquivo1', help='Primeiro arquivo contendo o grafo')
    parser.add_argument('arquivo2', help='Segundo arquivo contendo o grafo')

    parser.add_argument(
        '--algoritmo',
        choices=['isomorfismo', 'mcs', 'isomorfismo_molecular',
                 'distancia_proteina', 'hausdorff', 'automorfismos', 'todos'],
        default='todos',
        help='Algoritmo de comparação a ser usado (padrão: todos)'
    )

    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Exibe informações detalhadas durante a execução'
    )

    parser.add_argument(
        '--formato',
        choices=['auto', 'graphml', 'gml', 'txt'],
        default='auto',
        help='Formato do arquivo (padrão: auto-detect)'
    )

    args = parser.parse_args()

    # Banner inicial
    print("=" * 60)
    print("          SISTEMA DE COMPARAÇÃO DE GRAFOS")
    print("=" * 60)
    print(f"Arquivo 1: {args.arquivo1}")
    print(f"Arquivo 2: {args.arquivo2}")
    print(f"Algoritmo: {args.algoritmo}")
    print(f"Formato: {args.formato}")
    print(f"Verbose: {args.verbose}")
    print("-" * 60)

    try:
        # Carregar grafos
        if args.verbose:
            print("Carregando grafos...")

        grafo1 = carregar_grafo(args.arquivo1, args.formato)
        grafo2 = carregar_grafo(args.arquivo2, args.formato)

        # Estatísticas dos grafos
        if args.verbose:
            print(f"Grafo 1: {len(grafo1.vertices())} vértices, {len(grafo1.arestas())} arestas")
            print(f"Grafo 2: {len(grafo2.vertices())} vértices, {len(grafo2.arestas())} arestas")

            # Verificar se são árvores
            arvore1 = is_tree(grafo1)
            arvore2 = is_tree(grafo2)
            print(f"Grafo 1 é árvore: {arvore1}")
            print(f"Grafo 2 é árvore: {arvore2}")

            # Verificar se são outerplanar
            try:
                outer1 = grafo1.eh_outerplanar()
                outer2 = grafo2.eh_outerplanar()
                print(f"Grafo 1 é outerplanar: {outer1}")
                print(f"Grafo 2 é outerplanar: {outer2}")
            except:
                print("Não foi possível verificar outerplanaridade")

        # Comparar grafos
        if args.verbose:
            print(f"\nExecutando algoritmo: {args.algoritmo}")

        resultados = comparar_grafos(grafo1, grafo2, args.algoritmo, args.verbose)

        # Exibir resultados
        exibir_resultados(resultados, args.verbose)

    except FileNotFoundError as e:
        print(f"❌ ERRO: {e}")
        sys.exit(1)
    except ValueError as e:
        print(f"❌ ERRO: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"❌ ERRO inesperado: {e}")
        if args.verbose:
            traceback.print_exc()
        sys.exit(1)

    print("\n" + "=" * 60)
    print("Comparação concluída!")
    print("=" * 60)


if __name__ == "__main__":
    main()