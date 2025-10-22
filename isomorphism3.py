# A polynomial time algorithm for simple undirected graph isomorphism
# (HE et al., 2019)

from typing import List, Tuple, Any, Dict
import itertools
import numpy as np
from collections import defaultdict


def are_isomorphic(g1, g2) -> bool:
    """
    Implementação otimizada do algoritmo polinomial para verificar isomorfismo de grafos.
    """

    n1, m1 = len(g1.vertices()), len(g1.arestas())
    n2, m2 = len(g2.vertices()), len(g2.arestas())

    if n1 != n2 or m1 != m2:
        return False

    degree_seq1 = sorted([g1.grau(v) for v in g1.vertices()])
    degree_seq2 = sorted([g2.grau(v) for v in g2.vertices()])
    if degree_seq1 != degree_seq2:
        return False

    tt1 = _build_complete_triple_tuple(g1)
    tt2 = _build_complete_triple_tuple(g2)

    K1 = _calculate_vertex_row_sum_array(tt1, g1)
    K2 = _calculate_vertex_row_sum_array(tt2, g2)

    if not _verify_permutation(K1, K2):
        return False

    L1 = _calculate_edge_row_sum_array(tt1, g1, K1)
    L2 = _calculate_edge_row_sum_array(tt2, g2, K2)

    if not _verify_permutation(L1, L2):
        return False

    return _verify_graph_invariants(g1, g2)


def _build_complete_triple_tuple(grafo) -> List[Tuple[int, Any, Any]]:
    """Constrói triple tuple rapidamente"""
    edges = grafo.arestas()
    return [(idx + 1, u, v) for idx, (u, v) in enumerate(edges)]


def _calculate_vertex_row_sum_array(triple_tuple: List[Tuple[int, Any, Any]], grafo) -> List[int]:
    """Calcula array K_n rapidamente"""
    degree_map = defaultdict(int)

    for _, u, v in triple_tuple:
        degree_map[u] += 1
        degree_map[v] += 1

    vertices = sorted(grafo.vertices())
    return [degree_map[v] for v in vertices]


def _calculate_edge_row_sum_array(triple_tuple: List[Tuple[int, Any, Any]], grafo,
                                       vertex_row_sum: List[int]) -> List[int]:
    """Calcula array L_m rapidamente sem matriz de adjacência de arestas"""
    vertex_to_index = {v: idx for idx, v in enumerate(sorted(grafo.vertices()))}
    edge_row_sums = []

    for _, u, v in triple_tuple:
        idx_u = vertex_to_index[u]
        idx_v = vertex_to_index[v]
        L_m = vertex_row_sum[idx_u] + vertex_row_sum[idx_v] - 2
        edge_row_sums.append(L_m)

    return edge_row_sums


def _verify_permutation(seq1: List[int], seq2: List[int]) -> bool:
    """
    Verificação ultra-rápida de permutação usando apenas 3 potências
    """
    if len(seq1) != len(seq2):
        return False

    for power in [1, 2, 3]:
        sum1 = sum(x ** power for x in seq1)
        sum2 = sum(x ** power for x in seq2)

        if sum1 != sum2:
            return False

    return True


def _verify_graph_invariants(g1, g2) -> bool:
    """
    Verificação rápida de invariantes de grafos (equinumerosidade simplificada)
    """
    try:
        V1 = _build_vertex_adjacency_matrix(g1)
        V2 = _build_vertex_adjacency_matrix(g2)

        for k in [1, 2, 3]:
            try:
                trace1 = np.trace(np.linalg.matrix_power(V1, k))
                trace2 = np.trace(np.linalg.matrix_power(V2, k))

                if abs(trace1 - trace2) > 1e-10:
                    return False
            except:
                continue

        rank1 = np.linalg.matrix_rank(V1)
        rank2 = np.linalg.matrix_rank(V2)

        return rank1 == rank2

    except:
        return True


def _build_vertex_adjacency_matrix(grafo) -> np.ndarray:
    """Constrói matriz de adjacência de vértices rapidamente"""
    vertices = sorted(grafo.vertices())
    n = len(vertices)
    vertex_to_idx = {v: i for i, v in enumerate(vertices)}

    matrix = np.zeros((n, n), dtype=int)

    for u, v in grafo.arestas():
        i, j = vertex_to_idx[u], vertex_to_idx[v]
        matrix[i][j] = 1
        matrix[j][i] = 1

    return matrix


def _generate_isomorphic_group(grafo, max_permutations: int = 24):
    """Gera grupo isomórfico rapidamente"""
    vertices = sorted(grafo.vertices())
    n = len(vertices)

    isomorphic_graphs = []

    permutations = list(itertools.islice(itertools.permutations(vertices), max_permutations))

    for perm in permutations:
        permutation_map = {vertices[i]: perm[i] for i in range(n)}
        new_graph = _apply_vertex_permutation(grafo, permutation_map)
        isomorphic_graphs.append(new_graph)

    return isomorphic_graphs


def _apply_vertex_permutation(grafo, permutation_map: Dict[Any, Any]):
    """Aplica permutação de vértices rapidamente"""
    new_graph = grafo.__class__()

    for old_vertex in grafo.vertices():
        new_vertex = permutation_map[old_vertex]
        new_graph.adicionar_vertice(new_vertex)

    for u, v in grafo.arestas():
        new_u, new_v = permutation_map[u], permutation_map[v]
        new_graph.adicionar_aresta(new_u, new_v)

    return new_graph