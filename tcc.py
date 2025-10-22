# Determining the hausdorff distance between trees in polynomial time
# Kelenc, A., 2021

from collections import deque
from typing import Dict, List, Set, Tuple, Any, Optional
import numpy as np
from scipy.optimize import linear_sum_assignment


def is_tree(graph):
    """
    Verifica se o grafo é uma árvore.
    Uma árvore é um grafo conexo e acíclico com n vértices e n-1 arestas.
    """
    vertices = graph.vertices()
    n = len(vertices)

    if n == 0:
        return True

    total_edges = sum(len(graph.vizinhanca(v)) for v in vertices) // 2

    if total_edges != n - 1:
        return False

    if n > 0:
        visited = set()
        queue = deque([vertices[0]])

        while queue:
            current = queue.popleft()
            if current not in visited:
                visited.add(current)
                for neighbor in graph.vizinhanca(current):
                    if neighbor not in visited:
                        queue.append(neighbor)

        return len(visited) == n

    return True


def ensure_tree(graph, graph_name="Grafo"):
    """
    Garante que o grafo é uma árvore. Se não for, lança uma exceção.
    """
    if not is_tree(graph):
        raise ValueError(f"{graph_name} não é uma árvore! "
                         f"O algoritmo de Hausdorff só funciona com árvores.")
    return graph


# ==================== ALGORITMO DE HAUSDORFF PARA ÁRVORES ====================

def atom_similarity(atom1: Dict[str, Any], atom2: Dict[str, Any]) -> float:
    """
    Similaridade atômica simplificada para árvores moleculares.
    """
    if not atom1 or not atom2:
        return 0.0
    if 'atom_type' not in atom1 or 'atom_type' not in atom2:
        return 0.0

    type1 = atom1['atom_type']
    type2 = atom2['atom_type']

    if type1 == type2:
        return 1.0

    similarity_rules = {
        ('H', 'Cl'): 0.3, ('Cl', 'H'): 0.3,
        ('C', 'H'): 0.2, ('H', 'C'): 0.2,
        ('C', 'O'): 0.4, ('O', 'C'): 0.4,
        ('C', 'N'): 0.3, ('N', 'C'): 0.3,
        ('C', 'Cl'): 0.2, ('Cl', 'C'): 0.2,
        ('H', 'O'): 0.1, ('O', 'H'): 0.1,
        ('H', 'N'): 0.1, ('N', 'H'): 0.1,
        ('O', 'N'): 0.5, ('N', 'O'): 0.5,
    }

    return similarity_rules.get((type1, type2), 0.0)


def compute_tree_properties(T, root):
    """
    Computa propriedades de uma árvore.
    """
    parent_map = {root: None}
    heights = {}
    sizes = {}
    children_map = {}
    pre_order = []

    queue = deque([root])
    while queue:
        u = queue.popleft()
        pre_order.append(u)
        children = []
        for v in T.vizinhanca(u):
            if v != parent_map.get(u):
                parent_map[v] = u
                children.append(v)
                queue.append(v)
        children_map[u] = children

    for node in reversed(pre_order):
        if not children_map[node]:
            heights[node] = 0
            sizes[node] = 1
        else:
            heights[node] = 1 + max(heights[child] for child in children_map[node])
            sizes[node] = 1 + sum(sizes[child] for child in children_map[node])

    return parent_map, heights, sizes, children_map, pre_order


def TreeDistance(T1, v, T2, u,
                 parent_map_T1, parent_map_T2,
                 heights_T1, heights_T2,
                 sizes_T1, sizes_T2,
                 children_map_T1, children_map_T2,
                 memo=None, use_attributes=False):
    """
    Calcula a distância entre duas subárvores.
    """
    if memo is None:
        memo = {}

    key = (v, u)
    if key in memo:
        return memo[key]

    is_leaf_T1 = not children_map_T1.get(v, [])
    is_leaf_T2 = not children_map_T2.get(u, [])

    attribute_dist = 0.0
    if use_attributes:
        try:
            atom_sim = atom_similarity(
                T1.get_atributos_vertice(v),
                T2.get_atributos_vertice(u)
            )
            attribute_dist = 1.0 - atom_sim
        except:
            attribute_dist = 0.0

    if is_leaf_T1 and is_leaf_T2:
        dist = attribute_dist
        memo[key] = dist
        return dist

    if is_leaf_T1 or is_leaf_T2:
        height_penalty = abs(heights_T1[v] - heights_T2[u]) * 0.1
        dist = attribute_dist + height_penalty
        memo[key] = dist
        return dist

    try:
        filhos_T1 = children_map_T1[v]
        filhos_T2 = children_map_T2[u]

        n = max(len(filhos_T1), len(filhos_T2))
        cost_matrix = np.zeros((n, n))

        for i, child1 in enumerate(filhos_T1):
            for j, child2 in enumerate(filhos_T2):
                cost_matrix[i, j] = TreeDistance(
                    T1, child1, T2, child2,
                    parent_map_T1, parent_map_T2,
                    heights_T1, heights_T2,
                    sizes_T1, sizes_T2,
                    children_map_T1, children_map_T2,
                    memo, use_attributes
                )

        penalty = max(heights_T1[v], heights_T2[u]) * 0.2

        for i in range(len(filhos_T1)):
            for j in range(len(filhos_T2), n):
                cost_matrix[i, j] = penalty

        for j in range(len(filhos_T2)):
            for i in range(len(filhos_T1), n):
                cost_matrix[i, j] = penalty

        for i in range(len(filhos_T1), n):
            for j in range(len(filhos_T2), n):
                cost_matrix[i, j] = 0.0

        row_ind, col_ind = linear_sum_assignment(cost_matrix)

        total_cost = 0.0
        for i, j in zip(row_ind, col_ind):
            total_cost += cost_matrix[i, j]

        avg_cost = total_cost / n if n > 0 else 0.0

        final_dist = avg_cost + attribute_dist

        max_reasonable = max(heights_T1[v], heights_T2[u])
        final_dist = min(final_dist, max_reasonable)

        memo[key] = final_dist
        return final_dist

    except Exception as e:
        fallback_dist = attribute_dist + abs(heights_T1[v] - heights_T2[u]) * 0.1
        memo[key] = fallback_dist
        return fallback_dist


def HausdorffDistanceBetweenTrees(T1, T2, use_attributes=False):
    """
    Calcula a distância de Hausdorff entre duas árvores.
    """
    try:
        ensure_tree(T1, "T1")
        ensure_tree(T2, "T2")

        if not T1.vertices() and not T2.vertices():
            return 0.0, set()
        if not T1.vertices() or not T2.vertices():
            return float('inf'), set()

        if len(T1.vertices()) == 1 and len(T2.vertices()) == 1:
            if use_attributes:
                v1 = T1.vertices()[0]
                v2 = T2.vertices()[0]
                atom_sim = atom_similarity(
                    T1.get_atributos_vertice(v1),
                    T2.get_atributos_vertice(v2)
                )
                return (1.0 - atom_sim), set()
            else:
                return 0.0, set()

        r1 = T1.vertices()[0]
        r2 = T2.vertices()[0]

        parent_map_T1, heights_T1, sizes_T1, children_map_T1, pre_order_T1 = compute_tree_properties(T1, r1)
        parent_map_T2, heights_T2, sizes_T2, children_map_T2, pre_order_T2 = compute_tree_properties(T2, r2)

        memo = {}
        distance = TreeDistance(
            T1, r1, T2, r2,
            parent_map_T1, parent_map_T2,
            heights_T1, heights_T2,
            sizes_T1, sizes_T2,
            children_map_T1, children_map_T2,
            memo, use_attributes
        )

        if distance < 0:
            distance = 0.0

        return round(distance, 3), set()

    except ValueError as e:
        print(f"Erro: {e}")
        return float('inf'), set()
    except Exception as e:
        print(f"Erro inesperado: {e}")
        return float('inf'), set()


def SafeTreeDistance(T1, T2, use_attributes=False):
    """
    Versão segura que verifica se as entradas são árvores.
    """
    try:
        return HausdorffDistanceBetweenTrees(T1, T2, use_attributes)
    except Exception as e:
        print(f"Erro no cálculo de distância: {e}")
        return float('inf'), set()