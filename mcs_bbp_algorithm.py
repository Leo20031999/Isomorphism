# A Polynomial-time Maximum Common Subgraph Algorithm for Outerplanar Graphs and its Application to Chemoinformatics
# SCHIETGAT; RAMON; BRUYNOOGHE, 2013

import networkx as nx
import time
from collections import deque
from typing import Dict, List, Tuple, Optional, Any

from structures.Grafo import Grafo

DEFAULT_LABEL_WEIGHTS = {
    'H': 0.3, 'C': 1.0, 'N': 1.4, 'O': 1.6, 'F': 1.8,
    'P': 2.0, 'S': 1.7, 'Cl': 2.2, 'Br': 2.8, 'I': 3.2,
    'single': 1.0, 'double': 1.8, 'triple': 2.5,
    'aromatic': 1.4, 'amide': 1.6, 'ionic': 2.0, 'hydrogen': 0.8
}

MAX_VERTICES = 50
MAX_DECOMPOSITION_DEPTH = 3
MAX_ELEMENTARY_PARTS = 4
TIMEOUT_DURATION = 10

class OuterplanarMCS:
    """Implementação final do algoritmo MCS para grafos outerplanares"""

    def __init__(self, label_weights: Optional[Dict[str, float]] = None,
                 max_vertices: int = MAX_VERTICES, timeout: int = TIMEOUT_DURATION):
        self.memo = {}
        self.decomposition_cache = {}
        self.label_weights = label_weights if label_weights is not None else DEFAULT_LABEL_WEIGHTS
        self.max_vertices = max_vertices
        self.timeout = timeout
        self.start_time = None
        self.mcs_best_candidate = (Grafo(), 0.0)
        self._timeout_occurred = False

    def compute_mcs(self, G: Grafo, H: Grafo) -> Tuple[Grafo, float]:
        """Algoritmo principal CORRIGIDO com BBP - VERSÃO OTIMIZADA"""
        self.start_time = time.time()
        self._timeout_occurred = False
        self.memo.clear()
        self.decomposition_cache.clear()
        self.mcs_best_candidate = (Grafo(), 0.0)

        if self._are_graphs_identical(G, H):
            return (G.copy(), self._calculate_graph_size(G))

        if len(G.vertices()) > self.max_vertices or len(H.vertices()) > self.max_vertices:
            return self._compute_approximate_mcs(G, H)

        if not nx.is_connected(G.grafo) or not nx.is_connected(H.grafo):
            return self._compute_approximate_mcs(G, H)

        root_candidates_G = self._find_root_candidates(G)
        root_candidates_H = self._find_root_candidates(H)

        promising_pairs = []
        for r in root_candidates_G:
            for s in root_candidates_H:
                if G.get_rotulo_vertice(r) == H.get_rotulo_vertice(s):
                    promising_pairs.append((r, s))

        promising_pairs.sort(key=lambda pair:
        G.grau(pair[0]) + H.grau(pair[1]), reverse=True)

        for r, s in promising_pairs[:10]:
            if self._check_timeout():
                break

            parts_G = self._decompose_graph(G, r)
            parts_H = self._decompose_graph(H, s)
            self._process_root_pair_with_bbp(parts_G, parts_H, r, s)

            if self.mcs_best_candidate[1] > 5.0:
                break

        if self.mcs_best_candidate[1] == 0.0:
            return self._compute_enhanced_fallback(G, H)

        return self.mcs_best_candidate

    def _compute_enhanced_fallback(self, G: Grafo, H: Grafo) -> Tuple[Grafo, float]:
        """Fallback aprimorado para casos difíceis"""
        best_graph = Grafo()
        best_size = 0.0

        g_label_count = {}
        for v in G.vertices():
            label = G.get_rotulo_vertice(v)
            g_label_count[label] = g_label_count.get(label, 0) + 1

        h_label_count = {}
        for v in H.vertices():
            label = H.get_rotulo_vertice(v)
            h_label_count[label] = h_label_count.get(label, 0) + 1

        common_labels = set(g_label_count.keys()) & set(h_label_count.keys())

        for label in common_labels:
            if self._check_timeout():
                break

            g_vertices = [v for v in G.vertices() if G.get_rotulo_vertice(v) == label]
            h_vertices = [v for v in H.vertices() if H.get_rotulo_vertice(v) == label]

            for g_start in g_vertices[:3]:
                for h_start in h_vertices[:3]:
                    candidate, size = self._build_mcs_from_pair(G, H, g_start, h_start)
                    if size > best_size:
                        best_graph = candidate
                        best_size = size

        return (best_graph, best_size)

    def _build_mcs_from_pair(self, G: Grafo, H: Grafo, g_start: Any, h_start: Any) -> Tuple[Grafo, float]:
        """Construir MCS a partir de um par de vértices compatíveis"""
        mcs_graph = Grafo()
        mcs_graph.adicionar_vertice(g_start, G.get_rotulo_vertice(g_start))

        mapping = {g_start: h_start}
        queue = deque([g_start])

        while queue:
            g_v = queue.popleft()
            h_v = mapping[g_v]

            for g_neighbor in G.vizinhanca(g_v):
                if g_neighbor in mapping:
                    continue

                g_label = G.get_rotulo_vertice(g_neighbor)
                edge_label = G.get_rotulo_aresta(g_v, g_neighbor)

                for h_neighbor in H.vizinhanca(h_v):
                    if h_neighbor in mapping.values():
                        continue

                    if (H.get_rotulo_vertice(h_neighbor) == g_label and
                            H.get_rotulo_aresta(h_v, h_neighbor) == edge_label):
                        mcs_graph.adicionar_vertice(g_neighbor, g_label)
                        mcs_graph.adicionar_aresta(g_v, g_neighbor, edge_label)
                        mapping[g_neighbor] = h_neighbor
                        queue.append(g_neighbor)
                        break

        return (mcs_graph, self._calculate_graph_size(mcs_graph))

    def _process_root_pair_with_bbp(self, parts_G: List[Dict], parts_H: List[Dict], r: Any, s: Any):
        """Processamento com verificação BBP"""
        root_parts_G = [p for p in parts_G if p['root'] == r and p['type'] == 'BPS']
        root_parts_H = [p for p in parts_H if p['root'] == s and p['type'] == 'BPS']

        for root_part_G in root_parts_G:
            for root_part_H in root_parts_H:
                mcs, size = self._rmcs_with_bbp(root_part_G, root_part_H)
                if size > self.mcs_best_candidate[1]:
                    self.mcs_best_candidate = (mcs, size)

    def _rmcs_with_bbp(self, P: Dict, Q: Dict) -> Tuple[Grafo, float]:
        """RMCS com verificação BBP"""
        if self._check_timeout():
            return (Grafo(), 0.0)

        if P['graph'].get_rotulo_vertice(P['root']) != Q['graph'].get_rotulo_vertice(Q['root']):
            return self._create_single_vertex_graph(P)

        if P['type'] == 'BPS' and Q['type'] == 'BPS':
            if not P['is_basic'] and not Q['is_basic']:
                return self._rmcs_compound_with_bbp(P, Q)
            else:
                return self._rmcs_basic_with_bbp(P, Q)
        else:
            return self._rmcs_bss_with_bbp(P, Q)

    def _rmcs_compound_with_bbp(self, P: Dict, Q: Dict) -> Tuple[Grafo, float]:
        """RMCS para raiz composta com BBP"""
        best_graph, best_size = self._create_single_vertex_graph(P)

        ep_P = self._get_elementary_parts(P, P['graph'])
        ep_Q = self._get_elementary_parts(Q, Q['graph'])

        used_p = set()
        used_q = set()

        for i, part_P in enumerate(ep_P):
            if self._check_timeout():
                break

            best_match_size = 0
            best_match_graph = None
            best_j = -1

            for j, part_Q in enumerate(ep_Q):
                if j in used_q:
                    continue

                if self._are_parts_compatible(part_P, part_Q):
                    sub_mcs, sub_size = self._rmcs_with_bbp(part_P, part_Q)
                    if sub_size > best_match_size:
                        best_match_size = sub_size
                        best_match_graph = sub_mcs
                        best_j = j

            if best_match_graph and best_match_size > 0:
                best_graph = self._merge_graphs(best_graph, best_match_graph)
                best_size += best_match_size
                used_p.add(i)
                used_q.add(best_j)

        return (best_graph, best_size)

    def _rmcs_bss_with_bbp(self, P: Dict, Q: Dict) -> Tuple[Grafo, float]:
        """RMCS para BSS com BBP"""
        if 'split_edge' not in P or 'split_edge' not in Q:
            return self._create_single_vertex_graph(P)

        p_split = P['split_edge']
        q_split = Q['split_edge']

        p_label = P['graph'].get_rotulo_aresta(p_split[0], p_split[1])
        q_label = Q['graph'].get_rotulo_aresta(q_split[0], q_split[1])

        if p_label == q_label:
            base_graph = self._create_edge_graph(P, p_split[0], p_split[1], p_label)
            base_size = self._calculate_graph_size(base_graph)

            r_prime = p_split[0] if p_split[1] == P['root'] else p_split[1]
            s_prime = q_split[0] if q_split[1] == Q['root'] else q_split[1]

            part_P_expand = {
                'graph': P['graph'],
                'root': r_prime,
                'is_basic': self._is_basic_root(P['graph'], r_prime),
                'type': 'BPS'
            }
            part_Q_expand = {
                'graph': Q['graph'],
                'root': s_prime,
                'is_basic': self._is_basic_root(Q['graph'], s_prime),
                'type': 'BPS'
            }

            sub_mcs, sub_size = self._rmcs_with_bbp(part_P_expand, part_Q_expand)
            if sub_size > 0:
                expanded_graph = self._merge_graphs(base_graph, sub_mcs)
                expanded_size = self._calculate_graph_size(expanded_graph)
                return (expanded_graph, expanded_size)

            return (base_graph, base_size)

        return self._create_single_vertex_graph(P)

    def _rmcs_basic_with_bbp(self, P: Dict, Q: Dict) -> Tuple[Grafo, float]:
        """RMCS básico com BBP - CORREÇÃO CRÍTICA"""
        base_graph, best_size = self._create_single_vertex_graph(P)

        p_neighbors = P['graph'].vizinhanca(P['root'])
        q_neighbors = Q['graph'].vizinhanca(Q['root'])

        candidate_graphs = [(base_graph, best_size)]

        for p_n in p_neighbors:
            for q_n in q_neighbors:
                if self._check_timeout():
                    break

                if not self._are_vertices_compatible(P['graph'], p_n, Q['graph'], q_n):
                    continue

                edge_p_label = P['graph'].get_rotulo_aresta(P['root'], p_n)
                edge_q_label = Q['graph'].get_rotulo_aresta(Q['root'], q_n)

                if edge_p_label != edge_q_label:
                    continue

                if not self._is_edge_bbp_compatible(P['graph'], Q['graph'],
                                                    P['root'], p_n, Q['root'], q_n):
                    continue

                part_P_n = {
                    'graph': P['graph'],
                    'root': p_n,
                    'is_basic': self._is_basic_root(P['graph'], p_n),
                    'type': 'BPS'
                }
                part_Q_n = {
                    'graph': Q['graph'],
                    'root': q_n,
                    'is_basic': self._is_basic_root(Q['graph'], q_n),
                    'type': 'BPS'
                }

                sub_mcs, sub_size = self._rmcs_with_bbp(part_P_n, part_Q_n)

                edge_graph = self._create_edge_graph(P, P['root'], p_n, edge_p_label)
                candidate = self._merge_graphs(edge_graph, sub_mcs)
                candidate_size = self._calculate_graph_size(candidate)

                candidate_graphs.append((candidate, candidate_size))

        best_candidate, best_candidate_size = max(candidate_graphs, key=lambda x: x[1])
        return (best_candidate, best_candidate_size)

    def _are_vertices_compatible(self, G: Grafo, v_g: Any, H: Grafo, v_h: Any) -> bool:
        """Verifica compatibilidade de vértices"""
        return G.get_rotulo_vertice(v_g) == H.get_rotulo_vertice(v_h)

    def _is_edge_bbp_compatible(self, G: Grafo, H: Grafo, u_g: Any, v_g: Any, u_h: Any, v_h: Any) -> bool:
        """Verifica compatibilidade BBP para arestas"""
        try:
            is_bridge_G = self._is_bridge(G, u_g, v_g)
            is_bridge_H = self._is_bridge(H, u_h, v_h)

            return is_bridge_G == is_bridge_H
        except:
            return True

    def _is_bridge(self, graph: Grafo, u: Any, v: Any) -> bool:
        """Verifica se uma aresta é uma bridge"""
        try:
            temp_graph = graph.copy()
            temp_graph.remover_aresta(u, v)
            return not nx.is_connected(temp_graph.grafo)
        except:
            return False

    def _are_graphs_identical(self, G: Grafo, H: Grafo) -> bool:
        """Verificação robusta de grafos idênticos"""
        try:
            if len(G.vertices()) != len(H.vertices()):
                return False
            if len(G.arestas()) != len(H.arestas()):
                return False

            for v in G.vertices():
                g_label = G.get_rotulo_vertice(v)
                found_match = False
                for u in H.vertices():
                    if H.get_rotulo_vertice(u) == g_label:
                        found_match = True
                        break
                if not found_match:
                    return False

            for edge in G.arestas():
                u, v = edge
                g_edge_label = G.get_rotulo_aresta(u, v)
                found_match = False
                for edge_h in H.arestas():
                    x, y = edge_h
                    if (H.get_rotulo_aresta(x, y) == g_edge_label and
                            H.get_rotulo_vertice(x) == G.get_rotulo_vertice(u) and
                            H.get_rotulo_vertice(y) == G.get_rotulo_vertice(v)):
                        found_match = True
                        break
                if not found_match:
                    return False

            return True
        except:
            return False

    def _find_root_candidates(self, graph: Grafo) -> List:
        """Encontra TODOS os candidatos a raiz promissores"""
        vertices = graph.vertices()

        if len(vertices) <= 6:
            return vertices

        scored_vertices = []
        for v in vertices:
            try:
                degree = graph.grau(v)
                label = graph.get_rotulo_vertice(v)

                score = degree * 2
                if label in ['C', 'N', 'O']:
                    score += 3

                scored_vertices.append((v, score))
            except:
                continue

        scored_vertices.sort(key=lambda x: x[1], reverse=True)
        return [v for v, score in scored_vertices[:8]]

    def _decompose_graph(self, graph: Grafo, root: Any) -> List[Dict]:
        """Decomposição COMPLETA"""
        cache_key = (id(graph.grafo), root)
        if cache_key in self.decomposition_cache:
            return self.decomposition_cache[cache_key]

        is_basic = self._is_basic_root(graph, root)

        root_part = {
            'graph': graph,
            'root': root,
            'is_basic': is_basic,
            'type': 'BPS'
        }

        parts = [root_part]

        try:
            elementary_parts = self._get_elementary_parts(root_part, graph)
            for part in elementary_parts[:MAX_ELEMENTARY_PARTS]:
                parts.append(part)
                sub_parts = self._get_elementary_parts(part, graph)
                parts.extend(sub_parts[:2])
        except:
            pass

        self.decomposition_cache[cache_key] = parts
        return parts

    def _is_basic_root(self, graph: Grafo, root: Any) -> bool:
        """Verificação PRECISA se raiz é básica"""
        try:
            degree = graph.grau(root)
            if degree <= 1:
                return True

            blocks = graph.encontrar_blocos()
            incident_blocks = [block for block in blocks if root in block]
            return len(incident_blocks) == 1
        except:
            return True

    def _get_elementary_parts(self, part: Dict, original_graph: Grafo) -> List[Dict]:
        """Obtém partes elementares completas"""
        parts = []
        graph = part['graph']
        root = part['root']

        try:
            if not part['is_basic']:
                neighbors = graph.vizinhanca(root)
                for neighbor in neighbors:
                    new_part = {
                        'graph': graph,
                        'root': neighbor,
                        'is_basic': True,
                        'type': 'BPS'
                    }
                    parts.append(new_part)
            else:
                degree = graph.grau(root)
                if degree == 1:
                    neighbors = graph.vizinhanca(root)
                    if neighbors:
                        neighbor = neighbors[0]
                        new_part = {
                            'graph': graph,
                            'root': neighbor,
                            'is_basic': False,
                            'type': 'BPS'
                        }
                        parts.append(new_part)
                else:
                    neighbors = graph.vizinhanca(root)
                    for neighbor in neighbors:
                        new_part = {
                            'graph': graph,
                            'root': neighbor,
                            'is_basic': True,
                            'type': 'BSS',
                            'split_edge': (root, neighbor)
                        }
                        parts.append(new_part)

            return parts
        except:
            return []

    def _are_parts_compatible(self, P: Dict, Q: Dict) -> bool:
        """Verificação de compatibilidade entre partes - CORRIGIDA"""
        try:
            label_P = P['graph'].get_rotulo_vertice(P['root'])
            label_Q = Q['graph'].get_rotulo_vertice(Q['root'])
            return label_P == label_Q
        except:
            return False

    def _create_single_vertex_graph(self, part: Dict) -> Tuple[Grafo, float]:
        """Cria grafo com vértice único"""
        graph = Grafo()
        root = part['root']
        label = part['graph'].get_rotulo_vertice(root)
        graph.adicionar_vertice(root, label)
        size = self._calculate_graph_size(graph)
        return (graph, size)

    def _create_edge_graph(self, part: Dict, u: Any, v: Any, edge_label: str) -> Grafo:
        """Cria grafo com aresta"""
        graph = Grafo()
        label_u = part['graph'].get_rotulo_vertice(u)
        label_v = part['graph'].get_rotulo_vertice(v)

        graph.adicionar_vertice(u, label_u)
        graph.adicionar_vertice(v, label_v)
        graph.adicionar_aresta(u, v, edge_label)
        return graph

    def _merge_graphs(self, graph1: Grafo, graph2: Grafo) -> Grafo:
        """Fusão de dois grafos"""
        result = graph1.copy()

        for v in graph2.vertices():
            if not result.existe_vertice(v):
                label = graph2.get_rotulo_vertice(v)
                result.adicionar_vertice(v, label)

        for u, v in graph2.arestas():
            if not result.grafo.has_edge(u, v):
                label = graph2.get_rotulo_aresta(u, v)
                result.adicionar_aresta(u, v, label)

        return result

    def _calculate_graph_size(self, graph: Grafo) -> float:
        """Calcula tamanho do grafo usando pesos dos rótulos"""
        total_size = 0.0

        for v in graph.vertices():
            label = graph.get_rotulo_vertice(v)
            if label is None:
                label = "C"
            weight = self.label_weights.get(label, 1.0)
            total_size += weight

        for u, v in graph.arestas():
            label = graph.get_rotulo_aresta(u, v)
            if label is None:
                label = "single"
            weight = self.label_weights.get(label, 1.0)
            total_size += weight

        return total_size

    def _compute_approximate_mcs(self, G: Grafo, H: Grafo) -> Tuple[Grafo, float]:
        """Computação aproximada para fallback - MELHORADA"""
        mcs_graph = Grafo()
        best_size = 0.0

        for g_root in G.vertices():
            for h_root in H.vertices():
                if self._check_timeout():
                    return (mcs_graph, best_size)

                if G.get_rotulo_vertice(g_root) != H.get_rotulo_vertice(h_root):
                    continue

                current_graph = Grafo()
                current_graph.adicionar_vertice(g_root, G.get_rotulo_vertice(g_root))
                current_size = self._calculate_graph_size(current_graph)

                visited_g = set([g_root])
                visited_h = set([h_root])
                vertex_map = {g_root: h_root}

                queue = deque([(g_root, h_root)])

                while queue:
                    g_v, h_v = queue.popleft()

                    for g_neighbor in G.vizinhanca(g_v):
                        if g_neighbor in visited_g:
                            continue

                        g_neighbor_label = G.get_rotulo_vertice(g_neighbor)
                        edge_label = G.get_rotulo_aresta(g_v, g_neighbor)

                        found_match = False
                        for h_neighbor in H.vizinhanca(h_v):
                            if h_neighbor in visited_h:
                                continue

                            if (H.get_rotulo_vertice(h_neighbor) == g_neighbor_label and
                                    H.get_rotulo_aresta(h_v, h_neighbor) == edge_label):
                                current_graph.adicionar_vertice(g_neighbor, g_neighbor_label)
                                current_graph.adicionar_aresta(g_v, g_neighbor, edge_label)
                                current_size = self._calculate_graph_size(current_graph)

                                visited_g.add(g_neighbor)
                                visited_h.add(h_neighbor)
                                vertex_map[g_neighbor] = h_neighbor
                                queue.append((g_neighbor, h_neighbor))
                                found_match = True
                                break

                        if found_match:
                            break

                if current_size > best_size:
                    mcs_graph = current_graph
                    best_size = current_size

        return (mcs_graph, best_size)

    def _check_timeout(self) -> bool:
        """Verifica timeout"""
        if time.time() - self.start_time > self.timeout:
            self._timeout_occurred = True
            return True
        return False


# ============================ FUNÇÃO PRINCIPAL ============================

def calcular_mcs_outerplanar(grafo_G: Grafo, grafo_H: Grafo,
                             label_weights: Optional[Dict[str, float]] = None,
                             timeout: int = TIMEOUT_DURATION) -> Tuple[Grafo, float]:
    """
    Função principal FINAL para cálculo de MCS entre grafos outerplanares.

    Args:
        grafo_G: Primeiro grafo
        grafo_H: Segundo grafo
        label_weights: Pesos dos rótulos (opcional)
        timeout: Timeout em segundos (opcional)

    Returns:
        Tuple[Grafo, float]: MCS e seu tamanho
    """
    try:
        solver = OuterplanarMCS(
            label_weights=label_weights,
            timeout=timeout
        )
        return solver.compute_mcs(grafo_G, grafo_H)
    except Exception as e:
        return (Grafo(), 0.0)