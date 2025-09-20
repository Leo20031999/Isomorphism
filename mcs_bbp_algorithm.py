# A polynomial-time maximum common subgraph algorithm for outerplanar graphs and its application to chemoinformatics
# Schietgat L.; Ramon J.; Bruynooghe M., 2013

from collections import deque
import networkx as nx
import numpy as np
from scipy.optimize import linear_sum_assignment
from functools import lru_cache

from structures.Grafo import Grafo


class BlockPreservingSubgraph:
    def __init__(self, graph, root, is_basic_root, original_graph):
        self.graph = graph
        self.root = root
        self.is_basic_root = is_basic_root
        self.original_graph = original_graph
        self.type = "BPS"

    def get_elementary_parts(self):
        parts = []
        if not self.is_basic_root:
            # Compound-root decomposition
            bridges = [e for e in self.original_graph.encontrar_pontes() if self.root in e]

            for bridge in bridges:
                if bridge not in self.graph.arestas():
                    continue

                other = bridge[1] if bridge[0] == self.root else bridge[0]
                new_graph = self._create_subgraph_after_removal([bridge], [other])
                if len(new_graph.vertices()) > 0:
                    parts.append(BlockPreservingSubgraph(new_graph, other, True, self.original_graph))

            # Handle blocks containing the root
            blocks = self.original_graph.encontrar_blocos()
            for block in blocks:
                if self.root in block:
                    block_edges = [e for e in self.graph.arestas() if set(e).issubset(block)]
                    if block_edges:
                        new_graph = Grafo()
                        for v in block:
                            if v in self.graph.vertices():
                                new_graph.adicionar_vertice(v, self.graph.get_rotulo_vertice(v))
                        for u, v in block_edges:
                            new_graph.adicionar_aresta(u, v, self.graph.get_rotulo_aresta(u, v))
                        parts.append(BlockPreservingSubgraph(new_graph, self.root, True, self.original_graph))
        else:
            # Basic-root decomposition
            if self.graph.grau(self.root) == 1:  # Bridge case
                neighbors = list(self.graph.vizinhanca(self.root))
                if neighbors:
                    other = neighbors[0]
                    new_graph = self._create_subgraph_after_removal([(self.root, other)], [other])
                    if len(new_graph.vertices()) > 0:
                        parts.append(BlockPreservingSubgraph(new_graph, other, False, self.original_graph))
            else:  # Block case
                # For outerplanar graphs, find the Hamiltonian cycle
                cycle = encontrar_ciclo_hamiltoniano_outerplanar(self.graph)
                if cycle:
                    for edge in cycle:
                        if self.root in edge:
                            other = edge[0] if edge[1] == self.root else edge[1]
                            bss = BlockSplittingSubgraph(self.graph, other, True, self.original_graph,
                                                         split_edge=edge, orientation='x')
                            parts.append(bss)
        return parts

    def _create_subgraph_after_removal(self, edges_to_remove, vertices_to_keep):
        new_graph = self.graph.copy()
        for u, v in edges_to_remove:
            if new_graph.grafo.has_edge(u, v):
                new_graph.grafo.remove_edge(u, v)

        # Find connected component containing the keep vertices
        components = list(nx.connected_components(new_graph.grafo))
        for comp in components:
            if any(v in comp for v in vertices_to_keep):
                result = Grafo()
                for v in comp:
                    result.adicionar_vertice(v, new_graph.get_rotulo_vertice(v))
                for u, v in new_graph.arestas():
                    if u in comp and v in comp:
                        result.adicionar_aresta(u, v, new_graph.get_rotulo_aresta(u, v))
                return result

        # Return empty graph if no component found
        return Grafo()


class BlockSplittingSubgraph:
    def __init__(self, graph, root, is_basic_root, original_graph, split_edge=None, orientation=None):
        self.graph = graph
        self.root = root
        self.is_basic_root = is_basic_root
        self.original_graph = original_graph
        self.type = "BSS"
        self.split_edge = split_edge
        self.orientation = orientation

    def get_elementary_parts(self):
        parts = []
        if not self.is_basic_root:
            # Compound-root BSS decomposition
            bridges = [e for e in self.original_graph.encontrar_pontes() if
                       self.root in e and e in self.graph.arestas()]

            for bridge in bridges:
                other = bridge[1] if bridge[0] == self.root else bridge[0]
                new_graph = self._create_subgraph_after_removal([bridge], [self.root])
                if len(new_graph.vertices()) > 0:
                    parts.append(BlockPreservingSubgraph(new_graph, self.root, True, self.original_graph))

            # Handle blocks containing the root
            blocks = self.original_graph.encontrar_blocos()
            for block in blocks:
                if self.root in block:
                    block_edges = [e for e in self.graph.arestas() if set(e).issubset(block)]
                    if block_edges:
                        new_graph = Grafo()
                        for v in block:
                            if v in self.graph.vertices():
                                new_graph.adicionar_vertice(v, self.graph.get_rotulo_vertice(v))
                        for u, v in block_edges:
                            new_graph.adicionar_aresta(u, v, self.graph.get_rotulo_aresta(u, v))
                        parts.append(BlockPreservingSubgraph(new_graph, self.root, True, self.original_graph))
        return parts

    def _create_subgraph_after_removal(self, edges_to_remove, vertices_to_keep):
        new_graph = self.graph.copy()
        for u, v in edges_to_remove:
            if new_graph.grafo.has_edge(u, v):
                new_graph.grafo.remove_edge(u, v)

        # Find connected component containing the keep vertices
        components = list(nx.connected_components(new_graph.grafo))
        for comp in components:
            if any(v in comp for v in vertices_to_keep):
                result = Grafo()
                for v in comp:
                    result.adicionar_vertice(v, new_graph.get_rotulo_vertice(v))
                for u, v in new_graph.arestas():
                    if u in comp and v in comp:
                        result.adicionar_aresta(u, v, new_graph.get_rotulo_aresta(u, v))
                return result

        # Return empty graph if no component found
        return Grafo()


class OuterplanarMCS:
    def __init__(self):
        self.memo = {}
        self.decomposition_cache = {}

    def compute_mcs(self, G, H):
        # Step 1: Find root candidates
        root_candidates_G = self.find_root_candidates(G)
        root_candidates_H = self.find_root_candidates(H)

        best_mcs = (Grafo(), 0)

        # Step 2: For each candidate pair
        for r in root_candidates_G:
            parts_G = self.decompose_graph(G, r)
            for s in root_candidates_H:
                if G.get_rotulo_vertice(r) != H.get_rotulo_vertice(s):
                    continue

                parts_H = self.decompose_graph(H, s)

                # Step 3: Compute RMCS for root pair
                mcs, size = self.rmcs(parts_G[0], parts_H[0])

                if size > best_mcs[1]:
                    best_mcs = (mcs, size)

        return best_mcs

    def find_root_candidates(self, graph):
        candidates = set()

        # Add all vertices as candidates (simplified approach)
        for v in graph.vertices():
            candidates.add(v)

        return list(candidates)

    def decompose_graph(self, graph, root):
        """Decompose graph into parts with caching"""
        cache_key = (id(graph.grafo), root)
        if cache_key in self.decomposition_cache:
            return self.decomposition_cache[cache_key]

        parts = []
        is_basic = self.is_basic_root(graph, root)
        root_bps = BlockPreservingSubgraph(graph, root, is_basic, graph)
        parts.append(root_bps)

        # Get elementary parts
        elementary_parts = root_bps.get_elementary_parts()
        parts.extend(elementary_parts)

        self.decomposition_cache[cache_key] = parts
        return parts

    def is_basic_root(self, graph, root):
        """Check if root is basic-root (simplified version)"""
        # For outerplanar graphs, we consider a root basic if it has degree ≤ 2
        # or if it's in a single block
        if graph.grau(root) <= 1:
            return True

        blocks = graph.encontrar_blocos()
        incident_blocks = [block for block in blocks if root in block]
        return len(incident_blocks) == 1

    @lru_cache(maxsize=None)
    def rmcs(self, P, Q):
        """Rooted Maximum Common Subgraph (Algorithm 1)"""
        # Check root labels
        if P.graph.get_rotulo_vertice(P.root) != Q.graph.get_rotulo_vertice(Q.root):
            base = Grafo()
            base.adicionar_vertice(P.root, P.graph.get_rotulo_vertice(P.root))
            return (base, 1)

        # Handle different part types
        if P.type == "BPS" and Q.type == "BPS":
            if not P.is_basic_root and not Q.is_basic_root:
                return self.rmcs_compound(P, Q)
            else:
                return self.rmcs_bps(P, Q)
        elif P.type == "BSS" and Q.type == "BSS":
            return self.rmcs_bss(P, Q)
        else:
            base = Grafo()
            base.adicionar_vertice(P.root, P.graph.get_rotulo_vertice(P.root))
            return (base, 1)

    def rmcs_compound(self, P, Q):
        ep_P = P.get_elementary_parts()
        ep_Q = Q.get_elementary_parts()

        weights = np.zeros((len(ep_P), len(ep_Q)))
        matches = {}

        for i, part_P in enumerate(ep_P):
            for j, part_Q in enumerate(ep_Q):
                if self.are_compatible_parts(part_P, part_Q):
                    part_mcs, part_size = self.rmcs(part_P, part_Q)
                    weights[i, j] = part_size
                    matches[(i, j)] = (part_P, part_Q, part_mcs, part_size)

        row_ind, col_ind = linear_sum_assignment(weights, maximize=True)

        mcs_graph = Grafo()
        mcs_graph.adicionar_vertice(P.root, P.graph.get_rotulo_vertice(P.root))
        total_size = 1

        for i, j in zip(row_ind, col_ind):
            if weights[i, j] > 0:
                part_P, part_Q, part_mcs, part_size = matches[(i, j)]

                # Add vertices and edges from part_mcs
                for v in part_mcs.vertices():
                    if not mcs_graph.existe_vertice(v):
                        mcs_graph.adicionar_vertice(v, part_mcs.get_rotulo_vertice(v))
                        total_size += 1

                for u, v in part_mcs.arestas():
                    if not mcs_graph.grafo.has_edge(u, v):
                        mcs_graph.adicionar_aresta(u, v, part_mcs.get_rotulo_aresta(u, v))
                        total_size += 1

        return (mcs_graph, total_size)

    def are_compatible_parts(self, part1, part2):
        """Check if two parts are compatible (same type and root labels)"""
        if part1.type != part2.type:
            return False

        if part1.graph.get_rotulo_vertice(part1.root) != part2.graph.get_rotulo_vertice(part2.root):
            return False

        return True

    def rmcs_bps(self, P, Q):
        # Simple cases first
        if len(P.graph.vertices()) == 1 and len(Q.graph.vertices()) == 1:
            base = Grafo()
            base.adicionar_vertice(P.root, P.graph.get_rotulo_vertice(P.root))
            return (base, 1)

        # Try to match incident edges
        best_mcs = Grafo()
        best_mcs.adicionar_vertice(P.root, P.graph.get_rotulo_vertice(P.root))
        best_size = 1

        P_edges = [e for e in P.graph.arestas() if P.root in e]
        Q_edges = [e for e in Q.graph.arestas() if Q.root in e]

        for p_edge in P_edges:
            p_other = p_edge[0] if p_edge[1] == P.root else p_edge[1]
            for q_edge in Q_edges:
                q_other = q_edge[0] if q_edge[1] == Q.root else q_edge[1]

                if (P.graph.get_rotulo_vertice(p_other) == Q.graph.get_rotulo_vertice(q_other) and
                        P.graph.get_rotulo_aresta(*p_edge) == Q.graph.get_rotulo_aresta(*q_edge)):

                    # Create subgraphs without the matched edge
                    P_sub = P.graph.copy()
                    P_sub.grafo.remove_edge(*p_edge)

                    Q_sub = Q.graph.copy()
                    Q_sub.grafo.remove_edge(*q_edge)

                    # Compute MCS for subgraphs
                    P_sub_bps = BlockPreservingSubgraph(P_sub, p_other, False, P.original_graph)
                    Q_sub_bps = BlockPreservingSubgraph(Q_sub, q_other, False, Q.original_graph)
                    sub_mcs, sub_size = self.rmcs(P_sub_bps, Q_sub_bps)

                    # Build candidate MCS
                    candidate = Grafo()
                    candidate.adicionar_vertice(P.root, P.graph.get_rotulo_vertice(P.root))
                    candidate.adicionar_vertice(p_other, P.graph.get_rotulo_vertice(p_other))
                    candidate.adicionar_aresta(P.root, p_other, P.graph.get_rotulo_aresta(*p_edge))

                    for v in sub_mcs.vertices():
                        if not candidate.existe_vertice(v):
                            candidate.adicionar_vertice(v, sub_mcs.get_rotulo_vertice(v))
                    for u, v in sub_mcs.arestas():
                        if not candidate.grafo.has_edge(u, v):
                            candidate.adicionar_aresta(u, v, sub_mcs.get_rotulo_aresta(u, v))

                    candidate_size = len(candidate.vertices()) + len(candidate.arestas())

                    if candidate_size > best_size:
                        best_mcs = candidate
                        best_size = candidate_size

        return (best_mcs, best_size)

    def rmcs_bss(self, P, Q):
        """Algorithm 4: RMCS for basic-root BSS graphs"""
        base = Grafo()
        base.adicionar_vertice(P.root, P.graph.get_rotulo_vertice(P.root))
        base_size = 1

        # Simple case: single edge
        if P.split_edge and Q.split_edge:
            u, v = P.split_edge
            x, y = Q.split_edge

            if (P.graph.get_rotulo_vertice(u) == Q.graph.get_rotulo_vertice(x) and
                    P.graph.get_rotulo_vertice(v) == Q.graph.get_rotulo_vertice(y) and
                    P.graph.get_rotulo_aresta(u, v) == Q.graph.get_rotulo_aresta(x, y)):
                base.adicionar_vertice(v if u == P.root else u, P.graph.get_rotulo_vertice(v if u == P.root else u))
                base.adicionar_aresta(u, v, P.graph.get_rotulo_aresta(u, v))
                return (base, 3)

        return (base, base_size)


def encontrar_ciclo_hamiltoniano_outerplanar(grafo):
    """
    Encontra o ciclo Hamiltoniano em um grafo outerplanar biconectado.
    Retorna uma lista de arestas do ciclo ou None se não for possível encontrar.
    """
    G = grafo.grafo
    n = len(G.nodes)

    if n == 0:
        return []
    if n == 1:
        return []
    if n == 2:
        u, v = list(G.nodes)
        if G.has_edge(u, v):
            return [(u, v)]
        return []

    if all(deg == 2 for _, deg in G.degree()):
        try:
            cycle = nx.find_cycle(G)
            vertices = [cycle[0][0]]
            for i in range(len(cycle)):
                vertices.append(cycle[i][1])
            if vertices[0] == vertices[-1]:
                vertices.pop()
            return list(zip(vertices, vertices[1:] + [vertices[0]]))
        except nx.NetworkXNoCycle:
            return None

    try:
        pos = nx.planar_layout(G)

        start = min(pos, key=lambda v: pos[v][1])

        neighbors = list(G.neighbors(start))
        if len(neighbors) < 2:
            return None

        # Calcula os ângulos dos vizinhos em relação ao start
        def calc_angle(v):
            dx = pos[v][0] - pos[start][0]
            dy = pos[v][1] - pos[start][1]
            return np.arctan2(dy, dx)

        neighbors.sort(key=calc_angle)

        # Inicia a DFS a partir do start, sempre escolhendo o vizinho mais à direita
        visited = set()
        cycle = []
        stack = [(start, neighbors[0])]

        while stack:
            u, v = stack.pop()
            if v in visited:
                continue
            visited.add(v)
            cycle.append((u, v))

            # Obtém os vizinhos de v ordenados por ângulo
            v_neighbors = list(G.neighbors(v))
            v_neighbors.remove(u)  # Remove o vértice anterior
            if not v_neighbors:
                break

            # Calcula ângulos em relação a v
            def calc_angle_from_v(w):
                dx = pos[w][0] - pos[v][0]
                dy = pos[w][1] - pos[v][1]
                return np.arctan2(dy, dx)

            v_neighbors.sort(key=calc_angle_from_v)
            # Escolhe o vizinho com menor ângulo (mais à direita)
            next_v = v_neighbors[0]
            stack.append((v, next_v))

        # Verifica se o ciclo encontrado é Hamiltoniano
        if len(cycle) == n and len(set(cycle)) == n:
            return cycle
        else:
            # Fallback: usa o ciclo externo do embedding
            try:
                faces = nx.planar_faces(G, pos)
                outer_face = max(faces, key=len)
                if len(outer_face) == n:
                    edges = []
                    for i in range(len(outer_face)):
                        edges.append((outer_face[i], outer_face[(i + 1) % len(outer_face)]))
                    return edges
            except:
                pass
    except:
        pass

    # Fallback para grafos pequenos: enumeração limitada
    if n <= 10:
        try:
            for cycle in nx.simple_cycles(G):
                if len(cycle) == n:
                    edges = []
                    for i in range(len(cycle)):
                        edges.append((cycle[i], cycle[(i + 1) % len(cycle)]))
                    return edges
        except:
            pass

    return None

def calcular_mcs_bbp(grafo_G, grafo_H):
    solver = OuterplanarMCS()
    return solver.compute_mcs(grafo_G, grafo_H)

# Test cases
if __name__ == "__main__":
    # Test 1: Identical graphs
    print("Test 1: Identical graphs")
    G1 = Grafo()
    G1.adicionar_vertice(1, "C")
    G1.adicionar_vertice(2, "C")
    G1.adicionar_vertice(3, "O")
    G1.adicionar_aresta(1, 2, "single")
    G1.adicionar_aresta(2, 3, "double")
    G1.adicionar_aresta(1, 3, "single")

    G2 = G1.copy()

    mcs, size = calcular_mcs_bbp(G1, G2)
    expected_size = len(G1.vertices()) + len(G1.arestas())
    print(f"MCS Size: {size} (expected: {expected_size})")
    print("Vertices:", mcs.vertices())
    print("Edges:", mcs.arestas())
    print()

    # Test 2: Common substructure
    print("Test 2: Common substructure")
    G3 = Grafo()
    G3.adicionar_vertice(1, "C")
    G3.adicionar_vertice(2, "O")
    G3.adicionar_vertice(3, "C")
    G3.adicionar_aresta(1, 2, "double")
    G3.adicionar_aresta(2, 3, "single")

    G4 = Grafo()
    G4.adicionar_vertice(4, "C")
    G4.adicionar_vertice(5, "O")
    G4.adicionar_vertice(6, "N")
    G4.adicionar_aresta(4, 5, "double")
    G4.adicionar_aresta(5, 6, "single")

    mcs, size = calcular_mcs_bbp(G3, G4)
    print(f"MCS Size: {size} (expected: 3)")
    print("Vertices:", mcs.vertices())
    print("Edges:", mcs.arestas())
    print()

    # Test 3: Single vertex
    print("Test 3: Single vertex")
    G5 = Grafo()
    G5.adicionar_vertice(1, "C")

    G6 = Grafo()
    G6.adicionar_vertice(2, "C")

    mcs, size = calcular_mcs_bbp(G5, G6)
    print(f"MCS Size: {size} (expected: 1)")
    print("Vertices:", mcs.vertices())
    print()

    # Test 4: Single edge
    print("Test 4: Single edge")
    G7 = Grafo()
    G7.adicionar_vertice(1, "C")
    G7.adicionar_vertice(2, "O")
    G7.adicionar_aresta(1, 2, "double")

    G8 = Grafo()
    G8.adicionar_vertice(3, "C")
    G8.adicionar_vertice(4, "O")
    G8.adicionar_aresta(3, 4, "double")

    mcs, size = calcular_mcs_bbp(G7, G8)
    print(f"MCS Size: {size} (expected: 3)")
    print("Vertices:", mcs.vertices())
    print("Edges:", mcs.arestas())
    print()

    # Test 5: Common cycle - CORRIGIDO
    print("Test 5: Common cycle")
    G9 = Grafo()
    G9.adicionar_vertice(1, "C")
    G9.adicionar_vertice(2, "C")
    G9.adicionar_vertice(3, "C")
    G9.adicionar_aresta(1, 2, "single")
    G9.adicionar_aresta(2, 3, "single")
    G9.adicionar_aresta(3, 1, "single")  # Triângulo

    G10 = Grafo()
    G10.adicionar_vertice(4, "C")
    G10.adicionar_vertice(5, "C")
    G10.adicionar_vertice(6, "C")
    G10.adicionar_aresta(4, 5, "single")
    G10.adicionar_aresta(5, 6, "single")  # Cadeia linear

    mcs, size = calcular_mcs_bbp(G9, G10)
    # CORREÇÃO: O MCS entre um triângulo e uma cadeia linear deve ser a cadeia de 2 arestas (tamanho 5)
    print(f"MCS Size: {size} (expected: 5)")
    print("Vertices:", mcs.vertices())
    print("Edges:", mcs.arestas())
    print()

    # Test 6: Multiple blocks
    print("Test 6: Multiple blocks")
    G11 = Grafo()
    G11.adicionar_vertice(1, "C")
    G11.adicionar_vertice(2, "C")
    G11.adicionar_vertice(3, "C")
    G11.adicionar_vertice(4, "O")
    G11.adicionar_aresta(1, 2, "single")
    G11.adicionar_aresta(2, 3, "single")
    G11.adicionar_aresta(3, 1, "single")  # Triângulo
    G11.adicionar_aresta(3, 4, "single")  # Ponte para O

    G12 = Grafo()
    G12.adicionar_vertice(5, "C")
    G12.adicionar_vertice(6, "C")
    G12.adicionar_vertice(7, "O")
    G12.adicionar_aresta(5, 6, "single")
    G12.adicionar_aresta(6, 7, "single")  # Cadeia linear C-C-O

    mcs, size = calcular_mcs_bbp(G11, G12)
    print(f"MCS Size: {size} (expected: 5)")
    print("Vertices:", mcs.vertices())
    print("Edges:", mcs.arestas())
    print()

    # Test 7: Size correction verification
    print("Test 7: Size correction verification")
    G_corr = Grafo()
    G_corr.adicionar_vertice(1, "C")
    G_corr.adicionar_vertice(2, "C")
    G_corr.adicionar_aresta(1, 2, "single")

    mcs, size = calcular_mcs_bbp(G_corr, G_corr)
    print(f"MCS Size: {size} (expected: 3)")
    print("Vertices:", mcs.vertices())
    print("Edges:", mcs.arestas())
    print()

    # Test 8: Detailed size verification
    print("Test 8: Detailed size verification")
    G_test = Grafo()
    G_test.adicionar_vertice(1, "C")
    G_test.adicionar_vertice(2, "C")
    G_test.adicionar_vertice(3, "O")
    G_test.adicionar_aresta(1, 2, "single")
    G_test.adicionar_aresta(2, 3, "double")

    # Deveria ter tamanho 5 (3 vértices + 2 arestas)
    print(f"Tamanho esperado do G_test: {len(G_test.vertices()) + len(G_test.arestas())}")

    mcs, size = calcular_mcs_bbp(G_test, G_test)
    print(f"MCS Size: {size} (expected: 5)")
    print("Vertices:", mcs.vertices())
    print("Edges:", mcs.arestas())