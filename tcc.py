from structures.listaAdj import ListaAdjG
from collections import deque
import networkx as nx
import random


def OptimalTopDownCommonSubtree(T1, v, T2, u, Mp, memo=None):
    if memo is None:
        memo = {}

    if (v, u) in memo:
        return memo[(v, u)]

    if not T1.G.has_node(v) or not T2.G.has_node(u):
        raise ValueError(f"Vertex {v} or {u} not found in respective trees")

    if T1.is_leaf(v) or T2.is_leaf(u):
        result = (max(T1.altura().get(v, 0), T2.altura().get(u, 0)), Mp)
        memo[(v, u)] = result
        return result

    Gvu = ListaAdjG(orientado=False)
    Gvu.DefinirN(T1.n + T2.n + 2)

    for i in T1.G.neighbors(v):
        for j in T2.G.neighbors(u):
            Gvu.AdicionarAresta(i, j)

    dum_v = T1.n + 1
    dum_u = T1.n + 2

    for x in T1.G.neighbors(v):
        Gvu.AdicionarAresta(x, dum_u)
    for y in T2.G.neighbors(u):
        Gvu.AdicionarAresta(dum_v, y)

    heightsT1 = T1.altura()
    heightsT2 = T2.altura()

    edge_weights = {}

    for u_vertex in Gvu.G.nodes:
        for v_vertex in Gvu.G.nodes:
            if u_vertex == dum_v:
                edge_weights[(u_vertex, v_vertex)] = heightsT2.get(v_vertex, 0) + 1
            elif v_vertex == dum_u:
                edge_weights[(u_vertex, v_vertex)] = heightsT1.get(u_vertex, 0) + 1
            else:
                if u_vertex in T1.G.nodes and v_vertex in T2.G.nodes:
                    sub_distance, _ = OptimalTopDownCommonSubtree(T1, u_vertex, T2, v_vertex, Mp, memo)
                    edge_weights[(u_vertex, v_vertex)] = sub_distance
                else:
                    edge_weights[(u_vertex, v_vertex)] = 0

    Gvu.G.graph['edge_weights'] = edge_weights

    Mvu = solveOptimalPerfectMatching(Gvu)
    distance = max(edge_weights.get((u, v), 0) for u, v in Mvu)
    Mvu = {e for e in Mvu if not is_dummy_vertex(e[0]) and not is_dummy_vertex(e[1])}
    Mp.update(Mvu)

    result = (distance, Mp)
    memo[(v, u)] = result
    return result

def solveOptimalPerfectMatching(Gvu):
    if not isinstance(Gvu, ListaAdjG):
        raise TypeError("Gvu must be an instance of ListaAdjG")

    # Obter todos os pesos das arestas e ordená-los
    edge_weights = sorted(set(Gvu.getPeso(u, v) for u, v in Gvu.G.edges))
    print("Pesos das arestas:", edge_weights)

    left, right = 0, len(edge_weights) - 1
    best_matching = set()
    best_weight = float('inf')

    # Algoritmo de busca binária sobre os pesos
    while left <= right:
        mid = (left + right) // 2
        weight_threshold = edge_weights[mid]
        print(f"Threshold de peso: {weight_threshold}")

        # Criar subgrafo com as arestas válidas para o peso atual
        subgraph = nx.Graph()
        subgraph.add_nodes_from(range(1, Gvu.n + 1))

        for u, v in Gvu.G.edges:
            if Gvu.getPeso(u, v) <= weight_threshold:
                subgraph.add_edge(u, v)

        # Definindo corretamente os conjuntos esquerdo e direito
        left_set = {u for u in range(1, Gvu.n + 1) if u % 2 == 1}
        right_set = {u for u in range(1, Gvu.n + 1) if u % 2 == 0}

        print(f"Conjunto esquerdo: {left_set}")
        print(f"Conjunto direito: {right_set}")

        # Filtrando as arestas que conectam apenas os conjuntos esquerdo e direito
        valid_edges = [(u, v) for u, v in subgraph.edges if u in left_set and v in right_set]
        print(f"Arestas válidas entre os conjuntos: {valid_edges}")

        subgraph.clear_edges()
        subgraph.add_edges_from(valid_edges)

        # Encontrar o emparelhamento máximo bipartido
        matching = nx.bipartite.maximum_matching(subgraph, top_nodes=left_set)
        print(f"Emparelhamento encontrado: {matching}")

        # Ajustando o formato do emparelhamento (bipartite.maximum_matching retorna um dicionário)
        matching_set = set((u, v) for u, v in matching.items() if u in left_set)

        # Se o emparelhamento for perfeito (tamanho dos conjuntos esquerdo e direito coincidirem)
        print(f"Emparelhamento final: {matching_set}")
        print(f"Tamanho do emparelhamento: {len(matching_set)}")
        print(f"Tamanho do conjunto esquerdo: {len(left_set)}")

        if len(matching_set) == len(left_set) and all(v in right_set for u, v in matching_set):
            best_matching = matching_set
            best_weight = weight_threshold
            right = mid - 1
        else:
            left = mid + 1

    print("Emparelhamento ótimo perfeito:", best_matching)
    return best_matching if best_matching and len(best_matching) == len(left_set) else set()

def hopcroft_karp(graph, left_set, right_set): #funciona
    pair_u = {u: None for u in left_set}
    pair_v = {v: None for v in right_set}
    dist = {}

    def bfs():
        queue = deque()
        for u in left_set:
            if pair_u[u] is None:
                dist[u] = 0
                queue.append(u)
            else:
                dist[u] = float('inf')
        dist[None] = float('inf')

        while queue:
            u = queue.popleft()
            if dist[u] < dist[None]:
                for v in graph.N(u):
                    if v in pair_v:
                        if pair_v[v] is None:
                            dist[None] = dist[u] + 1
                        elif dist.get(pair_v[v], float('inf')) == float('inf'):
                            dist[pair_v[v]] = dist[u] + 1
                            queue.append(pair_v[v])
        return dist[None] != float('inf')

    def dfs(u):
        if u is not None:
            for v in graph.N(u):
                if v in pair_v:
                    if dist.get(pair_v[v], float('inf')) == dist[u] + 1:
                        if pair_v[v] is None or dfs(pair_v[v]):
                            pair_v[v] = u
                            pair_u[u] = v
                            return True
            dist[u] = float('inf')
            return False
        return True

    matching = set()
    while bfs():
        for u in left_set:
            if pair_u[u] is None and dfs(u):
                matching.add((u, pair_u[u]))

    return matching

# Criar um novo grafo
G = ListaAdjG()
G.DefinirN(6)

# Adicionar arestas com pesos
G.AdicionarAresta(1, 2, 5)
G.AdicionarAresta(1, 3, 2)
G.AdicionarAresta(1, 4, 6)
G.AdicionarAresta(2, 4, 1)
G.AdicionarAresta(2, 5, 3)
G.AdicionarAresta(3, 5, 4)
G.AdicionarAresta(3, 6, 7)
G.AdicionarAresta(4, 6, 8)

# Chamar o método solveOptimalPerfectMatching
resultado = solveOptimalPerfectMatching(G)

# Exibir o resultado
print("Resultado do emparelhamento ótimo:", resultado)