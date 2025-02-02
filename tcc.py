from structures.listaAdj import Grafo
from collections import deque
import networkx as nx

MAX_ITERACOES = 1000
MAX_DEPTH = 50

def OptimalTopDownCommonSubtree(T1, v, T2, u, Mp, memo=None, iteracoes=0, depth=0):
    if memo is None:
        memo = {}

    pilha = [(v, u, None, "início")]
    visited = set()

    while pilha:
        v, u, estado_parcial, op = pilha.pop()

        if (v, u) in memo:
            if estado_parcial:
                estado_parcial["distância"] = memo[(v, u)][0]
                estado_parcial["Mp"] = memo[(v, u)][1]
            continue

        if op == "início":
            print(f"Verificando se {u} (de T2) é folha...")
            if u not in T2.vertices():
                raise ValueError(f"Erro: Tentativa de acessar vértice inexistente {u} em T2.")

            if T1.is_leaf(v) or T2.is_leaf(u):
                altura_v = T1.altura().get(v, 0)
                altura_u = T2.altura().get(u, 0)
                memo[(v, u)] = (max(altura_v, altura_u), Mp)
                continue

            estado_parcial = {"distância": 0, "Mp": Mp.copy(), "pares_subárvores": []}
            pilha.append((v, u, estado_parcial, "resolvido"))

            offset = max(T1.vertices(), default=0) 
            vizinhos_u_orig = T2.vizinhanca(u)

            vizinhos_u = [j for j in vizinhos_u_orig if j in T2.vertices()]

            print(f"Offset utilizado: {offset}")
            print(f"Vizinhos de {u} em T2 antes do deslocamento: {vizinhos_u_orig}")
            print(f"Vizinhos de {u} em T2 após deslocamento: {vizinhos_u}")
            print(f"Todos os vértices de T2: {T2.vertices()}")

            Gvu = Grafo()
            total_nodes = T1.num_nos + T2.num_nos
            
            for i in range(1, total_nodes + 3):
                Gvu.grafo.add_node(i)

            vizinhos_v = T1.vizinhanca(v)
            vizinhos_u = [j + offset for j in vizinhos_u if (j + offset) in T2.vertices()]

            for i in vizinhos_v:
                for j in vizinhos_u:
                    Gvu.grafo.add_edge(i, j, weight=1)
                    iteracoes += 1
                    if iteracoes >= MAX_ITERACOES:
                        return memo.get((v, u), ("Limite de iterações atingido", Mp))

            dum_v, dum_u = total_nodes + 1, total_nodes + 2
            for x in vizinhos_v:
                Gvu.grafo.add_edge(x, dum_u, weight=1)
            for y in vizinhos_u:
                Gvu.grafo.add_edge(dum_v, y, weight=1)

            estado_parcial["grafo_aux"] = Gvu
            
            for i in vizinhos_v:
                for j in vizinhos_u:
                    if (i, j) not in visited:
                        pilha.append((i, j, estado_parcial, "início"))
                        estado_parcial["pares_subárvores"].append((i, j))
                        visited.add((i, j))

        elif op == "resolvido":
            Gvu = estado_parcial["grafo_aux"]
            dum_v, dum_u = total_nodes + 1, total_nodes + 2

            for u_vertex, v_vertex in Gvu.arestas():
                if u_vertex == dum_v:
                    Gvu.set_peso(u_vertex, v_vertex, T2.altura().get(v_vertex, 0) + 1)
                elif v_vertex == dum_u:
                    Gvu.set_peso(u_vertex, v_vertex, T1.altura().get(u_vertex, 0) + 1)
                elif (u_vertex, v_vertex) in memo:
                    Gvu.set_peso(u_vertex, v_vertex, memo[(u_vertex, v_vertex)][0])
                else:
                    Gvu.set_peso(u_vertex, v_vertex, 0)

            print(f"Arestas do Gvu: {list(Gvu.arestas())}")
            print(f"Pesos das arestas: {[Gvu.get_peso(u, v) for u, v in Gvu.arestas()]}")
            Mvu = solve_optimal_perfect_matching(Gvu)
            distância = max((Gvu.get_weight(u, v) for u, v in Mvu), default=0)
            Mvu_filtered = {(i, j) for i, j in Mvu if i <= T1.num_nos and j > T1.num_nos}

            Mp.update(Mvu_filtered)
            memo[(v, u)] = (distância, Mp)

    return memo[(v, u)], iteracoes

def solve_optimal_perfect_matching(gvu):
    if not isinstance(gvu, Grafo):
        raise TypeError("O parâmetro 'gvu' deve ser uma instância de Grafo.")
    
    todos_vertices = sorted(gvu.vertices())
    if not todos_vertices:
        return set()
    
    num_nos_T1 = max(v for v in todos_vertices if v <= len(todos_vertices) // 2)
    num_nos_T2 = len(todos_vertices) - num_nos_T1
    print(f"num_nos_T1: {num_nos_T1}, num_nos_T2: {num_nos_T2}")
    edge_weights = sorted(set(gvu.get_peso(u, v) for u, v in gvu.arestas()))
    if not edge_weights:
        return set() 

    def is_perfect_matching(weight_threshold):
        valid_edges = [(u, v) for u, v in gvu.arestas() if gvu.get_peso(u, v) <= weight_threshold]
        subgraph = nx.Graph()
        subgraph.add_edges_from(valid_edges)

        left_set = {u for u in gvu.vertices() if u <= num_nos_T1}  
        right_set = {u for u in gvu.vertices() if u > num_nos_T1 and u <= num_nos_T1 + num_nos_T2}

        print(f"left_set: {left_set}")
        print(f"right_set: {right_set}")

        subgraph.add_nodes_from(left_set)
        subgraph.add_nodes_from(right_set)

        for u in left_set:
            if subgraph.degree(u) == 0:
                return False, set()

        matching = nx.bipartite.maximum_matching(subgraph, top_nodes=left_set)
        matching_set = {(u, v) for u, v in matching.items() if u in left_set and v in right_set}

        is_perfect = len(matching_set) == len(left_set) and all(v in right_set for _, v in matching_set)
        return is_perfect, matching_set


    left, right = 0, len(edge_weights) - 1
    best_matching = set()
    best_weight = float('inf')
    while left <= right:
        mid = (left + right) // 2
        weight_threshold = edge_weights[mid]
        is_perfect, matching_set = is_perfect_matching(weight_threshold)
        if is_perfect:
            best_matching = matching_set
            best_weight = weight_threshold
            right = mid - 1
        else:
            left = mid + 1
    print(f"Melhor emparelhamento encontrado: {best_matching}")
    if best_matching and len(best_matching) == len({u for u in gvu.vertices() if u <= num_nos_T1}):
        return best_matching
    return set()

def hopcroft_karp(graph, left_set, right_set):
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
                if v in pair_v and dist.get(pair_v[v], float('inf')) == dist[u] + 1:
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

def teste_OptimalTopDownCommonSubtree_v4():
    T1 = Grafo()
    T2 = Grafo()

    arestas = [(1, 2), (2, 3), (2, 4), (4, 5)]
    T1.definir_n(5)
    T2.definir_n(5)
    
    for u, v in arestas:
        T1.adicionar_aresta(u, v)
        T2.adicionar_aresta(u, v)
    
    print("Iniciando teste...")
    Mp = set()
    
    try:
        resultado, iteracoes = OptimalTopDownCommonSubtree(T1, 1, T2, 1, Mp)
        assert resultado[0] > 0, "O emparelhamento deve ser maior que 0"
        print(f"Emparelhamento bem-sucedido! Distância: {resultado[0]}, Iterações: {iteracoes}")
    except Exception as e:
        print(f"Erro: {e}")
        raise

if __name__ == "__main__":
    teste_OptimalTopDownCommonSubtree_v4()