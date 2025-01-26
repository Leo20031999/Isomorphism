from structures.listaAdj import Grafo
from collections import deque
import networkx as nx

MAX_ITERACOES = 1000
MAX_DEPTH = 50

def OptimalTopDownCommonSubtree(T1, v, T2, u, Mp, memo=None, iteracoes=0, depth=0):
    if memo is None:
        memo = {}

    pilha = [(v, u, None, "início")]
    visited_T1 = set()
    visited_T2 = set()

    while pilha:
        v, u, estado_parcial, op = pilha.pop()

        if (v, u) in memo:
            if estado_parcial:
                estado_parcial["distância"] = memo[(v, u)][0]
                estado_parcial["Mp"] = memo[(v, u)][1]
            continue

        if op == "início":
            if T1.is_leaf(v) or T2.is_leaf(u):
                altura_v = T1.altura().get(v, 0)
                altura_u = T2.altura().get(u, 0)
                result = (max(altura_v, altura_u), Mp)
                memo[(v, u)] = result
                if estado_parcial:
                    estado_parcial["distância"] = result[0]
                    estado_parcial["Mp"] = result[1]
                continue

            estado_parcial = {
                "distância": 0,
                "Mp": Mp.copy(),
                "alturas": (T1.altura(), T2.altura()),
                "pares_subárvores": [],
                "grafo_aux": None,
            }

            pilha.append((v, u, estado_parcial, "resolvido"))

            Gvu = Grafo()
            for i in range(T1.num_nos + T2.num_nos + 2):
                Gvu.grafo.add_node(i)

            vizinhos_v = T1.vizinhanca(v)
            vizinhos_u = T2.vizinhanca(u)

            for i in vizinhos_v:
                for j in vizinhos_u:
                    if not Gvu.grafo.has_edge(i, j):
                        Gvu.grafo.add_edge(i, j)
                        iteracoes += 1
                    if iteracoes >= MAX_ITERACOES:
                        return memo.get((v, u), ("Limite de iterações atingido", Mp))

            dum_v, dum_u = T1.num_nos + 1, T1.num_nos + 2
            for x in vizinhos_v:
                Gvu.grafo.add_edge(x, dum_u)
            for y in vizinhos_u:
                Gvu.grafo.add_edge(dum_v, y)

            estado_parcial["grafo_aux"] = Gvu

            for i in vizinhos_v:
                for j in vizinhos_u:
                    if (i, j) not in visited_T1 and (i, j) not in visited_T2:
                        pilha.append((i, j, estado_parcial, "início"))
                        estado_parcial["pares_subárvores"].append((i, j))
                        visited_T1.add((i, j))
                        visited_T2.add((i, j))

        elif op == "resolvido":
            Gvu = estado_parcial["grafo_aux"]
            dum_v, dum_u = T1.num_nos + 1, T1.num_nos + 2
            alturasT1, alturasT2 = estado_parcial["alturas"]

            for u_vertex, v_vertex in Gvu.arestas():
                if u_vertex == dum_v:
                    Gvu.set_peso(u_vertex, v_vertex, alturasT2.get(v_vertex, 0) + 1)
                elif v_vertex == dum_u:
                    Gvu.set_peso(u_vertex, v_vertex, alturasT1.get(u_vertex, 0) + 1)
                elif (u_vertex, v_vertex) in memo:
                    Gvu.set_peso(u_vertex, v_vertex, memo[(u_vertex, v_vertex)][0])
                else:
                    Gvu.set_peso(u_vertex, v_vertex, 0)

            Mvu = solve_optimal_perfect_matching(Gvu)
            distância = max(Gvu.get_weight(u, v) for u, v in Mvu)

            Mvu_filtered = {(i, j) for i, j in Mvu if not is_dummy_vertex(i, dum_v, dum_u) and not is_dummy_vertex(j, dum_v, dum_u)}

            Mp.update(Mvu_filtered)

            result = (distância, Mp)
            memo[(v, u)] = result
            if estado_parcial:
                estado_parcial["distância"] = result[0]
                estado_parcial["Mp"] = result[1]

    return memo[(v, u)], iteracoes

def solve_optimal_perfect_matching(gvu):
    if not isinstance(gvu, Grafo):
        raise TypeError("O parâmetro 'gvu' deve ser uma instância de Grafo.")
    
    edge_weights = sorted(set(gvu.get_peso(u, v) for u, v in gvu.arestas()))
    
    if not edge_weights:
        return set() 

    def is_perfect_matching(weight_threshold):
        valid_edges = [(u, v) for u, v in gvu.arestas() if gvu.get_peso(u, v) <= weight_threshold]
        if not valid_edges:
            return False, set()
        subgraph = nx.Graph()
        subgraph.add_edges_from(valid_edges)

        left_set = {u for u in gvu.vertices() if u % 2 == 1}
        right_set = {u for u in gvu.vertices() if u % 2 == 0}

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

    if best_matching and len(best_matching) == len({u for u in range(1, gvu.n + 1) if u % 2 == 1}):
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

def teste_OptimalTopDownCommonSubtree_v2():
    T1 = Grafo()
    T2 = Grafo()

    T1.definir_n(4)
    T1.adicionar_aresta(1, 2)
    T1.adicionar_aresta(2, 3)
    T1.adicionar_aresta(2, 4)

    T2.definir_n(4)
    T2.adicionar_aresta(1, 2)
    T2.adicionar_aresta(2, 3)
    T2.adicionar_aresta(2, 4)

    Mp = set()
    print("Estado inicial de T1 e T2:")
    print("Vértices de T1:", list(T1.vertices()))
    print("Vértices de T2:", list(T2.vertices()))
    print("Arestas de T1:", T1.arestas())
    print("Arestas de T2:", T2.arestas())

    for v in T1.vertices():
        print(f"Vizinhança do vértice {v} em T1:", T1.vizinhanca(v))
    for v in T2.vertices():
        print(f"Vizinhança do vértice {v} em T2:", T2.vizinhanca(v))

    print("\nIniciando execução do OptimalTopDownCommonSubtree...")
    try:
        resultado, iteracoes = OptimalTopDownCommonSubtree(T1, 1, T2, 1, Mp)
        print("Resultado da execução:", resultado)
        print("Número de iterações: ", iteracoes)
    except Exception as e:
        print("Erro durante a execução do algoritmo:", e)
        raise

if __name__ == "__main__":
    print("Iniciando o teste do método OptimalTopDownCommonSubtree...")
    teste_OptimalTopDownCommonSubtree_v2()
    print("Teste finalizado.")

