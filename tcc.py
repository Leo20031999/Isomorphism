from structures.listaAdj import listaAdj
from collections import deque

def OptimalTopDownCommonSubtree(T1,v,T2,u,Mp, memo = None):
    if memo is None:
        memo = {}

    if (v,u) in memo:
        return memo[(v,u)]
    
    if v not in T1.adj_list or u not in T2.adj_list:
        raise ValueError(f"Vertex {v} or {u} not found in respective trees")
    
    if T1.is_leaf(v) or T2.is_leaf(u):
        result = (max(T1.heights().get(v,0), T2.heights().get(u,0)), Mp)
        memo[(v,u)] = result
        return result
    
    Gvu = listaAdj(orientado=False)
    Gvu.DefinirN(T1.n + T2.n + 2)

    for i in T1.adj_list.get(v,[]):
        for j in T2.adj_list.get(u, []):
            Gvu.AdicionarAresta(i,j)

    dum_v = T1.n + 1
    dum_u = T1.n + 2

    for x in T1.adj_list.get(v,[]):
        Gvu.AdicionarAresta(x,dum_u)
    for y in T2.adj_list.get(u,[]):
        Gvu.AdicionarAresta(dum_v, y)

    heightsT1 = T1.altura()
    heightsT2 = T2.altura()

    edge_weights = {}

    for u_vertex in Gvu.V():
        for v_vertex in Gvu.V():
            if u_vertex == dum_v:
                edge_weights[(u_vertex,v_vertex)] = heightsT2.get(v_vertex, 0) + 1
            elif v_vertex == dum_u:
                edge_weights[(u_vertex,v_vertex)] = heightsT1.get(u_vertex, 0) + 1
            else:
                if u_vertex in T1.adj_list and v_vertex in T2.adj_list:
                    sub_distance, _ = OptimalTopDownCommonSubtree(T1, u_vertex, T2, v_vertex, Mp, memo)
                    edge_weights[(u_vertex,v_vertex)] = sub_distance
                else:
                    edge_weights[(u_vertex,v_vertex)] = 0
    
    Gvu.adj_list = edge_weights

    Mvu = solveOptimalPerfectMatching(Gvu)
    distance = max(edge_weights.get((u,v), 0) for u,v in Mvu)
    Mvu = {e for e in Mvu if not (is_dummy_vertex(e[0]) or is_dummy_vertex(e[1]))}
    Mp.update(Mvu)

    result = (distance, Mp)
    memo[(v,u)] = result
    return result

def is_dummy_vertex(vertex):
    return isinstance(vertex, str) and vertex.startswith("dum")

def reconstructionOfMapping(T1, r1,r2, Mp, M): #FUNCIONA
    M.add((r1,r2))
    P1 = []
    def visit(v):
        P1.append(v)
    preorder(T1,r1,visit)
    for v in P1:
        for(vertex,u) in Mp:
            if vertex == v:
                u_parent = parent(T1,u)
                v_parent = parent(T1,v)
                if (u_parent, v_parent) in M:
                    M.add((v,u))
    return M

def parent(T, v):
    for u in T.N(v):
        if T.is_leaf(u):
            continue
        if u!=v:
            return u
    return None

def hopcroft_karp(graph, left_set, right_set): #FUNCIONA

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
                    if pair_v[v] is None:
                        if dist.get(pair_v[v], float('inf')) == float('inf'):
                            dist[pair_v[v]] = dist[u] + 1
                            queue.append(pair_v[v])
        return dist[None] != float('inf')
    
    def dfs(u):
        if u is not None:
            for v in graph.N(u):
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
                if pair_u[u] is not None:
                    matching.add((min(u, pair_u[u]), max(u, pair_u[u])))
    print("hopcroft_karp terminado")
    return matching


def solveOptimalPerfectMatching(Gvu):
    if not isinstance(Gvu,listaAdj):
        raise TypeError("Gvu must be an instance of listaAdj")

    edge_weights = sorted(set(Gvu.getPeso(u,v) for u,v in Gvu.E()))

    left, right = 0, len(edge_weights) - 1
    best_matching = set()
    best_weight = float('inf')

    while left<=right:
        mid = (left+right) // 2
        weight_threshold = edge_weights[mid]

        subgraph = listaAdj(orientado=False)
        subgraph.DefinirN(Gvu.n)
        for u,v in Gvu.E():
            if Gvu.getPeso(u,v) <= weight_threshold:
                subgraph.AdicionarAresta(u,v)

        print(f"Subgrafo com peso <= {weight_threshold}: {[ (u, v) for u, v in subgraph.E() ]}")
        left_set = set(subgraph.V())
        right_set = {v for u in left_set for v in subgraph.N(u) }

        for v in left_set:
            vizinhos = list(subgraph.N(v))
            if vizinhos:
                right_set.update(vizinhos)

        matching = hopcroft_karp(subgraph, left_set, right_set)
        print(f"Conjunto esquerdo: {left_set}")
        print(f"Conjunto direito: {right_set}")
        print(f"Emparelhamento encontrado com peso <= {weight_threshold}: {matching}")
        if len(matching) * 2 == len(left_set) and len(left_set) % 2 == 0:
            best_matching = matching
            best_weight = weight_threshold
            right = mid - 1
        else:
            left = mid + 1
    return best_matching

def calcularAlturas(T): #FUNCIONA
    alturas ={}
    def dfs(v, altura_atual, visitado):
        visitado[v] = True
        alturas[v] = altura_atual
        for vizinho in T.N(v):
            if not visitado[vizinho]:
                dfs(vizinho, altura_atual+1, visitado)
    visitado = [False] * (T.n + 1)
    dfs(1,0,visitado)
    return alturas

def hausdorffDistanceBetweenTrees(grafo1, grafo2):
    hd = float('inf')
    O = set()
    r1=grafo1.center()[0]
    r2 = None
    T1_heights = grafo1.heights()

    for u in grafo2.adj_list:
        Mp=set()
        T2_heights = grafo2.heights()
        distance = OptimalTopDownCommonSubtree(grafo1,r1,grafo2,u,Mp)
    
    if distance < hd:
        hd = distance
        r2 = u
        O = Mp

    if r2 is None:
        raise ValueError("R2 doesn't have any value")

    M = set()
    reconstructionOfMapping(grafo1,r1,r2,O,M)
    return hd, M

def preorder(graph,v,visit): #funciona
    pilha = [v]
    visitados = set()
    while pilha:
        atual = pilha.pop()
        if atual in visitados:
            continue
        visitados.add(atual)
        visit(atual)
        for w in list(graph.N(atual)):
            if w not in visitados:
                pilha.append(w)

def test_complex_perfect_matching():
    Gvu = listaAdj(orientado=False)
    Gvu.DefinirN(10)  # Definindo 10 vértices

    # Adicionando arestas com pesos
    Gvu.AdicionarAresta(1, 2, peso=1)
    Gvu.AdicionarAresta(2, 3, peso=1)
    Gvu.AdicionarAresta(3, 4, peso=2)
    Gvu.AdicionarAresta(4, 5, peso=2)
    Gvu.AdicionarAresta(5, 6, peso=1)
    Gvu.AdicionarAresta(6, 7, peso=1)
    Gvu.AdicionarAresta(7, 8, peso=2)
    Gvu.AdicionarAresta(8, 9, peso=2)
    Gvu.AdicionarAresta(9, 10, peso=1)
    Gvu.AdicionarAresta(1, 10, peso=3)  # Conexão extra, mas com peso alto

    resultado = solveOptimalPerfectMatching(Gvu)
    print(f"Emparelhamento perfeito ótimo: {resultado}")

# Executando o teste
test_complex_perfect_matching()