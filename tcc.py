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

    heightsT1 = T1.heights()
    heightsT2 = T2.heights()

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

    Mvu = solveOptimalPatternMatching(Gvu)
    distance = max(edge_weights.get((u,v), 0) for u,v in Mvu)
    Mvu = {e for e in Mvu if not (is_dummy_vertex(e[0]) or is_dummy_vertex(e[1]))}
    Mp.update(Mvu)

    result = (distance, Mp)
    memo[(v,u)] = result
    return result

def is_dummy_vertex(vertex):
    return isinstance(vertex, str) and vertex.startswith("dum")

def reconstructionOfMapping(T1, r1,r2, Mp, M):
    M.add((r1,r2))
    P1 = preorder(T1)
    for v in P1:
        for (vertex, u) in Mp:
            if vertex == v:
                if (parent(T1, v), parent(T1, u)) in M:
                    M.add(v,u)
    return M

#checar preorder e ver qual é a melhor assinatura para esse método
def preorder(T, vertex):
    traversal = [node]
    if node in T:
        for child in T[node]:
            traversal.extend(preorder(T,child))
    return traversal

def parent(T, vertex):
    for parent, children in T.items():
        if vertex in children:
            return parent
        return None

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
                for v in graph.adj_list.get(u,[]):
                    if dist[pair_v[v]] == float('inf'):
                        dist[pair_v[v]] = dist[u] + 1
                        queue.append(pair_v[v])
        return dist[None] != float('inf')
    
    def dfs(u):
        if u is None:
            for v in graph.adj_list.get(u,[]):
                if dist[pair_v[v]] == dist[u] + 1 and dfs(pair_v[v]):
                    pair_v[v] = u
                    pair_u[u] = v
                    return True
            dist[u] = float('inf')
            return True
        return False
    
    matching = set()
    while bfs():
        for u in left_set:
            if pair_u[u] is None and dfs(u):
                matching.add((u, pair_u[u]))
    return matching


def solveOptimalPatternMatching(Gvu):
    if not isinstance(Gvu,WeightedGraph):
        raise TypeError("Gvu must be an instance of WeightedGraph")

    edge_weights = sorted(set(Gvu.get_weight(u,v) for u,v in Gvu.get_edges()))

    left, right = 0, len(edge_weights) - 1
    best_matching = set()
    best_weight = float('inf')

    while left<=right:
        mid = (left+right) // 2
        weight_threshold = edge_weights[mid]

        subgraph = WeightedGraph()
        for u,v in Gvu.get_edges():
            if Gvu.get_weight(u,v) <= weight_threshold:
                subgraph.add_edge(u,v,Gvu.get_weight(u,v))
        left_set = set(u for u in subgraph.adj_list.keys())
        right_set = set(v for neighbors in subgraph.adj_list.values() for v in neighbors)

        matching = hopcroft_karp(subgraph, left_set, right_set)

        if len(matching) == len(left_set):
            best_matching = matching
            best_weight = weight_threshold
            right = mid - 1
        else:
            left = mid + 1
    
    return best_matching

def calcularAlturas(T):
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

T1 = listaAdj(orientado=False)
T2 = listaAdj(orientado=False) 


T1.DefinirN(5)
T2.DefinirN(4) 

T1.AdicionarAresta(1, 2)
T1.AdicionarAresta(1, 3)
T1.AdicionarAresta(2, 4)
T1.AdicionarAresta(2, 5)

T2.AdicionarAresta(1, 2)
T2.AdicionarAresta(1, 3)
T2.AdicionarAresta(3, 4)

T1.heights = lambda: calcularAlturas(T1)
T2.heights = lambda: calcularAlturas(T2)


T1.is_leaf = lambda v: len(list(T1.N(v))) == 0
T2.is_leaf = lambda v: len(list(T2.N(v))) == 0

print(calcularAlturas(T1))
print(calcularAlturas(T2))