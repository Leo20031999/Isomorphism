from structures import graph
from structures import weight

def OptimalTopDownCommonSubtree(T1,v,T2,u,Mp):
    if T1.is_leaf(v) or T2.is_leaf(u):
        return max(T1.heights()[v], T2.heights()[u]), Mp
    Gvu = graph.Graph()
    for i in T1.adj_list.get(v,[]):
        for j in T2.adj_list.get(u, []):
            Gvu.add_edge(i,j)
    dum_v = "dum_v"
    dum_u = "dum_u"
    for x in T1.adj_list.get(v,[]):
        Gvu.add_edge(x, dum_u)
    for y in T2.adj_list.get(u,[]):
        Gvu.add_edge(dum_v, y)

    heightsT1 = T1.heights()
    heightsT2 = T2.heights()

    for e in Gvu.adj_list:
        u, v = e
        if u == dum_v:
            Gvu.adj_list[e] = heightsT2[v] + 1
        elif v == dum_u:
            Gvu.adj_list[e] = heightsT1[u] + 1
        else:
            Gvu.adj_list[e] = OptimalTopDownCommonSubtree(T1,u,T2,v,Mp)[0]
    
    Mvu = solveOptimalPatternMatching(Gvu)
    distance = max(Gvu.adj_list[e] for e in Mvu)
    Mvu = {e for e in Mvu if not (is_dummy_vertex(e[0]) or is_dummy_vertex(e[1]))}
    Mp.update(Mvu)
    return distance, Mp

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

def hausdorffDistanceBetweenTrees(grafo1, grafo2):
    return 0   

graph = graph.Graph()
graph.add_edge(0, 1)
graph.add_edge(0, 2)
graph.add_edge(1, 3)
graph.add_edge(1, 4)
graph.add_edge(2, 5)
graph.add_edge(2, 6)
graph.add_edge(5, 7)
graph.add_edge(6, 8)

print("Representação do Grafo:")
graph.print_graph()

heights = graph.heights()
print("\nAlturas dos vértices:")
for vertex, height in heights.items():
    print(f"Vértice {vertex}: Altura {height}")