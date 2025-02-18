from structures.Grafo import Grafo
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

            estado_parcial = {"distância": 0, "Mp": Mp.copy(), "pares_subárvores": {}}
            pilha.append((v, u, estado_parcial, "resolvido"))

            L_orig = T1.vizinhanca(v)
            R_orig = T2.vizinhanca(u)
            L_map = {node: idx for idx, node in enumerate(L_orig, start=1)}
            R_map = {node: idx for idx, node in enumerate(R_orig, start=len(L_orig) + 1)}
            total = len(L_orig) + len(R_orig)
            
            G_aux = nx.Graph()
            for i in range(1, total + 3):
                G_aux.add_node(i)
            
            for orig_l, new_l in L_map.items():
                for orig_r, new_r in R_map.items():
                    G_aux.add_edge(new_l, new_r, weight=1)
                    iteracoes += 1
                    if iteracoes >= MAX_ITERACOES:
                        return memo.get((v, u), ("Limite de iterações atingido", Mp))
            
            dum_left, dum_right = total + 1, total + 2
            for orig_l, new_l in L_map.items():
                G_aux.add_edge(new_l, dum_right, weight=1)
            for orig_r, new_r in R_map.items():
                G_aux.add_edge(dum_left, new_r, weight=1)
            
            Gvu = Grafo()
            Gvu.grafo = G_aux 
            Gvu.left_size = len(L_map)
            Gvu.right_size = len(R_map)
            
            estado_parcial["L_map"] = L_map
            estado_parcial["R_map"] = R_map
            estado_parcial["grafo_aux"] = Gvu
            
            for orig_l in L_orig:
                for orig_r in R_orig:
                    if (orig_l, orig_r) not in visited:
                        pilha.append((orig_l, orig_r, estado_parcial, "início"))
                        estado_parcial["pares_subárvores"][(orig_l, orig_r)] = (L_map.get(orig_l), R_map.get(orig_r))
                        visited.add((orig_l, orig_r))
        
        elif op == "resolvido":
            Gvu = estado_parcial["grafo_aux"]
            total = len(estado_parcial["L_map"]) + len(estado_parcial["R_map"])
            dum_left, dum_right = total + 1, total + 2

            print(f"Arestas do Gvu: {list(Gvu.grafo.edges(data=True))}")
            print(f"Pesos das arestas: {[data.get('weight', None) for _,_,data in Gvu.grafo.edges(data=True)]}")
            
            Mvu = solve_optimal_perfect_matching(Gvu)
            distância = max((Gvu.get_peso(u_edge, v_edge) for u_edge, v_edge in Mvu), default=0)
            L_set = set(estado_parcial["L_map"].values())
            R_set = set(estado_parcial["R_map"].values())
            Mvu_filtered = {(u_edge, v_edge) for u_edge, v_edge in Mvu if u_edge in L_set and v_edge in R_set}
            
            Mp.update(Mvu_filtered)
            memo[(v, u)] = (distância, Mp)

    return memo[(v, u)], iteracoes

def solve_optimal_perfect_matching(gvu):
    if not isinstance(gvu, Grafo):
        raise TypeError("O parâmetro 'gvu' deve ser uma instância de Grafo.")
    
    G = gvu.grafo  
    todos = sorted(G.nodes())
    if not todos:
        return set()
    
    dummy = set(todos[-2:])
    valid_nodes = [v for v in todos if v not in dummy]
    
    k = gvu.left_size  
    valid_nodes_sorted = sorted(valid_nodes)
    left_set = set(valid_nodes_sorted[:k])
    right_set = set(valid_nodes_sorted[k:])
    print(f"solve_optimal_perfect_matching: left_set: {left_set}")
    print(f"solve_optimal_perfect_matching: right_set: {right_set}")
    
    matching = nx.bipartite.maximum_matching(G, top_nodes=left_set)
    matching_set = {(u, v) for u, v in matching.items() if u in left_set and v in right_set}
    print(f"Melhor emparelhamento encontrado: {matching_set}")
    return matching_set

def compute_parent_map(T, root):
    """
    Computa um dicionário que mapeia cada vértice ao seu pai na árvore enraizada,
    utilizando uma busca em largura (BFS).
    Assume que T.vizinhanca(v) retorna os vizinhos de v.
    """
    parent = {root: None}
    queue = [root]
    visited = {root}
    while queue:
        u = queue.pop(0)
        for v in T.vizinhanca(u):
            if v not in visited:
                visited.add(v)
                parent[v] = u
                queue.append(v)
    return parent

def ProcedureReconstructionOfMapping(T1, T2, r1, r2, M_prime, M):
    """
    Reconstrói o mapeamento de isomorfismo a partir do conjunto M_prime.

    Parâmetros:
      - T1: Instância de Grafo que representa a árvore 1. Deve possuir os métodos
            preorder(root) e vizinhanca(v).
      - T2: Instância de Grafo que representa a árvore 2.
      - r1: Vértice raiz de T1.
      - r2: Vértice raiz de T2.
      - M_prime: Conjunto de pares (v, w) que são soluções dos emparelhamentos perfeitos.
      - M: Conjunto (mapeamento corrente) que será atualizado.

    Retorna:
      - M: O mapeamento reconstruído.
      
    A lógica é a seguinte:
      1. Adiciona (r1, r2) em M.
      2. Calcula os mapas de pais para T1 e T2 usando compute_parent_map.
      3. Obtém a travessia em pré-ordem de T1 a partir de r1.
      4. Para cada vértice v em T1, se existe um par (v, w) em M_prime e se os pais de v e w
         já estão mapeados em M, adiciona (v, w) em M.
    """
    M.add((r1, r2))
    
    parent_map_T1 = compute_parent_map(T1, r1)
    parent_map_T2 = compute_parent_map(T2, r2)
    print("Parent map T1:", parent_map_T1)
    print("Parent map T2:", parent_map_T2)
    
    preorder = T1.preorder(r1)
    print("Preorder traversal of T1:", preorder)
    
    for v in preorder:
        for (v_prime, w) in M_prime:
            if v == v_prime:
                parent_v = parent_map_T1.get(v)
                parent_w = parent_map_T2.get(w)
                print(f"Verificando par ({v}, {w}), pais: ({parent_v}, {parent_w})")
                if parent_v is not None and parent_w is not None:
                    if (parent_v, parent_w) in M:
                        M.add((v, w))
                        print(f"Adicionado par: ({v}, {w})")
    return M

def HausdorffDistanceBetweenTrees(T1, T2):
    """
    Algorithm1: HausdorffDistanceBetweenTrees
    Input: Árvores arbitrárias T1 e T2, onde diam(T1) >= diam(T2).
    Output: A distância de Hausdorff entre T1 e T2 (hd) e a estrutura de subárvore comum (M).
    
    Requisitos:
      - T1.center() retorna um vértice central de T1.
      - T1.compute_heights(root) e T2.compute_heights(root) computam as alturas dos vértices da árvore,
        a partir da raiz dada.
      - T1.vertices() e T2.vertices() retornam as listas de vértices.
      - OptimalTopDownCommonSubtree(T1, r1, T2, u, M_prime) retorna uma tupla (distance, iteracoes).
      - ReconstructionOfMapping(T1, r1, r2, O, M) reconstrói o mapeamento final em M.
    """
    hd = float('inf')
    O = set()
    r1 = T1.center()
    
    T1.altura(r1)
    
    r2 = None
    for u in T2.vertices():
        M_prime = set()
        T2.altura(u)
        (result, _) = OptimalTopDownCommonSubtree(T1, r1, T2, u, M_prime)
        distance = result[0]
        if distance < hd:
            hd = distance
            r2 = u
            O = M_prime.copy()
    
    M = set()
    ProcedureReconstructionOfMapping(T1, T2, r1, r2, O, M)
    
    return hd, M

def test_ProcedureReconstructionOfMapping():
    # Cria os grafos T1 e T2 usando a classe Grafo
    T1 = Grafo()
    T2 = Grafo()
    
    # Construa T1 com 5 vértices e a estrutura:
    #         1
    #        / \
    #       2   3
    #       |
    #       4
    #       |
    #       5
    T1.definir_n(5)
    T1.adicionar_aresta(1, 2)
    T1.adicionar_aresta(1, 3)
    T1.adicionar_aresta(2, 4)
    T1.adicionar_aresta(4, 5)
    
    # Construa T2 com 4 vértices e a estrutura:
    #         10
    #        /  \
    #       20   30
    #       |
    #       40
    T2.definir_n(4)
    T2.grafo.clear()  # Limpa os nós criados por definir_n
    for v in [10, 20, 30, 40]:
        T2.grafo.add_node(v)
    T2.num_nos = 4
    T2.adicionar_aresta(10, 20)
    T2.adicionar_aresta(10, 30)
    T2.adicionar_aresta(20, 40)
    
    # Defina M_prime: suponha que a subárvore comum seja T1: {2,4,5} e T2: {20,40}.
    # Assim, M_prime conterá os pares (2,20) e (4,40).
    M_prime = {(2, 20), (4, 40)}
    M = set()
    
    # Para obter o mapeamento desejado, usamos as raízes desejadas: r1 = 1 e r2 = 10.
    mapping_reconstruido = ProcedureReconstructionOfMapping(T1, T2, 1, 10, M_prime, M)
    print("Mapping reconstruído:", mapping_reconstruido)

if __name__ == "__main__":
    test_ProcedureReconstructionOfMapping()