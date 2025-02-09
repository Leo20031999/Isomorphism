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
    utilizando uma DFS e o método T.vizinhanca(v) para obter os vizinhos.
    """
    parent = {root: None}
    stack = [root]
    visited = {root}
    while stack:
        u = stack.pop()
        for v in T.vizinhanca(u):
            if v not in visited:
                visited.add(v)
                parent[v] = u
                stack.append(v)
    return parent

def ProcedureReconstructionOfMapping(T1, T2, r1, r2, M_prime, M):
    """
    Reconstrói o mapeamento de isomorfismo a partir do conjunto M_prime.
    
    Parâmetros:
      - T1: Instância de Grafo (ou ListaAdjG) que representa a árvore 1.
            Deve ter os métodos preorder(root) e vizinhanca(v).
      - T2: Instância de Grafo que representa a árvore 2.
      - r1: Vértice raiz de T1.
      - r2: Vértice raiz de T2.
      - M_prime: Conjunto de pares (v, w) (soluções dos emparelhamentos perfeitos).
      - M: Conjunto (mapeamento corrente) que será atualizado.
    
    Retorna:
      - M: O mapeamento reconstruído.
      
    O procedimento realiza:
      1. Adiciona (r1, r2) em M.
      2. Computa os mapas de pais para T1 e T2.
      3. Obtém a travessia em pré-ordem de T1 a partir de r1.
      4. Para cada vértice v na ordem de pré-ordem, se existir um par (v, w) em M_prime e se
         os pais de v e w já estiverem mapeados em M, adiciona (v, w) em M.
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

# Método de teste utilizando a classe Grafo (ou ListaAdjG)
def test_ProcedureReconstructionOfMapping():
    # Crie os grafos T1 e T2 utilizando sua classe Grafo.
    T1 = Grafo()
    T2 = Grafo()
    
    # Construa T1: árvore com 5 vértices e a estrutura:
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
    
    # Construa T2: árvore com 4 vértices e a estrutura:
    #         10
    #        /  \
    #       20   30
    #       |
    #       40
    T2.definir_n(4)
    # Para T2, limpamos os nós criados e adicionamos os desejados manualmente:
    T2.grafo.clear()
    for v in [10, 20, 30, 40]:
        T2.grafo.add_node(v)
    T2.num_nos = 4
    T2.adicionar_aresta(10, 20)
    T2.adicionar_aresta(10, 30)
    T2.adicionar_aresta(20, 40)
    
    # Defina M_prime com os pares que representam a subárvore comum:
    # Suponha que a subárvore comum seja T1: {2,4,5} e T2: {20,40}
    M_prime = {(2, 20), (4, 40)}
    M = set()
    
    mapping_reconstruido = ProcedureReconstructionOfMapping(T1, T2, 1, 10, M_prime, M)
    print("Mapping reconstruído:", mapping_reconstruido)

if __name__ == "__main__":
    test_ProcedureReconstructionOfMapping()
