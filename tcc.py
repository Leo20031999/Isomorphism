from collections import deque
import networkx as nx
from networkx.algorithms import bipartite
from structures.Grafo import Grafo

def OptimalTopDownCommonSubtree(T1, v, T2, u, M_prime, parent_map_T1, parent_map_T2, memo=None):
    if memo is None:
        memo = {}
    if (v, u) in memo:
        return memo[(v, u)]
    
    if T1.is_leaf(v, parent_map_T1) and T2.is_leaf(u, parent_map_T2):
        dist = 0 
    elif T1.is_leaf(v, parent_map_T1) or T2.is_leaf(u, parent_map_T2):
        dist = max(T1.altura(v, parent_map_T1), T2.altura(u, parent_map_T2))
        memo[(v, u)] = dist
        M_prime.add((v, u))
        return dist
    
    L = sorted(T1.vizinhanca(v, parent_map_T1), key=lambda x: (-T1.altura(x, parent_map_T1), x))
    R = sorted(T2.vizinhanca(u, parent_map_T2), key=lambda x: (-T2.altura(x, parent_map_T2), x))
    
    L_prefixed = [f"T1_{l}" for l in L]
    R_prefixed = [f"T2_{r}" for r in R]
    
    left_partition = L_prefixed.copy()
    right_partition = R_prefixed.copy()
    
    dummy_left = []
    dummy_right = []
    dummy_prefix = f"dummy_{v}_{u}_"

    if len(L) > len(R):
        num_dummies = len(L) - len(R)
        dummy_right = [f"{dummy_prefix}R_{i}" for i in range(num_dummies)]
        right_partition += dummy_right
    elif len(R) > len(L):
        num_dummies = len(R) - len(L)
        dummy_left = [f"{dummy_prefix}L_{i}" for i in range(num_dummies)]
        left_partition += dummy_left
    
    B = nx.Graph()
    B.add_nodes_from(left_partition, bipartite=0)
    B.add_nodes_from(right_partition, bipartite=1)
    
    for i, l in enumerate(L):
        for j, r in enumerate(R):
            l_pref = L_prefixed[i]
            r_pref = R_prefixed[j]
            weight = OptimalTopDownCommonSubtree(T1, l, T2, r, M_prime, parent_map_T1, parent_map_T2, memo)
            B.add_edge(l_pref, r_pref, weight=weight)
    
    current_height_T1 = T1.altura(v, parent_map_T1)
    current_height_T2 = T2.altura(u, parent_map_T2)
    
    for l_pref in L_prefixed:
        for dr in dummy_right:
            B.add_edge(l_pref, dr, weight=current_height_T1)
    
    for r_pref in R_prefixed:
        for dl in dummy_left:
            B.add_edge(dl, r_pref, weight=current_height_T2)
    
    for dl in dummy_left:
        for dr in dummy_right:
            B.add_edge(dl, dr, weight=0)
    
    Mvu = SolveOptimalPerfectMatching(B, left_partition)
    
    real_pairs = [
        (int(x.split("_")[1]), int(y.split("_")[1]))
        for x, y in Mvu
        if x.startswith("T1_") and y.startswith("T2_")
    ]
    
    M_prime.update(real_pairs)
    
    max_dist = max(B[x][y]['weight'] for x, y in Mvu) if Mvu else 0
    height_diff = abs(current_height_T1 - current_height_T2)
    final_dist = max(max_dist, height_diff)
    
    memo[(v, u)] = final_dist
    return final_dist

def SolveOptimalPerfectMatching(G, left_nodes):
    edges = sorted(G.edges(data='weight'), key=lambda x: x[2])
    
    if not edges:
        return set()
    
    unique_weights = sorted({w for _, _, w in edges})
    best_matching = None
    low, high = 0, len(unique_weights) - 1
    
    while low <= high:
        mid = (low + high) // 2
        threshold = unique_weights[mid]
        subG = nx.Graph()
        
        for u, v, w in edges:
            if w <= threshold:
                subG.add_edge(u, v)
        
        left_sub = [n for n in left_nodes if n in subG]
        right_sub = [n for n in G.nodes if n not in left_nodes and n in subG]
        
        if not left_sub or not right_sub:
            low = mid + 1
            continue
        
        try:
            matching = bipartite.hopcroft_karp_matching(subG, left_sub)
            matched_left = {u for u in matching if u in left_sub}
            if len(matched_left) == len(left_sub):
                best_matching = matching
                high = mid - 1
            else:
                low = mid + 1
        except nx.NetworkXError:
            low = mid + 1
    
    if best_matching is None:
        return set()
    
    return {(u, v) for u, v in best_matching.items() 
            if u in left_sub and u.startswith("T1_") and v.startswith("T2_")}

def ProcedureReconstructionOfMapping(T1, T2, r1, r2, M_prime, parent_map_T1, parent_map_T2):
    mapping = {r1: r2}
    used_T2 = {r2} 
    queue = deque([r1])
    
    while queue:
        current = queue.popleft()
        for child in T1.vizinhanca(current, parent_map_T1):
            if child in mapping:
                continue
            parent_in_T2 = mapping[current]
            candidates = [
                (vp, wp) for (vp, wp) in M_prime 
                if vp == child and parent_map_T2.get(wp, None) == parent_in_T2
                and wp not in used_T2
            ]
            if candidates:
                candidate = min(candidates, key=lambda x: abs(x[0] - x[1]))
                mapping[child] = candidate[1]
                used_T2.add(candidate[1])
                queue.append(child)
    return {(v, u) for v, u in mapping.items()}

def compute_parent_map(T, root):
    if root not in T.vertices():
        raise ValueError(f"Root {root} não existe na árvore.")
    parent_map = {root: None}
    queue = deque([root])
    
    while queue:
        u = queue.popleft()
        for v in T.vizinhanca(u):
            if v != parent_map.get(u) and v not in parent_map:
                parent_map[v] = u
                queue.append(v)
    
    return parent_map

def HausdorffDistanceBetweenTrees(T1, T2):
    hd = float('inf')
    best_mapping = set()
    
    t1_centers = sorted(T1.center())
    r1 = t1_centers[0] if t1_centers else next(iter(T1.vertices()), None)
    if r1 is None:
        return hd, best_mapping
    
    t2_centers = sorted(T2.center())
    valid_t2_roots = [c for c in t2_centers if c in T2.vertices()]
    
    for r2 in valid_t2_roots:
        parent_map_T1 = compute_parent_map(T1, r1)
        parent_map_T2 = compute_parent_map(T2, r2)
        M_prime = set()
        distance = OptimalTopDownCommonSubtree(T1, r1, T2, r2, M_prime, parent_map_T1, parent_map_T2)
        
        if distance < hd:
            hd = distance
            best_r2 = r2
            best_M_prime = M_prime.copy()
    
    if hd != float('inf'):
        parent_map_T2 = compute_parent_map(T2, best_r2)
        best_mapping = ProcedureReconstructionOfMapping(T1, T2, r1, best_r2, best_M_prime, parent_map_T1, parent_map_T2)
    
    return hd, best_mapping

def test_HausdorffDistance_arvores_identicas():
    T1 = Grafo()
    T1.adicionar_aresta(1, 2)
    T1.adicionar_aresta(1, 3)

    T2 = Grafo()
    T2.adicionar_aresta(1, 2)
    T2.adicionar_aresta(1, 3)

    hd, mapping = HausdorffDistanceBetweenTrees(T1, T2)
    print("\nTeste 1 - Árvores Idênticas:")
    print(f"Distância: {hd}")  # Esperado: 0
    print(f"Mapeamento: {mapping}")
    assert hd == 0 and len(mapping) == 3

def test_HausdorffDistance_altura_diferente():
    T1 = Grafo()
    T1.adicionar_aresta(1, 2)
    T1.adicionar_aresta(1, 3)

    T2 = Grafo()
    T2.adicionar_aresta(1, 2)
    T2.adicionar_aresta(1, 3)
    T2.adicionar_aresta(2, 4)

    hd, mapping = HausdorffDistanceBetweenTrees(T1, T2)

    print("\nTeste 2 - Alturas Diferentes:")
    print(f"Distância: {hd}")  # Esperado: 1
    print(f"Mapeamento: {mapping}")

    assert hd == 1, f"Erro: Distância esperada 1, obtida {hd}"
    assert len(mapping) == 3, "Mapeamento incompleto"

def test_HausdorffDistance_assimetricas():
    T1 = Grafo()
    T1.adicionar_aresta(1, 2)
    T1.adicionar_aresta(1, 3)
    T1.adicionar_aresta(2, 4)  # Altura de 1: 3 (1-2-4)

    T2 = Grafo()
    T2.adicionar_aresta(5, 6)
    T2.adicionar_aresta(6, 7)
    T2.adicionar_aresta(7, 8)
    T2.adicionar_aresta(8, 9)  # Altura de 7: 3 (7-8-9)
    T2.adicionar_aresta(9, 10)  # Nova aresta para altura 4

    hd, mapping = HausdorffDistanceBetweenTrees(T1, T2)
    print("\nTeste 3 - Árvores Assimétricas Modificado:")
    print(f"Distância: {hd}")  # Esperado: 1
    print(f"Mapeamento: {mapping}")
    assert hd == 1, f"Erro: Distância esperada 1, obtida {hd}"
    assert len(mapping) >= 2, "Mapeamento deve incluir raiz e filhos imediatos"

def test_HausdorffDistance_folhas_nao_correspondentes():
    T1 = Grafo()
    T1.adicionar_aresta(1, 2)
    T1.adicionar_aresta(1, 3)

    T2 = Grafo()
    T2.adicionar_aresta(1, 4)
    T2.adicionar_aresta(1, 5)

    hd, mapping = HausdorffDistanceBetweenTrees(T1, T2)

    print("\nTeste 4 - Folhas Não Correspondentes:")
    print(f"Distância: {hd}")  # Esperado: 0
    print(f"Mapeamento: {mapping}")

    assert hd == 0, f"Erro: Distância esperada 0, obtida {hd}"
    assert len(mapping) == 3, "Mapeamento incompleto"

def test_HausdorffDistance_um_no():
    # Árvore 1: [1]   Árvore 2: [2]
    T1 = Grafo()
    T1.adicionar_aresta(1, 1)  # Adiciona nó 1 (self-loop)
    T1.remover_aresta(1, 1)    # Remove a aresta, deixando o nó 1 isolado

    T2 = Grafo()
    T2.adicionar_aresta(2, 2)  # Adiciona nó 2 (self-loop)
    T2.remover_aresta(2, 2)    # Remove a aresta, deixando o nó 2 isolado

    hd, mapping = HausdorffDistanceBetweenTrees(T1, T2)
    print("\nTeste 5 - Um Nó:")
    print(f"Distância: {hd}")  # Esperado: 0
    print(f"Mapeamento: {mapping}")  # Ex: {(1, 2)}
    assert hd == 0, f"Erro: Distância esperada 0, obtida {hd}"
    assert len(mapping) == 1, "Mapeamento incompleto"

if __name__ == "__main__":
    test_HausdorffDistance_arvores_identicas()
    test_HausdorffDistance_altura_diferente()
    test_HausdorffDistance_assimetricas()
    test_HausdorffDistance_folhas_nao_correspondentes()
    test_HausdorffDistance_um_no()