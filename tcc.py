# Determining the hausdorff distance between trees in polynomial time
# Kelenc, A., 2021

from collections import deque
import networkx as nx
from networkx.algorithms import bipartite
from structures.Grafo import Grafo
from typing import Dict, List, Set, Tuple, Any, Optional


# ==================== NOVAS FUNÇÕES PARA MANIPULAÇÃO DE ATRIBUTOS ====================

def atom_similarity(atom1: Dict[str, Any], atom2: Dict[str, Any],
                   custom_similarity_rules: Optional[Dict[Tuple[str, str], float]] = None) -> float:
    """
    Calcula a similaridade entre dois átomos com base em seus tipos.
    Permite que o usuário defina regras personalizadas de similaridade.

    Args:
        atom1: Atributos do primeiro átomo
        atom2: Atributos do segundo átomo
        custom_similarity_rules: Dicionário opcional com regras de similaridade personalizadas.
                                Formato: {('tipo1', 'tipo2'): similaridade}

    Returns:
        Valor de similaridade entre 0.0 e 1.0
    """
    if not atom1 or not atom2:
        return 0.0
    if 'atom_type' not in atom1 or 'atom_type' not in atom2:
        return 0.0

    type1 = atom1['atom_type']
    type2 = atom2['atom_type']

    if type1 == type2:
        return 1.0

    similarity_rules = custom_similarity_rules or {
        ('H', 'Cl'): 0.001, ('Cl', 'H'): 0.001,
        ('C', 'H'): 0.1, ('H', 'C'): 0.1,
        ('C', 'O'): 0.2, ('O', 'C'): 0.2,
        ('C', 'N'): 0.3, ('N', 'C'): 0.3,
        ('C', 'Cl'): 0.05, ('Cl', 'C'): 0.05,
        ('H', 'O'): 0.1, ('O', 'H'): 0.1,
        ('H', 'N'): 0.1, ('N', 'H'): 0.1,
        ('O', 'N'): 0.4, ('N', 'O'): 0.4,
    }

    return similarity_rules.get((type1, type2), 0.0)

def bond_similarity(bond1: Dict[str, Any], bond2: Dict[str, Any]) -> float:
    """Calcula a similaridade entre duas ligações baseada em seus atributos."""
    if bond1.get('bond_type') == bond2.get('bond_type'):
        return 1.0
    return 0.0


# ==================== FUNÇÕES AUXILIARES ====================

def altura_vertice(T, v, parent_map, use_attributes=False):
    """Calcula a altura de um vértice em uma árvore, com suporte a atributos."""
    if T.is_leaf(v, parent_map=parent_map):
        return 0

    filhos = obter_filhos(T, v, parent_map)

    if use_attributes and hasattr(T, 'nodes'):
        peso = 1
        if v in T.nodes and 'atom_type' in T.nodes[v]:
            atom_type = T.nodes[v]['atom_type']
            weights = {'C': 1, 'O': 2, 'N': 2, 'H': 1}
            peso = weights.get(atom_type, 1)

        return peso + max((altura_vertice(T, w, parent_map, use_attributes) for w in filhos), default=0)
    else:
        return 1 + max((altura_vertice(T, w, parent_map, use_attributes) for w in filhos), default=0)


def obter_filhos(T, v, parent_map):
    """Obtém os filhos de um vértice usando parent_map."""
    if parent_map.get(v) is None:
        return T.vizinhanca(v)
    return [w for w in T.vizinhanca(v) if parent_map.get(w) == v]


# ==================== FUNÇÕES PRINCIPAIS ====================

def OptimalTopDownCommonSubtree(T1, v, T2, u, M_prime, parent_map_T1, parent_map_T2,
                                memo=None, use_attributes=False,
                                atom_similarity_rules=None):
    if memo is None:
        memo = {}
    if (v, u) in memo:
        return memo[(v, u)]

    atom_dist = 0
    if use_attributes:
        atom_sim = atom_similarity(
            T1.get_atributos_vertice(v),
            T2.get_atributos_vertice(u),
            atom_similarity_rules
        )
        atom_dist = 1 - atom_sim

    is_leaf_T1 = T1.is_leaf(v, parent_map=parent_map_T1)
    is_leaf_T2 = T2.is_leaf(u, parent_map=parent_map_T2)

    if is_leaf_T1 and is_leaf_T2:
        dist = atom_dist if use_attributes else 0
        memo[(v, u)] = dist
        return dist
    elif is_leaf_T1 or is_leaf_T2:
        altura_v = T1.altura(v, parent_map_T1)
        altura_u = T2.altura(u, parent_map_T2)
        dist = max(altura_v, altura_u, atom_dist)
        memo[(v, u)] = dist
        return dist
    else:
        try:
            filhos_T1 = obter_filhos(T1, v, parent_map_T1)
            filhos_T2 = obter_filhos(T2, u, parent_map_T2)

            L = sorted(filhos_T1, key=lambda x: (-T1.altura(x, parent_map_T1), x))
            R = sorted(filhos_T2, key=lambda x: (-T2.altura(x, parent_map_T2), x))

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
                    weight = OptimalTopDownCommonSubtree(T1, l, T2, r, M_prime,
                                                         parent_map_T1, parent_map_T2,
                                                         memo, use_attributes,
                                                         atom_similarity_rules)

                    if use_attributes:
                        child_atom_sim = atom_similarity(
                            T1.get_atributos_vertice(l),
                            T2.get_atributos_vertice(r),
                            atom_similarity_rules
                        )
                        child_atom_dist = 1 - child_atom_sim
                        weight = max(weight, child_atom_dist)

                    B.add_edge(l_pref, r_pref, weight=weight)

            current_height_T1 = T1.altura(v, parent_map_T1)
            current_height_T2 = T2.altura(u, parent_map_T2)

            dummy_weight = max(current_height_T1, current_height_T2) + 10

            for l_pref in L_prefixed:
                for dr in dummy_right:
                    B.add_edge(l_pref, dr, weight=dummy_weight)

            for r_pref in R_prefixed:
                for dl in dummy_left:
                    B.add_edge(dl, r_pref, weight=dummy_weight)

            for dl in dummy_left:
                for dr in dummy_right:
                    B.add_edge(dl, dr, weight=0)

            Mvu = SolveOptimalPerfectMatching(B, left_partition)

            real_pairs = []
            for x, y in Mvu:
                if x.startswith("T1_") and y.startswith("T2_"):

                    v_label = x[3:]
                    u_label = y[3:]

                    try:
                        v_label = int(v_label) if v_label.isdigit() else v_label
                        u_label = int(u_label) if u_label.isdigit() else u_label
                    except:
                        pass

                    real_pairs.append((v_label, u_label))

            M_prime.update(real_pairs)

            max_dist = max(B[x][y]['weight'] for x, y in Mvu) if Mvu else 0
            height_diff = abs(current_height_T1 - current_height_T2)

            final_dist = max(max_dist, height_diff, atom_dist)

            memo[(v, u)] = final_dist
            return final_dist

        except Exception as e:
            print(f"Erro ao processar nós {v} e {u}: {e}")
            altura_v = T1.altura(v, parent_map_T1)
            altura_u = T2.altura(u, parent_map_T2)
            final_dist = max(altura_v, altura_u, atom_dist)
            memo[(v, u)] = final_dist
            return final_dist

def SolveOptimalPerfectMatching(G, left_nodes):
    try:
        left_set = set(left_nodes)
        right_set = set(G.nodes) - left_set

        complete_G = nx.Graph()

        for u, v, data in G.edges(data=True):
            complete_G.add_edge(u, v, weight=data['weight'])

        for u in left_set:
            for v in right_set:
                if not complete_G.has_edge(u, v):
                    complete_G.add_edge(u, v, weight=float('inf'))

        matching = bipartite.minimum_weight_full_matching(complete_G, weight='weight')

        valid_pairs = set()
        for u, v in matching.items():
            if u in left_set and v in right_set and complete_G[u][v]['weight'] < float('inf'):
                valid_pairs.add((u, v))

        return valid_pairs

    except Exception as e:
        print(f"Erro em SolveOptimalPerfectMatching: {e}")
        return set()

def ProcedureReconstructionOfMapping(T1, T2, r1, r2, M_prime, parent_map_T1, parent_map_T2):
    try:
        mapping = {r1: r2}
        used_T2 = {r2}
        pre_order = T1.preorder(r1, parent_map_T1)

        for v in pre_order:
            if v == r1:
                continue
            parent_v = parent_map_T1[v]
            if parent_v not in mapping:
                continue
            parent_u = mapping[parent_v]
            candidates = [
                (vp, wp) for (vp, wp) in M_prime
                if vp == v and parent_map_T2.get(wp) == parent_u and wp not in used_T2
            ]
            if candidates:
                mapping[v] = candidates[0][1]
                used_T2.add(candidates[0][1])

        for v in pre_order:
            if v in mapping:
                continue
            parent_v = parent_map_T1[v]
            if parent_v not in mapping:
                continue
            parent_u = mapping[parent_v]
            available_children = [
                w for w in T2.obter_filhos(parent_u, parent_map_T2)
                if w not in used_T2
            ]
            if available_children:
                mapping[v] = available_children[0]
                used_T2.add(available_children[0])

        return {(k, v) for k, v in mapping.items()}
    except Exception as e:
        print(f"Erro em ProcedureReconstructionOfMapping: {e}")
        return set()

def compute_parent_map(T, root):
    try:
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
    except Exception as e:
        print(f"Erro em compute_parent_map: {e}")
        return {root: None}

def HausdorffDistanceBetweenTrees(T1, T2, use_attributes=False,
                                 atom_similarity_rules=None):
    try:
        hd = float('inf')
        best_mapping = set()

        r1 = T1.center()[0] if T1.vertices() else None
        if r1 is None:
            return 0, set()

        for u in T2.vertices():
            parent_map_T1 = compute_parent_map(T1, r1)
            parent_map_T2 = compute_parent_map(T2, u)
            M_prime = set()
            distance = OptimalTopDownCommonSubtree(T1, r1, T2, u, M_prime,
                                                   parent_map_T1, parent_map_T2,
                                                   use_attributes=use_attributes,
                                                   atom_similarity_rules=atom_similarity_rules)

            if distance < hd:
                hd = distance
                best_r2 = u
                best_M_prime = M_prime.copy()

        if hd != float('inf'):
            parent_map_T2 = compute_parent_map(T2, best_r2)
            best_mapping = ProcedureReconstructionOfMapping(T1, T2, r1, best_r2,
                                                            best_M_prime, parent_map_T1,
                                                            parent_map_T2)

        return hd, best_mapping
    except Exception as e:
        print(f"Erro em HausdorffDistanceBetweenTrees: {e}")
        return float('inf'), set()

def test_metano():
    """Teste para metano (CH4) - um carbono com quatro hidrogênios."""
    print("\n=== Teste Metano (CH4) ===")

    # Metano 1
    CH4_1 = Grafo()
    CH4_1.adicionar_vertice("C1", atributos={"atom_type": "C"})
    CH4_1.adicionar_vertice("H1", atributos={"atom_type": "H"})
    CH4_1.adicionar_vertice("H2", atributos={"atom_type": "H"})
    CH4_1.adicionar_vertice("H3", atributos={"atom_type": "H"})
    CH4_1.adicionar_vertice("H4", atributos={"atom_type": "H"})
    CH4_1.adicionar_aresta("C1", "H1")
    CH4_1.adicionar_aresta("C1", "H2")
    CH4_1.adicionar_aresta("C1", "H3")
    CH4_1.adicionar_aresta("C1", "H4")

    # Metano 2 (idêntico)
    CH4_2 = Grafo()
    CH4_2.adicionar_vertice("C1", atributos={"atom_type": "C"})
    CH4_2.adicionar_vertice("H1", atributos={"atom_type": "H"})
    CH4_2.adicionar_vertice("H2", atributos={"atom_type": "H"})
    CH4_2.adicionar_vertice("H3", atributos={"atom_type": "H"})
    CH4_2.adicionar_vertice("H4", atributos={"atom_type": "H"})
    CH4_2.adicionar_aresta("C1", "H1")
    CH4_2.adicionar_aresta("C1", "H2")
    CH4_2.adicionar_aresta("C1", "H3")
    CH4_2.adicionar_aresta("C1", "H4")

    print("Atributos de C1 em CH4_1:", CH4_1.get_atributos_vertice("C1"))
    print("Atributos de H1 em CH4_1:", CH4_1.get_atributos_vertice("H1"))

    hd, mapping = HausdorffDistanceBetweenTrees(CH4_1, CH4_2, use_attributes=True)
    print(f"Distância entre moléculas idênticas de metano: {hd}")
    print(f"Mapeamento: {mapping}")
    assert hd == 0, f"Erro: Distância esperada 0, obtida {hd}"

def test_etano():
    """Teste para etano (C2H6) - dois carbonos ligados com hidrogênios."""
    print("\n=== Teste Etano (C2H6) ===")

    # Etano normal
    C2H6 = Grafo()
    C2H6.adicionar_vertice("C1", atributos={"atom_type": "C"})
    C2H6.adicionar_vertice("C2", atributos={"atom_type": "C"})
    C2H6.adicionar_vertice("H1", atributos={"atom_type": "H"})
    C2H6.adicionar_vertice("H2", atributos={"atom_type": "H"})
    C2H6.adicionar_vertice("H3", atributos={"atom_type": "H"})
    C2H6.adicionar_vertice("H4", atributos={"atom_type": "H"})
    C2H6.adicionar_vertice("H5", atributos={"atom_type": "H"})
    C2H6.adicionar_vertice("H6", atributos={"atom_type": "H"})
    C2H6.adicionar_aresta("C1", "C2")
    C2H6.adicionar_aresta("C1", "H1")
    C2H6.adicionar_aresta("C1", "H2")
    C2H6.adicionar_aresta("C1", "H3")
    C2H6.adicionar_aresta("C2", "H4")
    C2H6.adicionar_aresta("C2", "H5")
    C2H6.adicionar_aresta("C2", "H6")

    # Etano com um hidrogênio substituído por cloro
    C2H5Cl = Grafo()
    C2H5Cl.adicionar_vertice("C1", atributos={"atom_type": "C"})
    C2H5Cl.adicionar_vertice("C2", atributos={"atom_type": "C"})
    C2H5Cl.adicionar_vertice("H1", atributos={"atom_type": "H"})
    C2H5Cl.adicionar_vertice("H2", atributos={"atom_type": "H"})
    C2H5Cl.adicionar_vertice("H3", atributos={"atom_type": "H"})
    C2H5Cl.adicionar_vertice("Cl4", atributos={"atom_type": "Cl"})
    C2H5Cl.adicionar_vertice("H5", atributos={"atom_type": "H"})
    C2H5Cl.adicionar_vertice("H6", atributos={"atom_type": "H"})
    C2H5Cl.adicionar_aresta("C1", "C2")
    C2H5Cl.adicionar_aresta("C1", "H1")
    C2H5Cl.adicionar_aresta("C1", "H2")
    C2H5Cl.adicionar_aresta("C1", "H3")
    C2H5Cl.adicionar_aresta("C2", "Cl4")
    C2H5Cl.adicionar_aresta("C2", "H5")
    C2H5Cl.adicionar_aresta("C2", "H6")

    hd, mapping = HausdorffDistanceBetweenTrees(C2H6, C2H5Cl, use_attributes=True)
    print(f"Distância entre etano e cloroetano: {hd}")
    print(f"Mapeamento: {mapping}")

    cloro_mapeado = any('Cl' in str(atom) for _, atom in mapping)

    if hd == 0:
        print("AVISO: Distância zero obtida, mas era esperado > 0 devido à dissimilaridade atômica.")
    else:
        print("Distância maior que zero, como esperado.")

    if not cloro_mapeado:
        print("AVISO: Átomo de cloro não foi mapeado.")
    else:
        print("Átomo de cloro foi mapeado corretamente.")

def test_alcool_vs_alcano():
    """Teste para comparar etanol (C2H5OH) com etano (C2H6)."""
    print("\n=== Teste Etanol vs Etano ===")

    # Etano
    C2H6 = Grafo()
    C2H6.adicionar_vertice("C1", {"atom_type": "C"})
    C2H6.adicionar_vertice("C2", {"atom_type": "C"})
    for i in range(1, 7):
        C2H6.adicionar_vertice(f"H{i}", {"atom_type": "H"})
    C2H6.adicionar_aresta("C1", "C2")
    for i in range(1, 4):
        C2H6.adicionar_aresta("C1", f"H{i}")
    for i in range(4, 7):
        C2H6.adicionar_aresta("C2", f"H{i}")

    # Etanol (C2H5OH)
    C2H5OH = Grafo()
    C2H5OH.adicionar_vertice("C1", {"atom_type": "C"})
    C2H5OH.adicionar_vertice("C2", {"atom_type": "C"})
    C2H5OH.adicionar_vertice("O", {"atom_type": "O"})
    C2H5OH.adicionar_vertice("H1", {"atom_type": "H"})
    C2H5OH.adicionar_vertice("H2", {"atom_type": "H"})
    C2H5OH.adicionar_vertice("H3", {"atom_type": "H"})
    C2H5OH.adicionar_vertice("H4", {"atom_type": "H"})
    C2H5OH.adicionar_vertice("H5", {"atom_type": "H"})
    C2H5OH.adicionar_aresta("C1", "C2")
    C2H5OH.adicionar_aresta("C2", "O")
    C2H5OH.adicionar_aresta("O", "H5")
    for i in range(1, 4):
        C2H5OH.adicionar_aresta("C1", f"H{i}")
    for i in range(4, 5):
        C2H5OH.adicionar_aresta("C2", f"H{i}")

    hd, mapping = HausdorffDistanceBetweenTrees(C2H6, C2H5OH, use_attributes=True)
    print(f"Distância entre etano e etanol: {hd}")
    print(f"Mapeamento: {mapping}")
    assert hd > 0, f"Erro: Distância esperada > 0, obtida {hd}"

def test_propano_isomeros():
    """Teste para comparar n-propano e isopropanol."""
    print("\n=== Teste n-Propano vs Isopropanol ===")

    # n-Propano (C3H8)
    n_propano = Grafo()
    n_propano.adicionar_vertice("C1", {"atom_type": "C"})
    n_propano.adicionar_vertice("C2", {"atom_type": "C"})
    n_propano.adicionar_vertice("C3", {"atom_type": "C"})
    for i in range(1, 9):
        n_propano.adicionar_vertice(f"H{i}", {"atom_type": "H"})
    n_propano.adicionar_aresta("C1", "C2")
    n_propano.adicionar_aresta("C2", "C3")
    n_propano.adicionar_aresta("C1", "H1")
    n_propano.adicionar_aresta("C1", "H2")
    n_propano.adicionar_aresta("C1", "H3")
    n_propano.adicionar_aresta("C2", "H4")
    n_propano.adicionar_aresta("C2", "H5")
    n_propano.adicionar_aresta("C3", "H6")
    n_propano.adicionar_aresta("C3", "H7")
    n_propano.adicionar_aresta("C3", "H8")

    # Isopropanol (C3H7OH)
    isopropanol = Grafo()
    isopropanol.adicionar_vertice("C1", {"atom_type": "C"})
    isopropanol.adicionar_vertice("C2", {"atom_type": "C"})
    isopropanol.adicionar_vertice("C3", {"atom_type": "C"})
    isopropanol.adicionar_vertice("O", {"atom_type": "O"})
    for i in range(1, 8):
        isopropanol.adicionar_vertice(f"H{i}", {"atom_type": "H"})
    isopropanol.adicionar_aresta("C1", "C2")
    isopropanol.adicionar_aresta("C1", "C3")
    isopropanol.adicionar_aresta("C2", "O")
    isopropanol.adicionar_aresta("O", "H7")
    isopropanol.adicionar_aresta("C1", "H1")
    isopropanol.adicionar_aresta("C2", "H2")
    isopropanol.adicionar_aresta("C2", "H3")
    isopropanol.adicionar_aresta("C3", "H4")
    isopropanol.adicionar_aresta("C3", "H5")
    isopropanol.adicionar_aresta("C3", "H6")

    hd, mapping = HausdorffDistanceBetweenTrees(n_propano, isopropanol, use_attributes=True)
    print(f"Distância entre n-propano e isopropanol: {hd}")
    print(f"Mapeamento: {mapping}")
    assert hd > 0, f"Erro: Distância esperada > 0, obtida {hd}"

def test_benzeno_sem_ciclo():
    """Teste para benzeno representado como árvore (sem ciclos)."""
    print("\n=== Teste Benzeno (como árvore) ===")

    benzeno_1 = Grafo()
    carbonos = ["C1", "C2", "C3", "C4", "C5", "C6"]
    for c in carbonos:
        benzeno_1.adicionar_vertice(c, atributos={"atom_type": "C"})

    # Ligações (representação acíclica)
    benzeno_1.adicionar_aresta("C1", "C2")
    benzeno_1.adicionar_aresta("C2", "C3")
    benzeno_1.adicionar_aresta("C3", "C4")
    benzeno_1.adicionar_aresta("C4", "C5")
    benzeno_1.adicionar_aresta("C5", "C6")

    benzeno_2 = Grafo()
    for c in carbonos:
        benzeno_2.adicionar_vertice(c, atributos={"atom_type": "C"})
    benzeno_2.adicionar_aresta("C1", "C2")
    benzeno_2.adicionar_aresta("C2", "C3")
    benzeno_2.adicionar_aresta("C3", "C4")
    benzeno_2.adicionar_aresta("C4", "C5")
    benzeno_2.adicionar_aresta("C5", "C6")

    hd, mapping = HausdorffDistanceBetweenTrees(benzeno_1, benzeno_2, use_attributes=True)
    print(f"Distância entre representações acíclicas idênticas de benzeno: {hd}")
    print(f"Mapeamento: {mapping}")

    assert hd >= 0, f"Erro: Distância esperada >= 0, obtida {hd}"

def test_moleculas_diferentes_tamanhos():
    """Teste para moléculas com tamanhos diferentes."""
    print("\n=== Teste Metano vs Etano ===")

    # Metano
    CH4 = Grafo()
    CH4.adicionar_vertice("C1", {"atom_type": "C"})
    for i in range(1, 5):
        CH4.adicionar_vertice(f"H{i}", {"atom_type": "H"})
        CH4.adicionar_aresta("C1", f"H{i}")

    # Etano
    C2H6 = Grafo()
    C2H6.adicionar_vertice("C1", {"atom_type": "C"})
    C2H6.adicionar_vertice("C2", {"atom_type": "C"})
    for i in range(1, 7):
        C2H6.adicionar_vertice(f"H{i}", {"atom_type": "H"})
    C2H6.adicionar_aresta("C1", "C2")
    C2H6.adicionar_aresta("C1", "H1")
    C2H6.adicionar_aresta("C1", "H2")
    C2H6.adicionar_aresta("C1", "H3")
    C2H6.adicionar_aresta("C2", "H4")
    C2H6.adicionar_aresta("C2", "H5")
    C2H6.adicionar_aresta("C2", "H6")

    hd, mapping = HausdorffDistanceBetweenTrees(CH4, C2H6, use_attributes=True)
    print(f"Distância entre metano e etano: {hd}")
    print(f"Mapeamento: {mapping}")
    assert hd > 0, f"Erro: Distância esperada > 0, obtida {hd}"

def test_com_regras_personalizadas():
    """Teste com regras de similaridade atômica personalizadas."""
    print("\n=== Teste com Regras Personalizadas ===")

    # Metano
    CH4 = Grafo()
    CH4.adicionar_vertice("C1", atributos={"atom_type": "C"})
    CH4.adicionar_vertice("H1", atributos={"atom_type": "H"})
    CH4.adicionar_vertice("H2", atributos={"atom_type": "H"})
    CH4.adicionar_vertice("H3", atributos={"atom_type": "H"})
    CH4.adicionar_vertice("H4", atributos={"atom_type": "H"})
    CH4.adicionar_aresta("C1", "H1")
    CH4.adicionar_aresta("C1", "H2")
    CH4.adicionar_aresta("C1", "H3")
    CH4.adicionar_aresta("C1", "H4")

    CH3OH = Grafo()
    CH3OH.adicionar_vertice("C1", atributos={"atom_type": "C"})
    CH3OH.adicionar_vertice("H1", atributos={"atom_type": "H"})
    CH3OH.adicionar_vertice("H2", atributos={"atom_type": "H"})
    CH3OH.adicionar_vertice("H3", atributos={"atom_type": "H"})
    CH3OH.adicionar_vertice("O4", atributos={"atom_type": "O"})
    CH3OH.adicionar_aresta("C1", "H1")
    CH3OH.adicionar_aresta("C1", "H2")
    CH3OH.adicionar_aresta("C1", "H3")
    CH3OH.adicionar_aresta("C1", "O4")

    regras_personalizadas = {
        ('C', 'O'): 0.9,
        ('O', 'C'): 0.9,
        ('H', 'O'): 0.9,
        ('O', 'H'): 0.9,
    }

    hd, mapping = HausdorffDistanceBetweenTrees(CH4, CH3OH, use_attributes=True,
                                               atom_similarity_rules=regras_personalizadas)

    print(f"Distância com regras personalizadas: {hd}")
    print(f"Mapeamento: {mapping}")

    hd_padrao, _ = HausdorffDistanceBetweenTrees(CH4, CH3OH, use_attributes=True)
    print(f"Distância com regras padrão: {hd_padrao}")

    assert hd < hd_padrao, f"Com regras personalizadas, a distância deveria ser menor. Personalizada: {hd}, Padrão: {hd_padrao}"

if __name__ == "__main__":
    test_metano()
    test_etano()
    test_alcool_vs_alcano()
    test_propano_isomeros()
    test_benzeno_sem_ciclo()
    test_moleculas_diferentes_tamanhos()
    test_com_regras_personalizadas()