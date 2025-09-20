# A polynomial time algorithm for simple undirected graph isomorphism
# He, J., Chen, J., Huang, G., Cao, J., Zhang, Z., Zheng, H., ... & Van Zundert, A. (2021).

from structures.Grafo import Grafo
from typing import List, Tuple, Any, Dict
import itertools


def are_isomorphic_molecular(g1: Grafo, g2: Grafo) -> bool:
    """Verifica isomorfismo molecular entre g1 e g2 com verificação adicional de Weisfeiler-Lehman."""
    # Verificação básica
    if not (len(g1.vertices()) == len(g2.vertices()) and len(g1.arestas()) == len(g2.arestas())):
        return False

    vlabels1 = [_norm_label(g1.get_rotulo_vertice(v)) for v in g1.vertices()]
    vlabels2 = [_norm_label(g2.get_rotulo_vertice(v)) for v in g2.vertices()]
    if sorted(vlabels1) != sorted(vlabels2):
        return False

    elabels1 = [_norm_label(g1.get_rotulo_aresta(u, v)) for (u, v) in g1.arestas()]
    elabels2 = [_norm_label(g2.get_rotulo_aresta(u, v)) for (u, v) in g2.arestas()]
    if sorted(elabels1) != sorted(elabels2):
        return False

    degree_seq1 = sorted([g1.grau(v) for v in g1.vertices()])
    degree_seq2 = sorted([g2.grau(v) for v in g2.vertices()])
    if degree_seq1 != degree_seq2:
        return False

    atom_label_to_id = _build_label_map(set(vlabels1) | set(vlabels2))
    bond_label_to_id = _build_label_map(set(elabels1) | set(elabels2))

    tt1 = _build_triple_tuple(g1)
    tt2 = _build_triple_tuple(g2)

    K1 = _calculate_vertex_encoded_sequence(tt1, g1, atom_label_to_id)
    K2 = _calculate_vertex_encoded_sequence(tt2, g2, atom_label_to_id)
    if not _verify_permutation(K1, K2):
        return False

    Eadj1 = _build_edge_adjacency_matrix(tt1)
    Eadj2 = _build_edge_adjacency_matrix(tt2)
    Lraw1 = _row_sums_matrix(Eadj1)
    Lraw2 = _row_sums_matrix(Eadj2)

    L1 = _encode_edge_sequence_with_bond_labels(tt1, Lraw1, g1, bond_label_to_id)
    L2 = _encode_edge_sequence_with_bond_labels(tt2, Lraw2, g2, bond_label_to_id)

    if not _verify_permutation(L1, L2):
        return False

    if not weisfeiler_lehman_check(g1, g2, iterations=3):
        return False

    return True


def weisfeiler_lehman_check(g1: Grafo, g2: Grafo, iterations: int = 3) -> bool:
    """Verifica isomorfismo usando o algoritmo de Weisfeiler-Lehman."""
    color1 = {}
    color2 = {}

    for v in g1.vertices():
        color1[v] = (g1.get_rotulo_vertice(v), g1.grau(v))

    for v in g2.vertices():
        color2[v] = (g2.get_rotulo_vertice(v), g2.grau(v))

    if sorted(color1.values()) != sorted(color2.values()):
        return False

    for _ in range(iterations):
        all_colors = {}
        new_color1 = {}
        new_color2 = {}

        for v in g1.vertices():
            neighbors = g1.vizinhanca(v)
            neighbor_colors = []
            for u in neighbors:
                rotulo = g1.get_rotulo_aresta(v, u)
                neighbor_colors.append((rotulo, color1[u]))
            neighbor_colors.sort(key=lambda x: (str(x[0]), str(x[1])))
            new_color_repr = (color1[v], tuple(neighbor_colors))
            new_color1[v] = new_color_repr
            if new_color_repr not in all_colors:
                all_colors[new_color_repr] = len(all_colors) + 1

        for v in g2.vertices():
            neighbors = g2.vizinhanca(v)
            neighbor_colors = []
            for u in neighbors:
                rotulo = g2.get_rotulo_aresta(v, u)
                neighbor_colors.append((rotulo, color2[u]))
            neighbor_colors.sort(key=lambda x: (str(x[0]), str(x[1])))
            new_color_repr = (color2[v], tuple(neighbor_colors))
            new_color2[v] = new_color_repr
            if new_color_repr not in all_colors:
                all_colors[new_color_repr] = len(all_colors) + 1

        int_color1 = {}
        int_color2 = {}

        for v in g1.vertices():
            int_color1[v] = all_colors[new_color1[v]]

        for v in g2.vertices():
            int_color2[v] = all_colors[new_color2[v]]

        if sorted(int_color1.values()) != sorted(int_color2.values()):
            return False

        color1 = int_color1
        color2 = int_color2

    return True

# -------- Funções auxiliares existentes --------
def _norm_label(x: Any) -> str:
    return "__NONE__" if x is None else str(x)


def _build_label_map(labels: set) -> dict:
    return {lab: i + 1 for i, lab in enumerate(sorted(labels, key=lambda s: str(s)))}

def _build_triple_tuple(grafo: Grafo) -> List[Tuple[int, Any, Any]]:
    return [(idx, u, v) for idx, (u, v) in enumerate(grafo.arestas(), start=1)]

def _calculate_vertex_encoded_sequence(triple_tuple, grafo: Grafo, atom_label_to_id: dict) -> List[int]:
    verts = grafo.vertices()
    degree = {v: 0 for v in verts}
    for _, u, v in triple_tuple:
        degree[u] = degree.get(u, 0) + 1
        degree[v] = degree.get(v, 0) + 1

    max_atom_id = max(atom_label_to_id.values()) if atom_label_to_id else 1
    base = max_atom_id + 1

    seq = []
    for v in verts:
        lab = _norm_label(grafo.get_rotulo_vertice(v))
        lab_id = atom_label_to_id.get(lab, 0)
        encoded = degree.get(v, 0) * base + lab_id
        seq.append(int(encoded))
    return seq


def _build_edge_adjacency_matrix(triple_tuple) -> List[List[int]]:
    m = len(triple_tuple)
    edges = [(u, v) for (_, u, v) in triple_tuple]
    mat = [[0] * m for _ in range(m)]
    for i in range(m):
        ui, vi = edges[i]
        for j in range(i + 1, m):
            uj, vj = edges[j]
            if ui == uj or ui == vj or vi == uj or vi == vj:
                mat[i][j] = 1
                mat[j][i] = 1
    return mat


def _row_sums_matrix(mat: List[List[int]]) -> List[int]:
    return [sum(row) for row in mat]


def _encode_edge_sequence_with_bond_labels(triple_tuple, Lraw: List[int], grafo: Grafo, bond_label_to_id: dict) -> List[
    int]:
    max_bond_id = max(bond_label_to_id.values()) if bond_label_to_id else 1
    base = max_bond_id + 1
    encoded = []
    for idx, (_, u, v) in enumerate(triple_tuple):
        bond_lab = _norm_label(grafo.get_rotulo_aresta(u, v))
        bond_id = bond_label_to_id.get(bond_lab, 0)
        lraw = Lraw[idx] if idx < len(Lraw) else 0
        encoded.append(int(lraw * base + bond_id))
    return encoded


def _verify_permutation(seq1: List[int], seq2: List[int]) -> bool:
    if len(seq1) != len(seq2):
        return False
    if not ((sum(seq1) == sum(seq2)) and (sum(x * x for x in seq1) == sum(y * y for y in seq2))):
        return False
    return sorted(seq1) == sorted(seq2)

def build_methane():
    g = Grafo()
    g.adicionar_aresta("C", "H1", rotulo="single")
    g.adicionar_aresta("C", "H2", rotulo="single")
    g.adicionar_aresta("C", "H3", rotulo="single")
    g.adicionar_aresta("C", "H4", rotulo="single")
    g.set_rotulo_vertice("C", "C")
    for h in ["H1", "H2", "H3", "H4"]:
        g.set_rotulo_vertice(h, "H")
    return g


def build_ethane():
    g = Grafo()
    g.adicionar_aresta("C1", "C2", rotulo="single")
    for i in range(1, 4):
        g.adicionar_aresta("C1", f"H{i}", rotulo="single")
    for i in range(4, 7):
        g.adicionar_aresta("C2", f"H{i}", rotulo="single")
    g.set_rotulo_vertice("C1", "C")
    g.set_rotulo_vertice("C2", "C")
    for i in range(1, 7):
        g.set_rotulo_vertice(f"H{i}", "H")
    return g


def build_ethene():
    g = Grafo()
    g.adicionar_aresta("C1", "C2", rotulo="double")
    g.adicionar_aresta("C1", "H1", rotulo="single")
    g.adicionar_aresta("C1", "H2", rotulo="single")
    g.adicionar_aresta("C2", "H3", rotulo="single")
    g.adicionar_aresta("C2", "H4", rotulo="single")
    g.set_rotulo_vertice("C1", "C")
    g.set_rotulo_vertice("C2", "C")
    for i in range(1, 5):
        g.set_rotulo_vertice(f"H{i}", "H")
    return g


def build_benzene():
    g = Grafo()
    # ciclo C6
    for i in range(6):
        g.adicionar_aresta(f"C{i}", f"C{(i+1)%6}", rotulo="single" if i % 2 == 0 else "double")
    # hidrogênios
    for i in range(6):
        g.adicionar_aresta(f"C{i}", f"H{i}", rotulo="single")
    for i in range(6):
        g.set_rotulo_vertice(f"C{i}", "C")
        g.set_rotulo_vertice(f"H{i}", "H")
    return g

def build_butane():
    """
    Butano linear: C1 - C2 - C3 - C4
    Cada carbono tem os hidrogênios suficientes para completar 4 ligações.
    """
    g = Grafo()
    # cadeia principal
    g.adicionar_aresta("C1", "C2", rotulo="single")
    g.adicionar_aresta("C2", "C3", rotulo="single")
    g.adicionar_aresta("C3", "C4", rotulo="single")

    # hidrogênios de C1
    for i in range(1, 4):
        g.adicionar_aresta("C1", f"H1{i}", rotulo="single")
    # hidrogênios de C2
    for i in range(1, 3):
        g.adicionar_aresta("C2", f"H2{i}", rotulo="single")
    # hidrogênios de C3
    for i in range(1, 3):
        g.adicionar_aresta("C3", f"H3{i}", rotulo="single")
    # hidrogênios de C4
    for i in range(1, 4):
        g.adicionar_aresta("C4", f"H4{i}", rotulo="single")

    # rótulos
    for c in ["C1", "C2", "C3", "C4"]:
        g.set_rotulo_vertice(c, "C")
    for i in range(1, 4):
        g.set_rotulo_vertice(f"H1{i}", "H")
    for i in range(1, 3):
        g.set_rotulo_vertice(f"H2{i}", "H")
    for i in range(1, 3):
        g.set_rotulo_vertice(f"H3{i}", "H")
    for i in range(1, 4):
        g.set_rotulo_vertice(f"H4{i}", "H")

    return g


def build_isobutane():
    """
    Isobutano (metilpropano): carbono central ligado a 3 outros carbonos.
    """
    g = Grafo()
    # carbono central C1 ligado a C2, C3, C4
    g.adicionar_aresta("C1", "C2", rotulo="single")
    g.adicionar_aresta("C1", "C3", rotulo="single")
    g.adicionar_aresta("C1", "C4", rotulo="single")

    # hidrogênios do C1
    g.adicionar_aresta("C1", "H1", rotulo="single")

    # hidrogênios de C2
    for i in range(1, 4):
        g.adicionar_aresta("C2", f"H2{i}", rotulo="single")
    # hidrogênios de C3
    for i in range(1, 3):
        g.adicionar_aresta("C3", f"H3{i}", rotulo="single")
    # hidrogênios de C4
    for i in range(1, 3):
        g.adicionar_aresta("C4", f"H4{i}", rotulo="single")

    # rótulos
    for c in ["C1", "C2", "C3", "C4"]:
        g.set_rotulo_vertice(c, "C")
    g.set_rotulo_vertice("H1", "H")
    for i in range(1, 4):
        g.set_rotulo_vertice(f"H2{i}", "H")
    for i in range(1, 3):
        g.set_rotulo_vertice(f"H3{i}", "H")
    for i in range(1, 3):
        g.set_rotulo_vertice(f"H4{i}", "H")

    return g

def build_pentane():
    """
    Pentano linear: C1-C2-C3-C4-C5
    """
    g = Grafo()
    # cadeia principal
    for i in range(1, 5):
        g.adicionar_aresta(f"C{i}", f"C{i+1}", rotulo="single")

    # hidrogênios C1 e C5 (CH3)
    for i in range(1, 4):
        g.adicionar_aresta("C1", f"H1{i}", rotulo="single")
        g.adicionar_aresta("C5", f"H5{i}", rotulo="single")

    # hidrogênios C2, C3, C4 (CH2)
    for i in range(1, 3):
        g.adicionar_aresta("C2", f"H2{i}", rotulo="single")
        g.adicionar_aresta("C3", f"H3{i}", rotulo="single")
        g.adicionar_aresta("C4", f"H4{i}", rotulo="single")

    # rótulos
    for c in [f"C{i}" for i in range(1, 6)]:
        g.set_rotulo_vertice(c, "C")
    for i in range(1, 4):
        g.set_rotulo_vertice(f"H1{i}", "H")
        g.set_rotulo_vertice(f"H5{i}", "H")
    for i in range(1, 3):
        g.set_rotulo_vertice(f"H2{i}", "H")
        g.set_rotulo_vertice(f"H3{i}", "H")
        g.set_rotulo_vertice(f"H4{i}", "H")
    return g


def build_isopentane():
    """
    Isopentano (metilbutano): carbono central com ramificação.
    """
    g = Grafo()
    # cadeia principal C1-C2-C3-C4
    for i in range(1, 4):
        g.adicionar_aresta(f"C{i}", f"C{i+1}", rotulo="single")
    # ramificação em C2 -> C5
    g.adicionar_aresta("C2", "C5", rotulo="single")

    # hidrogênios C1 e C4 (CH3)
    for i in range(1, 4):
        g.adicionar_aresta("C1", f"H1{i}", rotulo="single")
        g.adicionar_aresta("C4", f"H4{i}", rotulo="single")
    # hidrogênios C2 (CH)
    g.adicionar_aresta("C2", "H2", rotulo="single")
    # hidrogênios C3 (CH2)
    for i in range(1, 3):
        g.adicionar_aresta("C3", f"H3{i}", rotulo="single")
    # hidrogênios C5 (CH3)
    for i in range(1, 4):
        g.adicionar_aresta("C5", f"H5{i}", rotulo="single")

    # rótulos
    for c in [f"C{i}" for i in range(1, 6)]:
        g.set_rotulo_vertice(c, "C")
    for i in range(1, 4):
        g.set_rotulo_vertice(f"H1{i}", "H")
        g.set_rotulo_vertice(f"H4{i}", "H")
        g.set_rotulo_vertice(f"H5{i}", "H")
    g.set_rotulo_vertice("H2", "H")
    for i in range(1, 3):
        g.set_rotulo_vertice(f"H3{i}", "H")
    return g


def build_hexane():
    """
    Hexano linear: C1-C2-C3-C4-C5-C6
    """
    g = Grafo()
    for i in range(1, 6):
        g.adicionar_aresta(f"C{i}", f"C{i+1}", rotulo="single")

    # extremos (C1, C6 -> CH3)
    for i in range(1, 4):
        g.adicionar_aresta("C1", f"H1{i}", rotulo="single")
        g.adicionar_aresta("C6", f"H6{i}", rotulo="single")

    # internos (C2..C5 -> CH2)
    for c in [2, 3, 4, 5]:
        for i in range(1, 3):
            g.adicionar_aresta(f"C{c}", f"H{c}{i}", rotulo="single")

    # rótulos
    for c in [f"C{i}" for i in range(1, 7)]:
        g.set_rotulo_vertice(c, "C")
    for c in [1, 6]:
        for i in range(1, 4):
            g.set_rotulo_vertice(f"H{c}{i}", "H")
    for c in [2, 3, 4, 5]:
        for i in range(1, 3):
            g.set_rotulo_vertice(f"H{c}{i}", "H")
    return g


def build_cyclohexane():
    """
    Ciclohexano: anel C6 com 2 hidrogênios por carbono.
    """
    g = Grafo()
    # ciclo
    for i in range(6):
        g.adicionar_aresta(f"C{i}", f"C{(i+1)%6}", rotulo="single")
    # hidrogênios
    for i in range(6):
        g.adicionar_aresta(f"C{i}", f"H{i}a", rotulo="single")
        g.adicionar_aresta(f"C{i}", f"H{i}b", rotulo="single")

    # rótulos
    for i in range(6):
        g.set_rotulo_vertice(f"C{i}", "C")
        g.set_rotulo_vertice(f"H{i}a", "H")
        g.set_rotulo_vertice(f"H{i}b", "H")
    return g

if __name__ == "__main__":
    methane1 = build_methane()
    methane2 = build_methane()
    ethane = build_ethane()
    ethene = build_ethene()
    benzene1 = build_benzene()
    benzene2 = build_benzene()
    butane = build_butane()
    isobutane = build_isobutane()
    pentane = build_pentane()
    isopentane = build_isopentane()
    hexane = build_hexane()
    cyclohexane = build_cyclohexane()

    tests = [
        ("Metano vs Metano", methane1, methane2, True),
        ("Metano vs Etano", methane1, ethane, False),
        ("Etano vs Etano", ethane, build_ethane(), True),
        ("Etano vs Eteno", ethane, ethene, False),
        ("Benzeno vs Benzeno", benzene1, benzene2, True),
        ("Benzeno vs Eteno", benzene1, ethene, False),
        ("Butano vs Isobutano", butane, isobutane, False),
        ("Butano vs Butano", butane, build_butane(), True),
        ("Isobutano vs Isobutano", isobutane, build_isobutane(), True),
        ("Pentano vs Isopentano", pentane, isopentane, False),
        ("Pentano vs Pentano", pentane, build_pentane(), True),
        ("Isopentano vs Isopentano", isopentane, build_isopentane(), True),
        ("Hexano vs Ciclohexano", hexane, cyclohexane, False),
        ("Ciclohexano vs Ciclohexano", cyclohexane, build_cyclohexane(), True),
    ]

    for desc, g1, g2, expected in tests:
        result = are_isomorphic_molecular(g1, g2)
        status = "SUCCESS" if result == expected else "FAILURE"
        print(f"\nTeste: {desc}\nEsperado: {expected} | Obtido: {result}\nStatus: {status}")
