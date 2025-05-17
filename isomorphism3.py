from structures.Grafo import Grafo
import numpy as np
import networkx as nx

def sao_isomorficos(g1: Grafo, g2: Grafo) -> bool:
    iso_g1 = len([v for v in g1.vertices() if g1.grau(v) == 0])
    iso_g2 = len([v for v in g2.vertices() if g2.grau(v) == 0])
    if iso_g1 != iso_g2:
        return False

    g1_clean = _filtrar_grafo(g1)
    g2_clean = _filtrar_grafo(g2)

    if not _verificacao_inicial(g1_clean, g2_clean):
        return False

    v1_sums, e1_sums = _calcular_propriedades(g1_clean)
    v2_sums, e2_sums = _calcular_propriedades(g2_clean)

    if not _sao_permutacao_otimizada(v1_sums, v2_sums):
        return False
    if not _sao_permutacao_otimizada(e1_sums, e2_sums):
        return False

    return _verificar_espectro_completo(g1_clean, g2_clean)

def _filtrar_grafo(grafo: Grafo) -> nx.Graph:
    g = grafo.grafo.copy()
    g.remove_nodes_from(list(nx.isolates(g)))
    return g

def _verificacao_inicial(g1: nx.Graph, g2: nx.Graph) -> bool:
    return (g1.number_of_nodes() == g2.number_of_nodes() and 
            g1.number_of_edges() == g2.number_of_edges() and
            sorted(d for _, d in g1.degree()) == sorted(d for _, d in g2.degree()))

def _calcular_propriedades(g: nx.Graph) -> tuple:
    vertex_sums = sorted(d for _, d in g.degree())
    edge_sums = sorted(g.degree(u) + g.degree(v) - 2 for u, v in g.edges())
    return (vertex_sums, edge_sums)

def _sao_permutacao_otimizada(arr1, arr2) -> bool:
    return sorted(arr1) == sorted(arr2)

def _verificar_espectro_completo(g1: nx.Graph, g2: nx.Graph) -> bool:
    try:
        nodes1 = sorted(g1.nodes())
        nodes2 = sorted(g2.nodes())
        
        adj1 = nx.to_numpy_array(g1, nodelist=nodes1, dtype=float, weight=None)
        adj2 = nx.to_numpy_array(g2, nodelist=nodes2, dtype=float, weight=None)
        
        eig1 = np.sort(np.round(np.linalg.eigvals(adj1), 5))
        eig2 = np.sort(np.round(np.linalg.eigvals(adj2), 5))
        
        _, s1, _ = np.linalg.svd(adj1)
        _, s2, _ = np.linalg.svd(adj2)
        
        return np.allclose(eig1, eig2, atol=1e-5) and np.allclose(s1, s2, atol=1e-5)
    except np.linalg.LinAlgError:
        return False

def construir_grafo(edges, isolated_vertices=[]):
    g = Grafo()
    for v in isolated_vertices:
        g.adicionar_aresta(v, v) 
        g.remover_aresta(v, v)   
    for u, v in edges:
        g.adicionar_aresta(u, v)
    return g

test_cases = [
    {
        "desc": "Grafos isomorfos simples (triângulos)",
        "g1_edges": [(1,2), (2,3), (3,1)],
        "g1_isolated": [],
        "g2_edges": [(4,5), (5,6), (6,4)],
        "g2_isolated": [],
        "esperado": True
    },
    {
        "desc": "Grafos não isomorfos (cadeia vs estrela)",
        "g1_edges": [(1,2), (2,3), (3,4)],
        "g1_isolated": [],
        "g2_edges": [(1,2), (1,3), (1,4)],
        "g2_isolated": [],
        "esperado": False
    },
    {
        "desc": "Grafos com vértices isolados",
        "g1_edges": [(1,2), (2,3)],
        "g1_isolated": [4],
        "g2_edges": [(5,6), (6,7)],
        "g2_isolated": [8],
        "esperado": True
    },
    {
        "desc": "Grafos completos K3",
        "g1_edges": [(1,2), (1,3), (2,3)],
        "g1_isolated": [],
        "g2_edges": [(4,5), (4,6), (5,6)],
        "g2_isolated": [],
        "esperado": True
    },
    {
        "desc": "Grafos vazios",
        "g1_edges": [],
        "g1_isolated": [],
        "g2_edges": [],
        "g2_isolated": [],
        "esperado": True
    },
    {
        "desc": "Grafos complexos isomorfos",
        "g1_edges": [(1,2), (1,5), (2,3), (2,4), (3,4), (4,5)],
        "g1_isolated": [],
        "g2_edges": [(5,4), (1,5), (3,4), (2,4), (2,3), (1,2)],
        "g2_isolated": [],
        "esperado": True
    }
]

# Execução dos testes
for case in test_cases:
    g1 = construir_grafo(case["g1_edges"], case["g1_isolated"])
    g2 = construir_grafo(case["g2_edges"], case["g2_isolated"])
    
    resultado = sao_isomorficos(g1, g2)
    status = "SUCESSO" if resultado == case["esperado"] else "FALHA"
    
    print(f"\nTeste: {case['desc']}")
    print(f"Vértices G1: {g1.vertices()}")
    print(f"Arestas G1: {g1.arestas()}")
    print(f"Vértices G2: {g2.vertices()}")
    print(f"Arestas G2: {g2.arestas()}")
    print(f"Resultado esperado: {case['esperado']}")
    print(f"Resultado obtido: {resultado}")
    print(f"Status: {status}")

print("\nTestes concluídos!")