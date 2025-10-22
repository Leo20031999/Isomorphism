# Isomorphism, Automorphism Partitioning, and Canonical Labeling Can Be Solved in Polynomial-Time for Molecular Graphs
# (JEAN-LOUP FAULON, 1998)
from structures.Grafo import Grafo
from collections import defaultdict
import networkx as nx

def transformar_para_grafo_simples(grafo_molecular, Z0=4):
    G_simples = Grafo()
    dummy_counter = 0

    def novo_dummy(tipo):
        nonlocal dummy_counter
        dummy_counter += 1
        return f"dummy_{tipo}_{dummy_counter}"

    for v in grafo_molecular.vertices():
        Z = grafo_molecular.get_rotulo_vertice(v)
        if Z is None:
            Z = 0

        G_simples.adicionar_vertice(v)

        for _ in range(Z + Z0):
            dummy = novo_dummy("atom")
            G_simples.adicionar_vertice(dummy)
            G_simples.adicionar_aresta(v, dummy)

    for (u, v) in grafo_molecular.arestas():
        ordem = grafo_molecular.get_rotulo_aresta(u, v) or 1

        if ordem == 1:
            G_simples.adicionar_aresta(u, v)
        else:
            prev = u

            for _ in range(1, ordem):
                dummy = novo_dummy("bond")
                G_simples.adicionar_vertice(dummy)
                G_simples.adicionar_aresta(prev, dummy)
                prev = dummy
            G_simples.adicionar_aresta(prev, v)

    return G_simples


def invariantes_equivalentes(G1, G2):
    """Verifica invariantes básicos entre grafos"""
    if len(G1.vertices()) != len(G2.vertices()):
        return False
    if len(G1.arestas()) != len(G2.arestas()):
        return False
    graus_G1 = sorted([G1.grau(v) for v in G1.vertices()])
    graus_G2 = sorted([G2.grau(v) for v in G2.vertices()])
    return graus_G1 == graus_G2


def calcular_invariantes(grafo):
    """Calcula invariantes para todos os vértices"""
    invariantes = {}
    for v in grafo.vertices():
        grau = grafo.grau(v)
        vizinhos = grafo.vizinhanca(v)
        graus_vizinhos = sorted([grafo.grau(n) for n in vizinhos])
        rotulo = grafo.get_rotulo_vertice(v) or 0
        invariante = (grau, tuple(graus_vizinhos), rotulo)
        invariantes[v] = hash(invariante)
    return invariantes


def isomorfismo_planar(G1, G2):
    """Verifica isomorfismo para grafos planares"""
    inv1 = calcular_invariantes(G1)
    inv2 = calcular_invariantes(G2)

    if sorted(inv1.values()) != sorted(inv2.values()):
        return False

    for v1 in G1.vertices():
        match = False
        for v2 in G2.vertices():
            if inv1[v1] != inv2[v2]:
                continue

            vizinhos1 = sorted(G1.vizinhanca(v1), key=lambda x: inv1.get(x, 0))
            vizinhos2 = sorted(G2.vizinhanca(v2), key=lambda x: inv2.get(x, 0))

            if len(vizinhos1) != len(vizinhos2):
                continue

            viz_match = True
            for i in range(len(vizinhos1)):
                if inv1.get(vizinhos1[i]) != inv2.get(vizinhos2[i]):
                    viz_match = False
                    break

            if viz_match:
                match = True
                break

        if not match:
            return False

    return True

def refinamento_cores(grafo, iteracoes=3):
    """Algoritmo de refinamento iterativo de cores"""
    cores = {}
    for v in grafo.vertices():
        grau = grafo.grau(v)
        rotulo = grafo.get_rotulo_vertice(v) or 0
        cores[v] = (rotulo, grau)

    for _ in range(iteracoes):
        novas_cores = {}
        mapa_cores = {}
        next_id = 0

        for v in grafo.vertices():
            cores_vizinhos = tuple(sorted(cores[u] for u in grafo.vizinhanca(v)))
            nova_cor = (cores[v], cores_vizinhos)

            if nova_cor not in mapa_cores:
                mapa_cores[nova_cor] = next_id
                next_id += 1

            novas_cores[v] = mapa_cores[nova_cor]

        cores = novas_cores

    return cores

def isomorfismo_valencia_limitada(G1, G2):
    """Algoritmo para grafos não planares"""
    cores1 = refinamento_cores(G1)
    cores2 = refinamento_cores(G2)

    contador1 = defaultdict(int)
    for cor in cores1.values():
        contador1[cor] += 1

    contador2 = defaultdict(int)
    for cor in cores2.values():
        contador2[cor] += 1

    return contador1 == contador2

def isomorfismo_molecular(G1, G2, Z0=4):
    """Verifica se dois grafos moleculares são isomorfos"""
    G1_simples = transformar_para_grafo_simples(G1, Z0)
    G2_simples = transformar_para_grafo_simples(G2, Z0)

    if not invariantes_equivalentes(G1_simples, G2_simples):
        return False

    try:
        planar1 = nx.is_planar(G1_simples.grafo)
        planar2 = nx.is_planar(G2_simples.grafo)
    except:
        planar1 = planar2 = False

    if planar1 and planar2:
        return isomorfismo_planar(G1_simples, G2_simples)
    else:
        return isomorfismo_valencia_limitada(G1_simples, G2_simples)

def automorfismos_moleculares(grafo_molecular, Z0=4):
    """Computa o grupo de automorfismo do grafo molecular"""
    G_simples = transformar_para_grafo_simples(grafo_molecular, Z0)
    invariantes = calcular_invariantes(G_simples)

    grupos = defaultdict(list)
    for v, inv in invariantes.items():
        grupos[inv].append(v)

    vertices_originais = grafo_molecular.vertices()
    automorfismos = []
    mapeamento_atual = {}
    mapeamento_inverso = {}

    def backtrack(idx):
        """Busca com retrocesso para encontrar automorfismos"""
        if idx == len(vertices_originais):
            valido = True
            for u, v in grafo_molecular.arestas():
                u_map = mapeamento_atual[u]
                v_map = mapeamento_atual[v]
                if not G_simples.grafo.has_edge(u_map, v_map):
                    valido = False
                    break
            if valido:
                automorfismos.append(mapeamento_atual.copy())
            return

        v = vertices_originais[idx]
        grupo = grupos[invariantes[v]]

        for candidato in grupo:
            if candidato in mapeamento_inverso:
                continue

            consistente = True
            for vizinho in grafo_molecular.vizinhanca(v):
                if vizinho in mapeamento_atual:
                    if not G_simples.grafo.has_edge(candidato, mapeamento_atual[vizinho]):
                        consistente = False
                        break

            if not consistente:
                continue

            mapeamento_atual[v] = candidato
            mapeamento_inverso[candidato] = v

            backtrack(idx + 1)

            del mapeamento_atual[v]
            del mapeamento_inverso[candidato]

    backtrack(0)
    return automorfismos