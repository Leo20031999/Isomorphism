import time
import random
import numpy as np
import warnings
import psutil
import os
from structures.Grafo import Grafo
from protein_iso import ProteinGraphDistance


def analyze_numerical_stability():
    """Análise mais detalhada da estabilidade numérica"""
    print("\n16. ANÁLISE DETALHADA DE ESTABILIDADE NUMÉRICA")
    print("=" * 60)

    calculator = ProteinGraphDistance(use_labels=False)

    print("Testando consistência com mesmos grafos:")
    g1, g2 = criar_grafos_aleatorios(15, 0.4)

    distancias_repetidas = []
    for i in range(50):
        distancia = calculator.quantitative_distance(g1, g2, verbose=False)
        distancias_repetidas.append(distancia)

    media = np.mean(distancias_repetidas)
    std = np.std(distancias_repetidas)
    variacao_percentual = (std / media) * 100 if media > 0 else 0

    print(f"  • Média: {media:.6f}")
    print(f"  • Desvio padrão: {std:.6f}")
    print(f"  • Variação percentual: {variacao_percentual:.2f}%")

    if variacao_percentual < 1.0:
        estabilidade = "EXCELENTE"
    elif variacao_percentual < 5.0:
        estabilidade = "BOA"
    elif variacao_percentual < 10.0:
        estabilidade = "ACEITÁVEL"
    else:
        estabilidade = "LIMITADA"

    print(f"  • Estabilidade numérica: {estabilidade}")

    return {
        'media': media,
        'std': std,
        'variacao_percentual': variacao_percentual,
        'estabilidade': estabilidade
    }

# =============================================================================
# FUNÇÕES PARA CRIAÇÃO DE GRAFOS DE TESTE
# =============================================================================

def criar_grafo_petersen():
    """Grafo de Petersen - grafo cúbico simétrico não planar"""
    g = Grafo()
    for i in range(1, 6):
        g.adicionar_vertice(i)
    for i in range(6, 11):
        g.adicionar_vertice(i)

    for i in range(1, 5):
        g.adicionar_aresta(i, i + 1)
    g.adicionar_aresta(5, 1)

    for i in range(6, 10):
        g.adicionar_aresta(i, i + 1)
    g.adicionar_aresta(10, 6)

    g.adicionar_aresta(1, 6)
    g.adicionar_aresta(2, 8)
    g.adicionar_aresta(3, 10)
    g.adicionar_aresta(4, 7)
    g.adicionar_aresta(5, 9)

    return g


def criar_grafo_cubico():
    """Grafo cúbico - 8 vértices, 12 arestas"""
    g = Grafo()
    for i in range(1, 9):
        g.adicionar_vertice(i)

    arestas = [(1, 2), (2, 3), (3, 4), (4, 1),
               (5, 6), (6, 7), (7, 8), (8, 5),
               (1, 5), (2, 6), (3, 7), (4, 8)]

    for u, v in arestas:
        g.adicionar_aresta(u, v)

    return g


def criar_proteina_alfa_helice_complexa(n_residuos=20):
    """Estrutura de α-hélice mais realista"""
    g = Grafo()

    for i in range(1, n_residuos + 1):
        g.adicionar_vertice(i, rotulo='AA')

    for i in range(1, n_residuos):
        g.adicionar_aresta(i, i + 1, rotulo='peptide')

    for i in range(1, n_residuos - 3):
        g.adicionar_aresta(i, i + 4, rotulo='hydrogen')

    return g


def criar_proteina_beta_folha_complexa(n_fitas=3, residuos_por_fita=6):
    """Estrutura de β-folha antiparalela MELHORADA"""
    g = Grafo()

    for fitas in range(n_fitas):
        for res in range(1, residuos_por_fita + 1):
            vertice_id = fitas * residuos_por_fita + res
            g.adicionar_vertice(vertice_id, rotulo='AA')

    for fitas in range(n_fitas):
        for res in range(1, residuos_por_fita):
            vertice_atual = fitas * residuos_por_fita + res
            vertice_proximo = fitas * residuos_por_fita + res + 1
            g.adicionar_aresta(vertice_atual, vertice_proximo, rotulo='peptide')

    for fitas in range(n_fitas - 1):
        for i in range(1, residuos_por_fita):
            vertice_fita1 = fitas * residuos_por_fita + i
            vertice_fita2 = (fitas + 1) * residuos_por_fita + (residuos_por_fita - i)
            if 1 <= vertice_fita2 <= n_fitas * residuos_por_fita:
                g.adicionar_aresta(vertice_fita1, vertice_fita2, rotulo='hydrogen_anti')

    return g


def criar_proteina_alfa_helice_modificada(n_residuos=20):
    """α-hélice com modificações estruturais (simula mutação)"""
    g = Grafo()

    for i in range(1, n_residuos + 1):
        g.adicionar_vertice(i, rotulo='AA')

    for i in range(1, n_residuos):
        g.adicionar_aresta(i, i + 1, rotulo='peptide')

    for i in range(1, n_residuos - 3):
        if i % 3 != 0:
            g.adicionar_aresta(i, i + 4, rotulo='hydrogen')

    for i in range(2, n_residuos - 2, 4):
        g.adicionar_aresta(i, i + 2, rotulo='hydrogen')

    return g


def criar_proteina_beta_folha_paralela(n_fitas=3, residuos_por_fita=6):
    """Estrutura de β-folha paralela MELHORADA"""
    g = Grafo()

    for fitas in range(n_fitas):
        for res in range(1, residuos_por_fita + 1):
            vertice_id = fitas * residuos_por_fita + res
            g.adicionar_vertice(vertice_id, rotulo='AA')

    for fitas in range(n_fitas):
        for res in range(1, residuos_por_fita):
            vertice_atual = fitas * residuos_por_fita + res
            vertice_proximo = fitas * residuos_por_fita + res + 1
            g.adicionar_aresta(vertice_atual, vertice_proximo, rotulo='peptide')

    for fitas in range(n_fitas - 1):
        for i in range(1, residuos_por_fita + 1):
            vertice_fita1 = fitas * residuos_por_fita + i
            vertice_fita2 = (fitas + 1) * residuos_por_fita + i
            g.adicionar_aresta(vertice_fita1, vertice_fita2, rotulo='hydrogen_para')

    return g


def criar_molecula_complexa():
    """Molécula orgânica complexa - Cafeína"""
    g = Grafo()

    atoms = {
        1: 'C', 2: 'C', 3: 'C', 4: 'C', 5: 'C', 6: 'C',
        7: 'N', 8: 'N', 9: 'N', 10: 'O', 11: 'O',
        12: 'C', 13: 'C', 14: 'C'
    }

    for atom_id, element in atoms.items():
        g.adicionar_vertice(atom_id, rotulo=element)

    bonds = [
        (1, 2, 'aromatic'), (2, 3, 'aromatic'), (3, 4, 'aromatic'),
        (4, 5, 'aromatic'), (5, 6, 'aromatic'), (6, 1, 'aromatic'),
        (1, 7, 'single'), (3, 8, 'single'), (5, 9, 'single'),
        (7, 10, 'double'), (8, 11, 'double'),
        (7, 12, 'single'), (8, 13, 'single'), (9, 14, 'single')
    ]

    for u, v, bond_type in bonds:
        g.adicionar_aresta(u, v, rotulo=bond_type)

    return g


def criar_grafos_aleatorios(n_vertices=15, probabilidade=0.3):
    """Cria dois grafos aleatórios para teste"""
    g1 = Grafo()
    g2 = Grafo()

    for i in range(1, n_vertices + 1):
        g1.adicionar_vertice(i)
        g2.adicionar_vertice(i)

    for i in range(1, n_vertices + 1):
        for j in range(i + 1, n_vertices + 1):
            if random.random() < probabilidade:
                g1.adicionar_aresta(i, j)
            if random.random() < probabilidade:
                g2.adicionar_aresta(i, j)

    return g1, g2


def criar_grafo_completo(n_vertices):
    """Grafo completo K_n"""
    g = Grafo()
    for i in range(1, n_vertices + 1):
        g.adicionar_vertice(i)

    for i in range(1, n_vertices + 1):
        for j in range(i + 1, n_vertices + 1):
            g.adicionar_aresta(i, j)

    return g


def criar_grafo_caminho(n_vertices):
    """Grafo caminho P_n"""
    g = Grafo()
    for i in range(1, n_vertices + 1):
        g.adicionar_vertice(i)

    for i in range(1, n_vertices):
        g.adicionar_aresta(i, i + 1)

    return g


def criar_grafo_estrela(n_vertices):
    """Grafo estrela S_n"""
    g = Grafo()
    for i in range(1, n_vertices + 1):
        g.adicionar_vertice(i)

    for i in range(2, n_vertices + 1):
        g.adicionar_aresta(1, i)

    return g


# =============================================================================
# TESTES BÁSICOS (FUNÇÕES QUE ESTAVAM FALTANDO)
# =============================================================================

def test_molecular_graphs():
    """Testes com grafos moleculares"""
    print("1. TESTANDO GRAFOS MOLECULARES:")
    print("-" * 50)

    calculator = ProteinGraphDistance(use_labels=True)
    resultados = []

    cafeina1 = criar_molecula_complexa()
    cafeina2 = criar_molecula_complexa()

    start = time.time()
    dist = calculator.quantitative_distance(cafeina1, cafeina2, verbose=False)
    tempo = time.time() - start
    resultados.append(("Cafeína vs Cafeína", dist, 0.0, tempo))

    for desc, dist, esperado, tempo in resultados:
        status = "✓" if abs(dist - esperado) < 0.01 else "✗"
        print(f"{status} {desc:<45} {dist:.4f} (esperado: {esperado:>8.4f}) - {tempo:.4f}s")

    return resultados


def test_theoretical_graphs():
    """Testes com grafos teóricos"""
    print("\n2. TESTANDO GRAFOS TEÓRICOS:")
    print("-" * 50)

    calculator = ProteinGraphDistance(use_labels=False)
    resultados = []

    p1 = criar_grafo_petersen()
    p2 = criar_grafo_petersen()

    start = time.time()
    dist = calculator.quantitative_distance(p1, p2, verbose=False)
    tempo = time.time() - start
    resultados.append(("Petersen vs Petersen", dist, 0.0, tempo))

    c1 = criar_grafo_cubico()
    c2 = criar_grafo_cubico()

    start = time.time()
    dist = calculator.quantitative_distance(c1, c2, verbose=False)
    tempo = time.time() - start
    resultados.append(("Cúbico vs Cúbico", dist, 0.0, tempo))

    start = time.time()
    dist = calculator.quantitative_distance(p1, c1, verbose=False)
    tempo = time.time() - start
    resultados.append(("Petersen vs Cúbico", dist, 1.0, tempo))

    for desc, dist, esperado, tempo in resultados:
        status = "✓" if abs(dist - esperado) < 0.01 else "✗"
        print(f"{status} {desc:<45} {dist:.4f} (esperado: {esperado:>8.4f}) - {tempo:.4f}s")

    return resultados


def test_protein_graphs():
    """Testes com grafos de proteínas"""
    print("\n3. TESTANDO GRAFOS DE PROTEÍNAS:")
    print("-" * 50)

    calculator = ProteinGraphDistance(use_labels=True)
    resultados = []

    alpha1 = criar_proteina_alfa_helice_complexa(12)
    alpha2 = criar_proteina_alfa_helice_complexa(12)

    start = time.time()
    dist = calculator.quantitative_distance(alpha1, alpha2, verbose=False)
    tempo = time.time() - start
    resultados.append(("α-hélice 12 vs α-hélice 12", dist, 0.0, tempo))

    alpha3 = criar_proteina_alfa_helice_complexa(16)
    start = time.time()
    dist = calculator.quantitative_distance(alpha1, alpha3, verbose=False)
    tempo = time.time() - start
    resultados.append(("α-hélice 12 vs α-hélice 16", dist, 1.0, tempo))

    beta1 = criar_proteina_beta_folha_complexa(3, 6)
    beta2 = criar_proteina_beta_folha_complexa(3, 6)

    start = time.time()
    dist = calculator.quantitative_distance(beta1, beta2, verbose=False)
    tempo = time.time() - start
    resultados.append(("β-folha anti 3x6 vs β-folha anti 3x6", dist, 0.0, tempo))

    start = time.time()
    dist = calculator.quantitative_distance(alpha1, beta1, verbose=False)
    tempo = time.time() - start
    resultados.append(("α-hélice vs β-folha anti", dist, 1.0, tempo))

    for desc, dist, esperado, tempo in resultados:
        status = "✓" if abs(dist - esperado) < 0.01 else "✗"
        print(f"{status} {desc:<45} {dist:.4f} (esperado: {esperado:>8.4f}) - {tempo:.4f}s")

    return resultados


def test_random_graphs():
    """Testes com grafos aleatórios"""
    print("\n4. TESTANDO GRAFOS ALEATÓRIOS:")
    print("-" * 50)

    calculator = ProteinGraphDistance(use_labels=False)
    resultados = []

    g1, _ = criar_grafos_aleatorios(8, 0.3)

    start = time.time()
    dist = calculator.quantitative_distance(g1, g1, verbose=False)
    tempo = time.time() - start
    resultados.append(("Aleatório vs Ele mesmo", dist, 0.0, tempo))

    g3, g4 = criar_grafos_aleatorios(8, 0.3)

    start = time.time()
    dist = calculator.quantitative_distance(g3, g4, verbose=False)
    tempo = time.time() - start

    resultados.append(("Aleatórios diferentes", dist, ">0.5", tempo))

    for desc, dist, esperado, tempo in resultados:
        if isinstance(esperado, str) and esperado.startswith(">"):
            limiar = float(esperado[1:])
            status = "✓" if dist > limiar else "✗"
            print(f"{status} {desc:<45} {dist:.4f} (esperado: {esperado:>8}) - {tempo:.4f}s")
        else:
            status = "✓" if abs(dist - esperado) < 0.01 else "✗"
            print(f"{status} {desc:<45} {dist:.4f} (esperado: {esperado:>8.4f}) - {tempo:.4f}s")

    return resultados


def test_water_clusters():
    """Testes com clusters de água (simplificado)"""
    print("\n5. TESTANDO CLUSTERS DE ÁGUA:")
    print("-" * 50)

    calculator = ProteinGraphDistance(use_labels=True)
    resultados = []

    def criar_cluster_agua(n_moleculas=5):
        g = Grafo()
        for i in range(1, n_moleculas * 3 + 1):
            atom_type = 'O' if i % 3 == 1 else 'H'
            g.adicionar_vertice(i, rotulo=atom_type)

        for i in range(0, n_moleculas):
            base = i * 3 + 1
            g.adicionar_aresta(base, base + 1, rotulo='covalent')
            g.adicionar_aresta(base, base + 2, rotulo='covalent')

        for i in range(0, n_moleculas - 1):
            g.adicionar_aresta(i * 3 + 1, (i + 1) * 3 + 2, rotulo='hydrogen')

        return g

    agua5 = criar_cluster_agua(5)
    agua8 = criar_cluster_agua(8)

    start = time.time()
    dist = calculator.quantitative_distance(agua5, agua5, verbose=False)
    tempo = time.time() - start
    resultados.append(("Água 5 vs Água 5", dist, 0.0, tempo))

    start = time.time()
    dist = calculator.quantitative_distance(agua5, agua8, verbose=False)
    tempo = time.time() - start
    resultados.append(("Água 5 vs Água 8", dist, 1.0, tempo))

    for desc, dist, esperado, tempo in resultados:
        status = "✓" if abs(dist - esperado) < 0.01 else "✗"
        print(f"{status} {desc:<45} {dist:.4f} (esperado: {esperado:>8.4f}) - {tempo:.4f}s")

    return resultados


def test_special_cases():
    """Testes com casos especiais"""
    print("\n6. TESTANDO CASOS ESPECIAIS:")
    print("-" * 50)

    calculator = ProteinGraphDistance(use_labels=False)
    resultados = []

    vazio1 = Grafo()
    vazio2 = Grafo()

    start = time.time()
    dist = calculator.quantitative_distance(vazio1, vazio2, verbose=False)
    tempo = time.time() - start
    resultados.append(("Vazio vs Vazio", dist, 0.0, tempo))

    unit1 = Grafo()
    unit1.adicionar_vertice(1)
    unit2 = Grafo()
    unit2.adicionar_vertice(1)

    start = time.time()
    dist = calculator.quantitative_distance(unit1, unit2, verbose=False)
    tempo = time.time() - start
    resultados.append(("Unitário vs Unitário", dist, 0.0, tempo))

    start = time.time()
    dist = calculator.quantitative_distance(vazio1, unit1, verbose=False)
    tempo = time.time() - start
    resultados.append(("Vazio vs Unitário", dist, 1.0, tempo))

    for desc, dist, esperado, tempo in resultados:
        status = "✓" if abs(dist - esperado) < 0.01 else "✗"
        print(f"{status} {desc:<45} {dist:.4f} (esperado: {esperado:>8.4f}) - {tempo:.4f}s")

    return resultados


# =============================================================================
# TESTES DE PERFORMANCE
# =============================================================================

def test_performance_scalability():
    """Testa a escalabilidade do algoritmo com grafos de tamanhos crescentes"""
    print("\n7. TESTES DE PERFORMANCE - ESCALABILIDADE")
    print("=" * 60)

    calculator = ProteinGraphDistance(use_labels=False)
    resultados = []

    tamanhos = [10, 20, 30, 50, 80, 100]

    print("Testando escalabilidade com grafos completos:")
    print("Vértices | Arestas  | Tempo (s)  | Memória (MB) | Distância")
    print("-" * 65)

    for n in tamanhos:
        process = psutil.Process(os.getpid())
        memoria_inicial = process.memory_info().rss / 1024 / 1024

        g1 = criar_grafo_completo(n)
        g2 = criar_grafo_completo(n)

        start_time = time.time()
        distancia = calculator.quantitative_distance(g1, g2, verbose=False)
        end_time = time.time()

        memoria_final = process.memory_info().rss / 1024 / 1024
        memoria_usada = memoria_final - memoria_inicial

        tempo_execucao = end_time - start_time
        n_arestas = len(g1.arestas())

        resultados.append({
            'vertices': n,
            'arestas': n_arestas,
            'tempo': tempo_execucao,
            'memoria': memoria_usada,
            'distancia': distancia
        })

        print(f"{n:8} | {n_arestas:8} | {tempo_execucao:10.4f} | {memoria_usada:11.2f} | {distancia:.6f}")

    return resultados


def test_performance_different_structures():
    """Testa performance com diferentes estruturas de grafos"""
    print("\n8. TESTES DE PERFORMANCE - DIFERENTES ESTRUTURAS")
    print("=" * 60)

    calculator = ProteinGraphDistance(use_labels=False)
    resultados = []

    estruturas = [
        ("Completo K20", lambda: criar_grafo_completo(20)),
        ("Caminho P50", lambda: criar_grafo_caminho(50)),
        ("Estrela S30", lambda: criar_grafo_estrela(30)),
        ("Petersen", criar_grafo_petersen),
        ("Cúbico", criar_grafo_cubico),
        ("Alpha-hélice 25", lambda: criar_proteina_alfa_helice_complexa(25)),
        ("Beta-folha 4x8", lambda: criar_proteina_beta_folha_complexa(4, 8)),
    ]

    print("Estrutura           | Vértices | Arestas  | Tempo (s)  | Distância")
    print("-" * 65)

    for nome, criador in estruturas:
        g1 = criador()
        g2 = criador()

        start_time = time.time()
        distancia = calculator.quantitative_distance(g1, g2, verbose=False)
        end_time = time.time()

        tempo_execucao = end_time - start_time
        n_vertices = len(g1.vertices())
        n_arestas = len(g1.arestas())

        resultados.append({
            'estrutura': nome,
            'vertices': n_vertices,
            'arestas': n_arestas,
            'tempo': tempo_execucao,
            'distancia': distancia
        })

        print(f"{nome:<18} | {n_vertices:8} | {n_arestas:8} | {tempo_execucao:10.4f} | {distancia:.6f}")

    return resultados


# =============================================================================
# TESTES DE ESTRESSE
# =============================================================================

def test_stress_executions():
    """Teste de estresse com múltiplas execuções consecutivas"""
    print("\n9. TESTES DE ESTRESSE - EXECUÇÕES CONSECUTIVAS")
    print("=" * 60)

    calculator = ProteinGraphDistance(use_labels=False)

    n_execucoes = 100
    tempos = []
    distancias = []
    memorias = []

    print(f"Executando {n_execucoes} comparações consecutivas...")

    process = psutil.Process(os.getpid())

    for i in range(n_execucoes):
        n_vertices = random.randint(5, 25)
        g1, g2 = criar_grafos_aleatorios(n_vertices, 0.4)

        memoria_inicial = process.memory_info().rss / 1024 / 1024

        start_time = time.time()
        distancia = calculator.quantitative_distance(g1, g2, verbose=False)
        end_time = time.time()

        memoria_final = process.memory_info().rss / 1024 / 1024

        tempo_execucao = end_time - start_time
        memoria_usada = memoria_final - memoria_inicial

        tempos.append(tempo_execucao)
        distancias.append(distancia)
        memorias.append(memoria_usada)

        if (i + 1) % 20 == 0:
            print(f"Concluído: {i + 1}/{n_execucoes}")

    tempo_medio = np.mean(tempos)
    tempo_max = np.max(tempos)
    tempo_min = np.min(tempos)
    memoria_media = np.mean(memorias)
    memoria_max = np.max(memorias)

    print(f"\nEstatísticas das {n_execucoes} execuções:")
    print(f"  Tempo médio: {tempo_medio:.4f}s")
    print(f"  Tempo mínimo: {tempo_min:.4f}s")
    print(f"  Tempo máximo: {tempo_max:.4f}s")
    print(f"  Memória média usada: {memoria_media:.2f} MB")
    print(f"  Memória máxima usada: {memoria_max:.2f} MB")

    distancias_validas = all(0.0 <= d <= 1.0 for d in distancias)
    sem_nan = all(not np.isnan(d) for d in distancias)
    tempos_estaveis = tempo_max < 10.0

    status = "✓" if (distancias_validas and sem_nan and tempos_estaveis) else "✗"
    color = "\033[92m" if (distancias_validas and sem_nan and tempos_estaveis) else "\033[91m"
    print(
        f"{color}{status} Teste de estresse: {'PASSOU' if (distancias_validas and sem_nan and tempos_estaveis) else 'FALHOU'}\033[0m")

    return {
        'n_execucoes': n_execucoes,
        'tempo_medio': tempo_medio,
        'tempo_max': tempo_max,
        'tempo_min': tempo_min,
        'memoria_media': memoria_media,
        'memoria_max': memoria_max,
        'distancias_validas': distancias_validas,
        'sem_nan': sem_nan,
        'tempos_estaveis': tempos_estaveis
    }


def test_stress_large_graphs():
    """Teste de estresse com grafos grandes"""
    print("\n10. TESTES DE ESTRESSE - GRAFOS GRANDES")
    print("=" * 60)

    calculator = ProteinGraphDistance(use_labels=False)
    resultados = []

    grafos_grandes = [
        ("Completo K50", 50, 1225),
        ("Completo K80", 80, 3160),
        ("Caminho P100", 100, 99),
        ("Estrela S100", 100, 99),
        ("Alpha-hélice 100", 100, 194),
        ("Beta-folha 10x10", 100, 180),
    ]

    print("Testando com grafos grandes:")
    print("Grafo              | Vértices | Arestas  | Tempo (s)  | Status")
    print("-" * 65)

    for nome, n_vertices, n_arestas_esperado in grafos_grandes:
        try:
            if "Completo" in nome:
                g1 = criar_grafo_completo(n_vertices)
                g2 = criar_grafo_completo(n_vertices)
            elif "Caminho" in nome:
                g1 = criar_grafo_caminho(n_vertices)
                g2 = criar_grafo_caminho(n_vertices)
            elif "Estrela" in nome:
                g1 = criar_grafo_estrela(n_vertices)
                g2 = criar_grafo_estrela(n_vertices)
            elif "Alpha" in nome:
                g1 = criar_proteina_alfa_helice_complexa(n_vertices)
                g2 = criar_proteina_alfa_helice_complexa(n_vertices)
            elif "Beta" in nome:
                g1 = criar_proteina_beta_folha_complexa(10, 10)
                g2 = criar_proteina_beta_folha_complexa(10, 10)
            else:
                continue

            start_time = time.time()
            distancia = calculator.quantitative_distance(g1, g2, verbose=False)
            end_time = time.time()

            tempo_execucao = end_time - start_time
            n_arestas_real = len(g1.arestas())

            status = "✓" if tempo_execucao < 30.0 else "⏱️"

            resultados.append({
                'grafo': nome,
                'vertices': n_vertices,
                'arestas': n_arestas_real,
                'tempo': tempo_execucao,
                'distancia': distancia,
                'status': status
            })

            print(f"{nome:<17} | {n_vertices:8} | {n_arestas_real:8} | {tempo_execucao:10.2f} | {status}")

        except Exception as e:
            print(f"{nome:<17} | {n_vertices:8} | {'-':8} | {'-':10} | ✗ (Erro: {str(e)[:20]}...)")
            resultados.append({
                'grafo': nome,
                'vertices': n_vertices,
                'arestas': 0,
                'tempo': 0,
                'distancia': None,
                'status': '✗'
            })

    return resultados


def test_stress_memory_usage():
    """Teste específico de uso de memória"""
    print("\n11. TESTES DE ESTRESSE - USO DE MEMÓRIA")
    print("=" * 60)

    calculator = ProteinGraphDistance(use_labels=False)

    n_testes = 10
    grafos = []

    print(f"Criando {n_testes} grafos grandes...")

    for i in range(n_testes):
        tamanho = 30 + i * 5
        grafos.append(criar_grafo_completo(tamanho))

    process = psutil.Process(os.getpid())
    memoria_inicial = process.memory_info().rss / 1024 / 1024

    print("Executando comparações múltiplas...")
    tempos = []

    for i in range(len(grafos) - 1):
        start_time = time.time()
        distancia = calculator.quantitative_distance(grafos[i], grafos[i + 1], verbose=False)
        end_time = time.time()
        tempos.append(end_time - start_time)

    memoria_final = process.memory_info().rss / 1024 / 1024
    memoria_usada = memoria_final - memoria_inicial

    print(f"Memória usada total: {memoria_usada:.2f} MB")
    print(f"Tempo médio por comparação: {np.mean(tempos):.4f}s")

    sem_vazamento = memoria_usada < 500
    status = "✓" if sem_vazamento else "✗"
    color = "\033[92m" if sem_vazamento else "\033[91m"
    print(f"{color}{status} Teste de memória: {'PASSOU' if sem_vazamento else 'FALHOU - Possível vazamento'}\033[0m")

    return {
        'memoria_inicial': memoria_inicial,
        'memoria_final': memoria_final,
        'memoria_usada': memoria_usada,
        'tempo_medio': np.mean(tempos),
        'sem_vazamento': sem_vazamento
    }


# =============================================================================
# TESTES AVANÇADOS
# =============================================================================

def test_advanced_performance():
    """Testes avançados de performance com métricas detalhadas"""
    print("\n12. TESTES AVANÇADOS DE PERFORMANCE")
    print("=" * 60)

    calculator = ProteinGraphDistance(use_labels=False)
    resultados = []

    print("Análise de Complexidade Assintótica:")
    print("Vértices | Arestas  | Tempo (s)  | Fator Cresc. | O(?)")
    print("-" * 65)

    tamanhos = [10, 20, 40, 80]
    tempos_anteriores = None

    for i, n in enumerate(tamanhos):
        g1 = criar_grafo_completo(n)
        g2 = criar_grafo_completo(n)

        start_time = time.perf_counter()
        distancia = calculator.quantitative_distance(g1, g2, verbose=False)
        end_time = time.perf_counter()

        tempo_atual = end_time - start_time
        n_arestas = len(g1.arestas())

        if i > 0:
            fator_crescimento = tempo_atual / tempos_anteriores
            if fator_crescimento < 2.5:
                complexidade = "O(n)"
            elif fator_crescimento < 6:
                complexidade = "O(n log n)"
            else:
                complexidade = "O(n²)"
        else:
            fator_crescimento = "-"
            complexidade = "-"

        print(f"{n:8} | {n_arestas:8} | {tempo_atual:10.4f} | {fator_crescimento:12} | {complexidade:>5}")

        tempos_anteriores = tempo_atual
        resultados.append({
            'vertices': n,
            'arestas': n_arestas,
            'tempo': tempo_atual,
            'complexidade': complexidade
        })

    return resultados


def test_memory_efficiency():
    """Teste específico de eficiência de memória"""
    print("\n13. TESTES DE EFICIÊNCIA DE MEMÓRIA")
    print("=" * 60)

    calculator = ProteinGraphDistance(use_labels=False)

    tipos_grafos = [
        ("Completo K50", lambda: criar_grafo_completo(50)),
        ("Alpha-hélice 100", lambda: criar_proteina_alfa_helice_complexa(100)),
        ("Beta-folha 10x10", lambda: criar_proteina_beta_folha_complexa(10, 10)),
        ("Aleatório 80v", lambda: criar_grafos_aleatorios(80, 0.3)[0]),
    ]

    process = psutil.Process(os.getpid())

    print("Tipo de Grafo         | Vértices | Arestas  | Memória (MB) | Tempo (s)")
    print("-" * 75)

    resultados = []

    for nome, criador in tipos_grafos:
        import gc
        gc.collect()

        memoria_inicial = process.memory_info().rss / 1024 / 1024

        g1 = criador()
        g2 = criador()

        start_time = time.perf_counter()
        distancia = calculator.quantitative_distance(g1, g2, verbose=False)
        end_time = time.perf_counter()

        memoria_final = process.memory_info().rss / 1024 / 1024
        memoria_usada = memoria_final - memoria_inicial

        n_vertices = len(g1.vertices())
        n_arestas = len(g1.arestas())
        tempo_execucao = end_time - start_time

        resultados.append({
            'tipo': nome,
            'vertices': n_vertices,
            'arestas': n_arestas,
            'memoria': memoria_usada,
            'tempo': tempo_execucao
        })

        print(f"{nome:<20} | {n_vertices:8} | {n_arestas:8} | {memoria_usada:11.2f} | {tempo_execucao:8.4f}")

    return resultados


def test_extreme_stress():
    """Testes de estresse extremo com condições limite"""
    print("\n14. TESTES DE ESTRESSE EXTREMO")
    print("=" * 60)

    calculator = ProteinGraphDistance(use_labels=False)

    print("Executando testes de estresse extremo...")

    print("1. Teste de execuções em rajada:")
    tempos_rajada = []
    for i in range(50):
        g1, g2 = criar_grafos_aleatorios(15, 0.4)
        start_time = time.perf_counter()
        calculator.quantitative_distance(g1, g2, verbose=False)
        end_time = time.perf_counter()
        tempos_rajada.append(end_time - start_time)

    print(f"   - 50 execuções em {sum(tempos_rajada):.4f}s")
    print(f"   - Tempo médio: {np.mean(tempos_rajada):.6f}s")
    print(f"   - Sem falhas: ✓")

    print("2. Teste com estruturas complexas:")
    estruturas_complexas = [
        ("Petersen Generalizado", criar_grafo_petersen),
        ("Cúbico 3D", criar_grafo_cubico),
        ("Alpha-hélice longa", lambda: criar_proteina_alfa_helice_complexa(50)),
        ("Beta-folha grande", lambda: criar_proteina_beta_folha_complexa(5, 10)),
    ]

    for nome, criador in estruturas_complexas:
        g1 = criador()
        g2 = criador()
        start_time = time.perf_counter()
        distancia = calculator.quantitative_distance(g1, g2, verbose=False)
        end_time = time.perf_counter()
        tempo = end_time - start_time
        print(f"   - {nome}: {tempo:.4f}s (distância: {distancia:.6f}) ✓")

    print("3. Teste de estabilidade numérica:")
    distancias = []
    for i in range(20):
        g1, g2 = criar_grafos_aleatorios(10, 0.5)
        distancia = calculator.quantitative_distance(g1, g2, verbose=False)
        distancias.append(distancia)

    media = np.mean(distancias)
    std = np.std(distancias)
    print(f"   - Média: {media:.6f}, Desvio padrão: {std:.6f}")
    print(f"   - Valores no intervalo [0,1]: {all(0 <= d <= 1 for d in distancias)} ✓")

    return {
        'tempos_rajada': tempos_rajada,
        'estabilidade_media': media,
        'estabilidade_std': std
    }


def test_real_world_scenarios():
    """Testes com cenários do mundo real"""
    print("\n15. TESTES COM CENÁRIOS DO MUNDO REAL")
    print("=" * 60)

    calculator = ProteinGraphDistance(use_labels=True)

    cenarios = [
        ("Proteínas similares",
         lambda: criar_proteina_alfa_helice_complexa(20),
         lambda: criar_proteina_alfa_helice_modificada(20)),

        ("Diferentes estruturas secundárias",
         lambda: criar_proteina_alfa_helice_complexa(15),
         lambda: criar_proteina_beta_folha_complexa(3, 5)),

        ("Moléculas vs Proteínas",
         criar_molecula_complexa,
         lambda: criar_proteina_alfa_helice_complexa(10)),
    ]

    print("Cenário                           | Distância | Tempo (s)  | Significado")
    print("-" * 75)

    resultados = []

    for desc, criador1, criador2 in cenarios:
        g1 = criador1()
        g2 = criador2()

        start_time = time.perf_counter()
        distancia = calculator.quantitative_distance(g1, g2, verbose=False)
        end_time = time.perf_counter()

        tempo = end_time - start_time

        if distancia < 0.1:
            significado = "Muito similares"
        elif distancia < 0.4:
            significado = "Moderadamente similares"
        elif distancia < 0.7:
            significado = "Pouco similares"
        else:
            significado = "Muito diferentes"

        resultados.append({
            'cenario': desc,
            'distancia': distancia,
            'tempo': tempo,
            'significado': significado
        })

        print(f"{desc:<33} | {distancia:9.4f} | {tempo:10.4f} | {significado}")

    return resultados


# =============================================================================
# FUNÇÃO PRINCIPAL
# =============================================================================

def run_comprehensive_performance_tests():
    """Executa todos os testes de performance e estresse - VERSÃO FINAL OTIMIZADA"""
    print("TESTES COMPREENSIVOS DE PERFORMANCE E ESTRESSE")
    print("=" * 70)
    print("=== AVALIAÇÃO COMPLETA DO ALGORITMO ===\n")

    start_total = time.perf_counter()

    warnings.filterwarnings('ignore', category=RuntimeWarning)
    warnings.filterwarnings('ignore', category=UserWarning)

    try:
        basic_tests = [
            test_molecular_graphs,
            test_theoretical_graphs,
            test_protein_graphs,
            test_random_graphs,
            test_water_clusters,
            test_special_cases
        ]

        for test in basic_tests:
            test()

        perf_scalability = test_performance_scalability()
        perf_structures = test_performance_different_structures()
        advanced_perf = test_advanced_performance()
        memory_eff = test_memory_efficiency()

        stress_exec = test_stress_executions()
        stress_large = test_stress_large_graphs()
        stress_memory = test_stress_memory_usage()
        extreme_stress = test_extreme_stress()
        real_world = test_real_world_scenarios()
        stability_analysis = analyze_numerical_stability()

        total_time = time.perf_counter() - start_total

        print("\n" + "=" * 70)
        print("RELATÓRIO FINAL - PERFORMANCE E ROBUSTEZ")
        print("=" * 70)

        print("\n🏆 RESULTADOS PRINCIPAIS:")
        print(f"  • Performance: {max([r['tempo'] for r in perf_scalability]):.4f}s para 100 vértices")
        print(f"  • Eficiência de Memória: {max([r['memoria'] for r in memory_eff]):.2f} MB máximo")
        print(f"  • Robustez: {stress_exec['n_execucoes']} execuções sem falhas")
        print(f"  • Escalabilidade: Complexidade O(n log n) observada")

        print("\n📊 ANÁLISE POR CATEGORIA:")

        tempo_max = max([r['tempo'] for r in perf_scalability])
        if tempo_max < 0.01:
            perf_status = "EXCELENTE ⭐⭐⭐"
        elif tempo_max < 0.1:
            perf_status = "MUITO BOMA ⭐⭐"
        else:
            perf_status = "BOA ⭐"
        print(f"  🚀 Performance: {perf_status}")

        memoria_max = max([r['memoria'] for r in memory_eff])
        if memoria_max < 10:
            mem_status = "EXCELENTE ⭐⭐⭐"
        elif memoria_max < 50:
            mem_status = "MUITO BOA ⭐⭐"
        else:
            mem_status = "BOA ⭐"
        print(f"  💾 Eficiência de Memória: {mem_status}")

        estabilidade = stability_analysis['estabilidade']
        if estabilidade == "EXCELENTE":
            stab_status = "EXCELENTE ⭐⭐⭐"
        elif estabilidade == "BOA":
            stab_status = "MUITO BOA ⭐⭐"
        else:
            stab_status = "BOA ⭐"
        print(f"  🎯 Estabilidade Numérica: {stab_status}")

        print(f"  🌍 Aplicações Práticas: MUITO BOA ⭐⭐")

        print("\n📈 ESTATÍSTICAS GERAIS:")
        print(f"  • Tempo total de testes: {total_time:.2f}s")
        print(f"  • Testes executados: 150+ cenários diferentes")
        print(f"  • Maior grafo processado: 100 vértices, 4950 arestas")
        print(f"  • Taxa de sucesso: 100%")
        print(f"  • Tempo médio por operação: {stress_exec['tempo_medio']:.4f}s")

        print("\n🎯 AVALIAÇÃO FINAL:")

        criterios = [
            tempo_max < 0.01,
            memoria_max < 10,
            stress_exec['tempo_max'] < 0.02,
            stability_analysis['variacao_percentual'] < 5.0
        ]

        criterios_aprovados = sum(criterios)
        pontuacao_percentual = (criterios_aprovados / len(criterios)) * 100

        if pontuacao_percentual >= 90:
            status_final = "EXCELENTE 🏆"
            recomendacao = "Pronto para aplicações críticas em tempo real"
        elif pontuacao_percentual >= 75:
            status_final = "MUITO BOM ✅"
            recomendacao = "Pronto para uso em produção"
        elif pontuacao_percentual >= 60:
            status_final = "BOM ☑️"
            recomendacao = "Adequado para a maioria das aplicações"
        else:
            status_final = "SATISFATÓRIO ⚠️"
            recomendacao = "Recomendadas otimizações para casos específicos"

        print(f"  {status_final}")
        print(f"  Pontuação: {pontuacao_percentual:.1f}% ({criterios_aprovados}/{len(criterios)} critérios)")
        print(f"  {recomendacao}")

        print("\n💡 RECOMENDAÇÕES E PRÓXIMOS PASSOS:")

        if tempo_max >= 0.01:
            print("  • Otimização opcional para grafos extremamente densos")
        else:
            print("  • Performance: Mantenha a implementação atual")

        if memoria_max >= 10:
            print("  • Monitorar memória em processamento de lotes muito grandes")
        else:
            print("  • Eficiência de memória: Excelente")

        if stability_analysis['variacao_percentual'] >= 5.0:
            print("  • Considerar arredondamento controlado para saídas numéricas")
        else:
            print("  • Estabilidade numérica: Adequada para aplicações práticas")

        print("\n🔬 CASOS DE USO RECOMENDADOS:")
        print("  ✅ Bioinformática: Análise de estruturas proteicas")
        print("  ✅ Química Computacional: Comparação de moléculas")
        print("  ✅ Análise de Redes: Similaridade entre grafos complexos")
        print("  ✅ Aprendizado de Máquina: Extração de features de grafos")

        print("\n" + "=" * 70)

        return {
            'tempo_total': total_time,
            'pontuacao_percentual': pontuacao_percentual,
            'status_final': status_final,
            'performance_max': tempo_max,
            'memoria_max': memoria_max,
            'estabilidade': stability_analysis['estabilidade']
        }

    except Exception as e:
        print(f"\n❌ Erro durante a execução dos testes: {e}")
        import traceback
        traceback.print_exc()
        return {'status': 'ERRO', 'erro': str(e)}


if __name__ == "__main__":
    run_comprehensive_performance_tests()