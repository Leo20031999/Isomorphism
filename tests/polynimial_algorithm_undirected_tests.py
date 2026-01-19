import gc
import time
import random
import psutil
import os
import math
from algorithms.polynomial_algorithm_undirected import are_isomorphic, _generate_isomorphic_group
from structures.Grafo import Grafo


def criar_grafo_petersen():
    """Grafo de Petersen - grafo cúbico simétrico não planar"""
    g = Grafo()
    for i in range(1, 11):
        g.adicionar_vertice(i)

    for i in range(1, 5):
        g.adicionar_aresta(i, i + 1)
    g.adicionar_aresta(5, 1)

    for i in range(6, 10):
        g.adicionar_aresta(i, i + 1)
    g.adicionar_aresta(10, 6)

    for i in range(1, 6):
        g.adicionar_aresta(i, i + 5)

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


def criar_grafos_aleatorios(n_vertices=15, probabilidade=0.3, seed=None):
    """Cria dois grafos aleatórios para teste - CORRIGIDA"""
    if seed is not None:
        random.seed(seed)

    g1 = Grafo()
    g2 = Grafo()

    for i in range(1, n_vertices + 1):
        g1.adicionar_vertice(i)
        g2.adicionar_vertice(i)

    random_values = []
    for i in range(1, n_vertices + 1):
        for j in range(i + 1, n_vertices + 1):
            random_values.append(random.random())

    value_index = 0
    for i in range(1, n_vertices + 1):
        for j in range(i + 1, n_vertices + 1):
            r = random_values[value_index]
            if r < probabilidade:
                g1.adicionar_aresta(i, j)
            if r < probabilidade:
                g2.adicionar_aresta(i, j)
            value_index += 1

    return g1, g2


def build_methane():
    """Constrói grafo do metano"""
    g = Grafo()
    g.adicionar_aresta("C", "H1")
    g.adicionar_aresta("C", "H2")
    g.adicionar_aresta("C", "H3")
    g.adicionar_aresta("C", "H4")
    return g


def build_ethane():
    """Constrói grafo do etano"""
    g = Grafo()
    g.adicionar_aresta("C1", "C2")
    for i in range(1, 4):
        g.adicionar_aresta("C1", f"H{i}")
    for i in range(4, 7):
        g.adicionar_aresta("C2", f"H{i}")
    return g


def build_butane():
    """Constrói grafo do butano"""
    g = Grafo()
    g.adicionar_aresta("C1", "C2")
    g.adicionar_aresta("C2", "C3")
    g.adicionar_aresta("C3", "C4")

    for i in range(1, 4):
        g.adicionar_aresta("C1", f"H1{i}")
    for i in range(1, 3):
        g.adicionar_aresta("C2", f"H2{i}")
    for i in range(1, 3):
        g.adicionar_aresta("C3", f"H3{i}")
    for i in range(1, 4):
        g.adicionar_aresta("C4", f"H4{i}")

    return g


def build_benzene():
    """Constrói grafo do benzeno (anel aromático)"""
    g = Grafo()
    for i in range(1, 7):
        g.adicionar_aresta(f"C{i}", f"C{(i % 6) + 1}")
    for i in range(1, 7):
        g.adicionar_aresta(f"C{i}", f"H{i}")
    return g


def build_water_cluster(n=5):
    """Constrói cluster de água (grafo mais complexo)"""
    g = Grafo()
    for i in range(n):
        g.adicionar_aresta(f"O{i}", f"H{i}_1")
        g.adicionar_aresta(f"O{i}", f"H{i}_2")

    for i in range(n - 1):
        g.adicionar_aresta(f"H{i}_1", f"O{i + 1}")

    return g


def build_cycle_graph(n):
    """Constrói grafo ciclo C_n"""
    g = Grafo()
    for i in range(n):
        g.adicionar_vertice(i)
    for i in range(n):
        g.adicionar_aresta(i, (i + 1) % n)
    return g


def build_path_graph(n):
    """Constrói grafo caminho P_n"""
    g = Grafo()
    for i in range(n):
        g.adicionar_vertice(i)
    for i in range(n - 1):
        g.adicionar_aresta(i, i + 1)
    return g


def build_complete_bipartite_graph(n, m):
    """Constrói grafo bipartido completo K_{n,m}"""
    g = Grafo()
    for i in range(n + m):
        g.adicionar_vertice(i)
    for i in range(n):
        for j in range(n, n + m):
            g.adicionar_aresta(i, j)
    return g


def build_complete_graph(n):
    """Constrói grafo completo K_n"""
    g = Grafo()
    for i in range(n):
        g.adicionar_vertice(i)
    for i in range(n):
        for j in range(i + 1, n):
            g.adicionar_aresta(i, j)
    return g


def build_petersen_graph():
    """Constrói grafo de Petersen (famoso contraexemplo)"""
    g = Grafo()
    for i in range(1, 11):
        g.adicionar_vertice(i)

    for i in range(1, 5):
        g.adicionar_aresta(i, i + 1)
    g.adicionar_aresta(5, 1)

    for i in range(6, 10):
        g.adicionar_aresta(i, i + 1)
    g.adicionar_aresta(10, 6)

    for i in range(1, 6):
        g.adicionar_aresta(i, i + 5)

    return g

def build_wheel_graph(n):
    """Constrói grafo roda W_n (ciclo + vértice central)"""
    g = Grafo()
    g.adicionar_vertice(0)
    for i in range(1, n):
        g.adicionar_vertice(i)
    for i in range(1, n):
        g.adicionar_aresta(i, (i % (n - 1)) + 1)
    for i in range(1, n):
        g.adicionar_aresta(0, i)
    return g


def build_star_graph(n):
    """Constrói grafo estrela S_n"""
    g = Grafo()
    g.adicionar_vertice(0)
    for i in range(1, n):
        g.adicionar_vertice(i)
        g.adicionar_aresta(0, i)
    return g


def build_small_non_isomorphic_regular_graphs():
    """Constrói pares de grafos regulares não isomórficos conhecidos"""

    g1 = build_cycle_graph(6)

    g2 = Grafo()
    for i in range(3):
        g2.adicionar_aresta(i, (i + 1) % 3)
    for i in range(3, 6):
        g2.adicionar_aresta(i, 3 + (i - 2) % 3)

    g3 = build_complete_bipartite_graph(3, 3)

    g4 = build_petersen_graph()

    return [("C6 vs 2xC3", g1, g2, False),
            ("K3,3 vs Petersen", g3, g4, False)]


def get_cpu_metrics():
    """Obtém métricas de CPU mais precisas"""
    try:
        process = psutil.Process(os.getpid())

        cpu_percent = psutil.cpu_percent(interval=0.1)
        cpu_times = psutil.cpu_times_percent(interval=0.1)
        cpu_freq = psutil.cpu_freq()
        load_avg = os.getloadavg() if hasattr(os, 'getloadavg') else (0, 0, 0)

        return {
            'system_cpu_percent': cpu_percent,
            'system_user_percent': cpu_times.user,
            'system_system_percent': cpu_times.system,
            'cpu_frequency': cpu_freq.current if cpu_freq else None,
            'cpu_cores': psutil.cpu_count(logical=False),
            'cpu_threads': psutil.cpu_count(logical=True),
            'load_avg_1min': load_avg[0],
            'load_avg_5min': load_avg[1],
            'load_avg_15min': load_avg[2]
        }
    except Exception as e:
        return {'error': f"Erro ao obter métricas CPU: {str(e)}"}


def measure_cpu_usage(func, *args, **kwargs):
    """Mede o uso de CPU de forma precisa e confiável"""
    process = psutil.Process(os.getpid())

    start_cpu_times = process.cpu_times()
    start_perf = time.perf_counter()
    start_process = time.process_time()

    result = func(*args, **kwargs)

    end_perf = time.perf_counter()
    end_process = time.process_time()
    end_cpu_times = process.cpu_times()

    real_time = end_perf - start_perf
    process_cpu_time = end_process - start_process

    user_time = end_cpu_times.user - start_cpu_times.user
    system_time = end_cpu_times.system - start_cpu_times.system
    total_cpu_time = user_time + system_time

    if real_time <= 0:
        real_time = 0.000001

    cpu_percent = (total_cpu_time / real_time) * 100
    efficiency = (process_cpu_time / real_time) * 100

    cpu_percent = min(cpu_percent, 1000)  
    efficiency = min(efficiency, 1000)

    return result, {
        'execution_time': real_time,
        'cpu_time_user': user_time,
        'cpu_time_system': system_time,
        'cpu_time_used': total_cpu_time,
        'process_cpu_time': process_cpu_time,
        'cpu_percent': cpu_percent,
        'efficiency': efficiency
    }


def run_cpu_analysis():
    """Executa análise de CPU com métricas corrigidas e estrutura compatível"""
    print("\n" + "=" * 100)
    print("ANÁLISE DE CPU - MÉTRICAS CORRIGIDAS")
    print("=" * 100)

    def test_with_larger_graph():
        g1 = build_complete_graph(20)
        g2 = build_complete_graph(20)
        return are_isomorphic(g1, g2)

    print("Executando teste com gráfico maior para medições precisas...")

    all_metrics = []
    for i in range(5):
        result, metrics = measure_cpu_usage(test_with_larger_graph)
        all_metrics.append(metrics)
        time.sleep(0.1)

    avg_metrics = {}
    for key in all_metrics[0].keys():
        values = [m[key] for m in all_metrics]
        avg_metrics[key] = sum(values) / len(values)

    print_cpu_metrics(avg_metrics, "Média de 5 execuções (K20)")

    system_metrics = get_cpu_metrics()
    print(f"\n💻 MÉTRICAS DO SISTEMA:")
    print(f"  • Núcleos físicos: {system_metrics['cpu_cores']}")
    print(f"  • Threads: {system_metrics['cpu_threads']}")
    print(f"  • Uso geral do CPU: {system_metrics['system_cpu_percent']}%")
    if system_metrics['cpu_frequency']:
        print(f"  • Frequência CPU: {system_metrics['cpu_frequency']} MHz")

    return {
        'total_time': avg_metrics['execution_time'],
        'total_cpu_user': avg_metrics['cpu_time_user'],
        'total_cpu_system': avg_metrics['cpu_time_system'],
        'cpu_metrics': {
            'avg_cpu_percent': avg_metrics['cpu_percent'],
            'efficiency': avg_metrics['efficiency']
        }
    }

def print_cpu_metrics(metrics, test_name=""):
    """Imprime métricas de CPU de forma mais precisa"""
    print(f"\n📊 MÉTRICAS DE CPU - {test_name}:")
    print(f"  ⏱️  Tempo de execução real: {metrics['execution_time']:.6f}s")
    print(f"  🔄 Tempo de CPU (usuário): {metrics['cpu_time_user']:.6f}s")
    print(f"  🔄 Tempo de CPU (sistema): {metrics['cpu_time_system']:.6f}s")
    print(f"  🔄 Tempo total de CPU: {metrics['cpu_time_used']:.6f}s")
    print(f"  📈 Uso de CPU: {metrics['cpu_percent']:.2f}%")
    print(f"  🎯 Eficiência: {metrics['efficiency']:.1f}%")

    if metrics['execution_time'] < 0.001:
        print("  ⚠️  Medição muito rápida - pode ser imprecisa")
    elif metrics['cpu_time_used'] < 0.0001:
        print("  ⚠️  Tempo de CPU muito baixo - possível imprecisão")


def test_basic_isomorphism():
    """Testes básicos de isomorfismo"""
    print("=== TESTES BÁSICOS DE ISOMORFISMO ===\n")

    methane1 = build_methane()
    methane2 = build_methane()

    ethane = build_ethane()
    butane = build_butane()

    benzene1 = build_benzene()
    benzene2 = build_benzene()

    tests = [
        ("Metano vs Metano", methane1, methane2, True),
        ("Metano vs Etano", methane1, ethane, False),
        ("Metano vs Butano", methane1, butane, False),
        ("Etano vs Butano", ethane, butane, False),
        ("Benzeno vs Benzeno", benzene1, benzene2, True),
    ]

    print("Teste                   | Esperado | Obtido   | Tempo (s)  | Status")
    print("--------------------------------------------------------------")

    all_passed = True
    for desc, g1, g2, expected in tests:
        start_time = time.perf_counter()
        result = are_isomorphic(g1, g2)
        end_time = time.perf_counter()
        elapsed = end_time - start_time

        status = "PASS" if result == expected else "FAIL"
        color = "\033[92m" if result == expected else "\033[91m"
        reset = "\033[0m"

        print(f"{desc:22} | {str(expected):8} | {str(result):8} | {elapsed:8.6f} | {color}{status}{reset}")

        if result != expected:
            all_passed = False

    return all_passed


def test_isomorphic_graphs():
    """Testes com grafos isomórficos gerados por permutação"""
    print("\n=== TESTES COM GRAFOS ISOMÓRFICOS GERADOS ===\n")

    original_graphs = [
        ("Metano", build_methane()),
        ("Etano", build_ethane()),
        ("Butano", build_butane()),
        ("Benzeno", build_benzene()),
    ]

    print("Grafo           | Permutações | Corretas | Tempo (s)  | Status")
    print("-------------------------------------------------------------")

    all_passed = True
    for name, original in original_graphs:
        isomorphic_graphs = _generate_isomorphic_group(original, max_permutations=5)

        start_time = time.perf_counter()
        correct = 0
        for graph in isomorphic_graphs:
            if are_isomorphic(original, graph):
                correct += 1
        end_time = time.perf_counter()
        elapsed = end_time - start_time

        total = len(isomorphic_graphs)
        status = "PASS" if correct == total else "FAIL"
        color = "\033[92m" if correct == total else "\033[91m"
        reset = "\033[0m"

        print(f"{name:15} | {total:11} | {correct:8} | {elapsed:8.6f} | {color}{status}{reset}")

        if correct != total:
            all_passed = False

    return all_passed


def test_non_isomorphic_graphs():
    """Testes com grafos não isomórficos"""
    print("\n=== TESTES COM GRAFOS NÃO ISOMÓRFICOS ===\n")

    g1 = build_cycle_graph(4)
    g2 = build_star_graph(5)

    regular_tests = build_small_non_isomorphic_regular_graphs()

    g5 = build_path_graph(5)
    g6 = build_cycle_graph(5)

    tests = [
                ("C4 vs K1,4", g1, g2, False),
                ("P5 vs C5", g5, g6, False),
            ] + regular_tests

    print("Teste                   | Esperado | Obtido   | Status")
    print("-------------------------------------------------------")

    all_passed = True
    for desc, g1, g2, expected in tests:
        result = are_isomorphic(g1, g2)
        status = "PASS" if result == expected else "FAIL"
        color = "\033[92m" if result == expected else "\033[91m"
        reset = "\033[0m"

        print(f"{desc:22} | {str(expected):8} | {str(result):8} | {color}{status}{reset}")

        if result != expected:
            all_passed = False
            print(f"  AVISO: Algoritmo falhou em distinguir grafos não isomórficos!")
            print(f"  Grafo 1: {len(g1.vertices())} vértices, {len(g1.arestas())} arestas")
            print(f"  Grafo 2: {len(g2.vertices())} vértices, {len(g2.arestas())} arestas")

    return all_passed


def test_regular_graphs():
    """Testes específicos para grafos regulares"""
    print("\n=== TESTES ESPECÍFICOS PARA GRAFOS REGULARES ===\n")

    regular_graphs = [
        ("C3 (triângulo)", build_cycle_graph(3)),
        ("C4", build_cycle_graph(4)),
        ("C5", build_cycle_graph(5)),
        ("C6", build_cycle_graph(6)),
        ("K4", build_complete_graph(4)),
        ("K5", build_complete_graph(5)),
        ("K3,3", build_complete_bipartite_graph(3, 3)),
        ("Petersen", build_petersen_graph()),
        ("W6 (roda)", build_wheel_graph(6)),
    ]

    print("Grafo              | Vértices | Grau | Auto-isomorfismo | Status")
    print("---------------------------------------------------------------")

    all_passed = True
    for name, graph in regular_graphs:
        result = are_isomorphic(graph, graph)
        n = len(graph.vertices())
        vertices_list = list(graph.vertices())
        k = graph.grau(vertices_list[0]) if n > 0 else 0

        status = "PASS" if result else "FAIL"
        color = "\033[92m" if result else "\033[91m"
        reset = "\033[0m"

        print(f"{name:18} | {n:8} | {k:4} | {str(result):15} | {color}{status}{reset}")

        if not result:
            all_passed = False

    return all_passed


def test_special_cases():
    """Testes com casos especiais e extremos"""
    print("\n=== TESTES COM CASOS ESPECIAIS ===\n")

    empty1 = Grafo()
    empty2 = Grafo()

    single1 = Grafo()
    single1.adicionar_vertice("A")
    single2 = Grafo()
    single2.adicionar_vertice("B")

    disconnected1 = Grafo()
    disconnected1.adicionar_vertice("X")
    disconnected1.adicionar_vertice("Y")

    disconnected2 = Grafo()
    disconnected2.adicionar_vertice("P")
    disconnected2.adicionar_vertice("Q")

    connected1 = Grafo()
    connected1.adicionar_aresta("A", "B")

    connected2 = Grafo()
    connected2.adicionar_aresta("X", "Y")

    tests = [
        ("Grafo vazio vs Grafo vazio", empty1, empty2, True),
        ("Vértice único vs Vértice único", single1, single2, True),
        ("2 vértices desconectados vs 2 vértices desconectados", disconnected1, disconnected2, True),
        ("2 vértices conectados vs 2 vértices conectados", connected1, connected2, True),
        ("Vértice único vs Grafo vazio", single1, empty1, False),
        ("2 vértices desconectados vs 2 vértices conectados", disconnected1, connected1, False),
    ]

    print("Teste                                       | Esperado | Obtido   | Status")
    print("------------------------------------------------------------------------")

    all_passed = True
    for desc, g1, g2, expected in tests:
        result = are_isomorphic(g1, g2)
        status = "PASS" if result == expected else "FAIL"
        color = "\033[92m" if result == expected else "\033[91m"
        reset = "\033[0m"

        print(f"{desc:43} | {str(expected):8} | {str(result):8} | {color}{status}{reset}")

        if result != expected:
            all_passed = False

    return all_passed


def get_memory_usage():
    """Retorna o uso de memória atual em MB com maior precisão"""
    gc.collect()
    process = psutil.Process(os.getpid())
    try:
        mem = process.memory_full_info().uss
    except AttributeError:
        mem = process.memory_info().rss
    return mem / 1024 / 1024


def linear_regression(x, y):
    """Calcula regressão linear sem numpy."""
    n = len(x)
    if n < 2:
        return 0, 0

    sum_x = sum(x)
    sum_y = sum(y)
    sum_xy = sum(xi * yi for xi, yi in zip(x, y))
    sum_x2 = sum(xi ** 2 for xi in x)

    denominator = n * sum_x2 - sum_x ** 2
    if denominator == 0:
        return 0, 0

    slope = (n * sum_xy - sum_x * sum_y) / denominator
    intercept = (sum_y - slope * sum_x) / n

    return slope, intercept


def analyze_complexity(times, sizes):
    """Analisa a complexidade com base nos tempos de execução usando regressão linear"""
    if len(times) < 2:
        return "", 1.0, 0.0

    min_time = 0.0001
    filtered_times = [max(t, min_time) for t in times]

    log_sizes = [math.log(s) for s in sizes]
    log_times = [math.log(t) for t in filtered_times]

    slope, intercept = linear_regression(log_sizes, log_times)

    if slope <= 0.3:
        complexity = "O(1)"
    elif slope <= 1.2:
        complexity = "O(n)"
    elif slope <= 1.8:
        complexity = "O(n log n)"
    elif slope <= 2.2:
        complexity = "O(n²)"
    elif slope <= 3.2:
        complexity = "O(n³)"
    else:
        complexity = f"O(n^{slope:.2f})"

    growth_factors = []
    for i in range(1, len(filtered_times)):
        if filtered_times[i - 1] > min_time:
            growth_factor = filtered_times[i] / filtered_times[i - 1]
            growth_factors.append(growth_factor)

    avg_growth = sum(growth_factors) / len(growth_factors) if growth_factors else 1.0

    return complexity, avg_growth, slope

def test_asymptotic_complexity():
    """Testes de complexidade assintótica com análise detalhada - ATÉ 100 VÉRTICES"""
    print("\n" + "=" * 80)
    print("TESTES DE COMPLEXIDADE ASSINTÓTICA - ATÉ 100 VÉRTICES")
    print("=" + "=" * 79)

    test_cases = [
        ("Path Graphs", lambda n: build_path_graph(n)),
        ("Cycle Graphs", lambda n: build_cycle_graph(n)),
        ("Complete Graphs", lambda n: build_complete_graph(n)),
        ("Star Graphs", lambda n: build_star_graph(n)),
        ("Bipartite Complete", lambda n: build_complete_bipartite_graph(n // 2, n // 2)),
        ("Wheel Graphs", lambda n: build_wheel_graph(n)),
    ]

    sizes = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

    print("\n📊 ANÁLISE DE COMPLEXIDADE ASSINTÓTICA (ATÉ 100 VÉRTICES):")
    print("=" * 100)

    all_results = {}

    for graph_type, builder in test_cases:
        print(f"\n🔍 ANALISANDO {graph_type.upper()} (ATÉ 100 VÉRTICES):")
        print("Vértices | Arestas  | Tempo (s)  | Fator Tempo | Complexidade")
        print("-" * 75)

        times = []
        vertices_counts = []

        for n in sizes:
            if "Bipartite" in graph_type:
                actual_n = (n // 2) * 2
                if actual_n < 4:
                    actual_n = 4
            else:
                actual_n = n

            if "Complete" in graph_type and actual_n > 60:
                continue

            g1 = builder(actual_n)
            g2 = builder(actual_n)

            time_measurements = []

            for _ in range(3):
                start_time = time.perf_counter()
                result = are_isomorphic(g1, g2)
                end_time = time.perf_counter()
                time_measurements.append(end_time - start_time)

            time_measurements.sort()
            elapsed = time_measurements[len(time_measurements) // 2]

            times.append(elapsed)
            vertices_counts.append(actual_n)

            n_edges = len(g1.arestas())

            if len(times) > 1:
                time_factor = times[-1] / times[-2]
            else:
                time_factor = 1.0

            if len(times) >= 3:
                current_complexity, _, _ = analyze_complexity(times, vertices_counts)
            else:
                current_complexity = ""

            status = "✓" if result else "✗"
            color = "\033[92m" if result else "\033[91m"
            reset = "\033[0m"

            print(f"{actual_n:8} | {n_edges:8} | {elapsed:10.6f} | "
                  f"{time_factor:10.2f} | {current_complexity:8} {color}{status}{reset}")

        complexity, avg_growth, exponent = analyze_complexity(times, vertices_counts)
        all_results[graph_type] = {
            'complexity': complexity,
            'avg_growth': avg_growth,
            'exponent': exponent,
            'max_time': max(times),
            'times': times,
            'sizes': vertices_counts
        }

        print(f"📈 Complexidade Final: {complexity}")
        print(f"📊 Expoente: {exponent:.2f}, Fator médio: {avg_growth:.2f}x")

    print("\n" + "=" * 80)
    print("📋 RESUMO DE COMPLEXIDADE POR TIPO DE GRAFO (ATÉ 100 VÉRTICES)")
    print("=" + "=" * 79)
    print("Tipo de Grafo        | Complexidade     | Expoente | Fator Cresc. | Tempo Máx (s)")
    print("-" * 90)

    for graph_type, result in all_results.items():
        complexity = result['complexity']
        exponent = result['exponent']
        avg_growth = result['avg_growth']
        max_time = result['max_time']

        if "O(1)" in complexity:
            color = "\033[92m"
        elif "O(n)" in complexity:
            color = "\033[96m"
        elif "O(n log n)" in complexity:
            color = "\033[94m"
        elif "O(n²)" in complexity:
            color = "\033[93m"
        else:
            color = "\033[91m"

        reset = "\033[0m"

        print(
            f"{graph_type:20} | {color}{complexity:16}{reset} | {exponent:8.2f} | {avg_growth:12.2f} | {max_time:12.6f}")

    return all_results

def test_memory_complexity():
    """Teste específico para análise de complexidade de memória - VERSÃO RÁPIDA"""
    print("\n" + "=" * 80)
    print("TESTE DE COMPLEXIDADE DE MEMÓRIA - VERSÃO RÁPIDA")
    print("=" + "=" * 79)

    sizes = [10, 15, 20, 25, 30]  

    print("\n📊 ANALISANDO CONSUMO DE MEMÓRIA (VERSÃO RÁPIDA):")
    print("Vértices | Tempo (s)  | Memória (MB) | Variação  | Fator Mem | Status")
    print("-" * 80)

    memory_diffs = []
    times = []
    previous_memory = 0.001

    for n in sizes:
        gc.collect()
        initial_memory = get_memory_usage()

        start_time = time.perf_counter()

        g1 = build_complete_graph(n)
        g2 = build_complete_graph(n)
        result = are_isomorphic(g1, g2)

        end_time = time.perf_counter()

        gc.collect()
        final_memory = get_memory_usage()

        elapsed = end_time - start_time
        memory_diff = max(0.0, final_memory - initial_memory)

        memory_diffs.append(memory_diff)
        times.append(elapsed)

        growth_factor = memory_diff / previous_memory if previous_memory > 0.001 else 1.0
        previous_memory = memory_diff

        status = "✓" if result else "✗"
        color = "\033[92m" if result else "\033[91m"
        reset = "\033[0m"

        print(
            f"{n:8} | {elapsed:10.6f} | {memory_diff:12.4f} | {memory_diff:9.4f} | {growth_factor:9.2f} | {color}{status}{reset}")

    print(f"\n💾 ANÁLISE DE COMPLEXIDADE DE MEMÓRIA:")

    if len(memory_diffs) >= 3:
        mem_complexity, mem_exponent = analyze_memory_complexity(memory_diffs, sizes)
        print(f"  • Complexidade de Memória: {mem_complexity}")
        print(f"  • Expoente: {mem_exponent:.2f}")
    else:
        mem_complexity = "Dados insuficientes"
        mem_exponent = 0

    max_mem = max(memory_diffs) if memory_diffs else 0
    if max_mem < 5.0:
        trend = "USO DE MEMÓRIA CONSTANTE/BAIXO"
        efficiency = "EXCELENTE"
        color = "\033[92m"
    elif max_mem < 20.0:
        trend = "USO DE MEMÓRIA LINEAR/MODERADO"
        efficiency = "MUITO BOA"
        color = "\033[96m"
    elif max_mem < 50.0:
        trend = "USO DE MEMÓRIA QUADRÁTICO"
        efficiency = "BOA"
        color = "\033[93m"
    else:
        trend = "USO DE MEMÓRIA ALTO"
        efficiency = "MODERADA"
        color = "\033[91m"

    reset = "\033[0m"
    print(f"  • Tendência: {color}{trend}{reset}")
    print(f"  • Eficiência: {color}{efficiency}{reset}")
    print(f"  • Consumo máximo: {max_mem:.4f}MB")

    return memory_diffs, times

def analyze_memory_complexity(memory_diffs, sizes):
    """Análise específica para complexidade de memória com critérios realistas"""
    if len(memory_diffs) < 3:
        return "Dados insuficientes", 1.0, 0.0

    filtered_memory = [max(m, 0.001) for m in memory_diffs]

    log_sizes = [math.log(s) for s in sizes]
    log_memory = [math.log(m) for m in filtered_memory]

    slope, intercept = linear_regression(log_sizes, log_memory)

    if slope <= 0.3:
        complexity = "O(1) - Constante"
    elif slope <= 1.2:
        complexity = "O(n) - Linear"
    elif slope <= 1.8:
        complexity = "O(n log n) - Log-linear"
    elif slope <= 2.2:
        complexity = "O(n²) - Quadrática"
    else:
        complexity = f"O(n^{slope:.2f}) - Polinomial"

    return complexity, slope


def run_comprehensive_cpu_analysis_isomorphism():
    """Executa análise completa de CPU para Isomorfismo - VERSÃO CORRIGIDA"""
    print("\n" + "=" * 100)
    print("ANÁLISE COMPREENSIVA DE CICLO DE PROCESSADOR - ISOMORFISMO")
    print("=" * 100)

    start_total = time.perf_counter()
    total_cpu_user = 0
    total_cpu_system = 0
    cpu_percentages = []

    def test_cpu_intensive():
        g1 = build_complete_graph(15)
        g2 = build_complete_graph(15)
        return are_isomorphic(g1, g2)

    result, cpu_metrics = measure_cpu_usage(test_cpu_intensive)
    print_cpu_metrics(cpu_metrics, "Grafos Completos K15")

    total_cpu_user += cpu_metrics['cpu_time_used']
    cpu_percentages.append(cpu_metrics['cpu_percent'])

    print("\n📈 ESCALABILIDADE DE CPU:")
    print("Tipo de Grafo       | Vértices | Tempo Real (s) | Tempo CPU (s) | Uso CPU (%) | Eficiência (%)")
    print("-" * 100)

    test_cases = [
        ("Completo", build_complete_graph),
        ("Ciclo", build_cycle_graph),
        ("Estrela", build_star_graph),
        ("Caminho", build_path_graph)
    ]

    for nome, builder in test_cases:
        for n in [10, 15, 20]:
            def test_escalabilidade():
                g1 = builder(n)
                g2 = builder(n)
                return are_isomorphic(g1, g2)

            result, metrics = measure_cpu_usage(test_escalabilidade)
            print(f"{nome:<18} | {n:8} | {metrics['execution_time']:13.6f} | {metrics['cpu_time_used']:12.6f} | "
                  f"{metrics['cpu_percent']:11.2f} | {metrics['efficiency']:13.1f}")

            total_cpu_user += metrics['cpu_time_used']
            cpu_percentages.append(metrics['cpu_percent'])

    print("\n🔬 TESTE COM GRAFOS COMPLEXOS:")
    print("Tipo de Grafo       | Vértices | Tempo Real (s) | Tempo CPU (s) | Uso CPU (%) | Eficiência (%)")
    print("-" * 100)

    complex_tests = [
        ("Petersen", criar_grafo_petersen()),
        ("Cúbico", criar_grafo_cubico()),
        ("Alpha-hélice 20", criar_proteina_alfa_helice_complexa(20)),
        ("Cafeína", criar_molecula_complexa())
    ]

    for nome, grafo in complex_tests:
        def test_complexo():
            return are_isomorphic(grafo, grafo)

        result, metrics = measure_cpu_usage(test_complexo)
        print(
            f"{nome:<18} | {len(grafo.vertices()):8} | {metrics['execution_time']:13.6f} | {metrics['cpu_time_used']:12.6f} | "
            f"{metrics['cpu_percent']:11.2f} | {metrics['efficiency']:13.1f}")

        total_cpu_user += metrics['cpu_time_used']
        cpu_percentages.append(metrics['cpu_percent'])

    end_total = time.perf_counter()
    total_time = end_total - start_total

    avg_cpu_percent = sum(cpu_percentages) / len(cpu_percentages) if cpu_percentages else 0

    print(f"\n📊 ESTATÍSTICAS GERAIS DE CPU:")
    print(f"  • Tempo total de execução: {total_time:.2f}s")
    print(f"  • Tempo total de CPU (usuário): {total_cpu_user:.4f}s")
    print(f"  • Uso médio de CPU: {avg_cpu_percent:.1f}%")

    cpu_efficiency = (total_cpu_user / total_time) * 100 if total_time > 0 else 0

    print(f"\n🎯 ANÁLISE DE EFICIÊNCIA:")
    if cpu_efficiency > 70:
        print("  • Alta eficiência de CPU - algoritmo bem otimizado")
    elif cpu_efficiency > 40:
        print("  • Eficiência moderada de CPU - desempenho aceitável")
    else:
        print("  • Baixa eficiência de CPU - possível gargalo de I/O ou espera")

    return {
        'total_time': total_time,
        'total_cpu_user': total_cpu_user,
        'total_cpu_system': 0,  
        'cpu_metrics': {
            'avg_cpu_percent': avg_cpu_percent,
            'efficiency': cpu_efficiency
        }
    }

def run_comprehensive_complexity_analysis():
    """Executa análise completa de complexidade"""
    print("\n" + "=" * 100)
    print("ANÁLISE COMPLETA DE COMPLEXIDADE ASSINTÓTICA DO ALGORITMO")
    print("=" + "=" * 99)

    print("""
    OBJETIVOS DA ANÁLISE:
    • Determinar complexidade temporal assintótica usando regressão linear
    • Analisar consumo de memória em função do tamanho da entrada
    • Identificar gargalos de performance
    • Fornecer recomendações de otimização baseadas em dados
    """)

    start_time = time.perf_counter()

    print("1. EXECUTANDO ANÁLISE DE COMPLEXIDADE TEMPORAL...")
    time_complexity_results = test_asymptotic_complexity()

    print("\n2. EXECUTANDO ANÁLISE DE COMPLEXIDADE DE MEMÓRIA...")
    memory_diffs, memory_times = test_memory_complexity()

    end_time = time.perf_counter()
    total_analysis_time = end_time - start_time

    print("\n" + "=" * 100)
    print("RELATÓRIO FINAL - COMPLEXIDADE ASSINTÓTICA")
    print("=" + "=" * 99)

    all_complexities = [result['complexity'] for result in time_complexity_results.values()]
    all_exponents = [result['exponent'] for result in time_complexity_results.values()]
    all_growths = [result['avg_growth'] for result in time_complexity_results.values()]

    complexity_count = {
        "O(1)": sum(1 for c in all_complexities if "O(1)" in c),
        "O(n)": sum(1 for c in all_complexities if "O(n)" in c and "log" not in c),
        "O(n log n)": sum(1 for c in all_complexities if "log" in c),
        "O(n²)": sum(1 for c in all_complexities if "O(n²)" in c or "O(n^2)" in c),
        "Outras": sum(
            1 for c in all_complexities if "O(1)" not in c and "O(n)" not in c and "log" not in c and "O(n²)" not in c)
    }

    predominant_complexity = max(complexity_count.items(), key=lambda x: x[1])
    avg_exponent = sum(all_exponents) / len(all_exponents)
    overall_avg_growth = sum(all_growths) / len(all_growths)

    print(f"📈 ESTATÍSTICAS GERAIS:")
    print(f"  • Tipos de grafos analisados: {len(time_complexity_results)}")
    print(
        f"  • Complexidade predominante: {predominant_complexity[0]} ({predominant_complexity[1]}/{len(all_complexities)})")
    print(f"  • Expoente médio: {avg_exponent:.2f}")
    print(f"  • Fator médio de crescimento: {overall_avg_growth:.2f}x")
    print(f"  • Tempo total de análise: {total_analysis_time:.2f}s")

    print(f"\n🎯 CLASSIFICAÇÃO DE PERFORMANCE:")
    if avg_exponent <= 1.0 and overall_avg_growth <= 1.5:
        rating = "EXCELENTE 🏆"
        color = "\033[92m"
        explanation = "Algoritmo altamente eficiente para todas as aplicações"
    elif avg_exponent <= 1.5 and overall_avg_growth <= 2.0:
        rating = "MUITO BOA 🥈"
        color = "\033[96m"
        explanation = "Adequado para aplicações em tempo real e grande escala"
    elif avg_exponent <= 2.0 and overall_avg_growth <= 3.0:
        rating = "BOA 🥉"
        color = "\033[93m"
        explanation = "Recomendado para a maioria dos casos de uso práticos"
    else:
        rating = "PRECISA DE OTIMIZAÇÃO ⚠️"
        color = "\033[91m"
        explanation = "Considerar otimizações para grafos muito grandes"

    reset = "\033[0m"
    print(f"  {color}{rating}{reset}")
    print(f"  {explanation}")

    print(f"\n🔮 PREVISÕES DE PERFORMANCE:")
    base_time = 0.001

    if avg_exponent <= 1.0:
        for vertices in [100, 500, 1000]:
            estimated_time = base_time * (vertices / 10) ** avg_exponent
            print(f"  • {vertices:4} vértices: ~{estimated_time:.4f}s")
        print(f"  • Performance mantida mesmo para grafos muito grandes")
    elif avg_exponent <= 1.5:
        for vertices in [100, 500]:
            estimated_time = base_time * (vertices / 10) ** avg_exponent
            print(f"  • {vertices:4} vértices: ~{estimated_time:.4f}s")
        print(f"  • Para 1000+ vértices: ~{base_time * (1000 / 10) ** avg_exponent:.4f}s (ainda viável)")
    else:
        for vertices in [100, 500]:
            estimated_time = base_time * (vertices / 10) ** avg_exponent
            print(f"  • {vertices:4} vértices: ~{estimated_time:.4f}s")
        print(f"  ⚠️  1000+ vértices: considerar otimizações ou heurísticas")

    print(f"\n💡 RECOMENDAÇÕES ESPECÍFICAS:")

    max_mem_used = max(memory_diffs) if memory_diffs else 0
    if max_mem_used < 5.0:
        print("  ✅ Excelente eficiência de memória - adequado para dispositivos com restrições")
    elif max_mem_used < 20.0:
        print("  ✅ Boa eficiência de memória - consumo moderado mesmo para grafos grandes")
    else:
        print("  ⚠️  Consumo de memória elevado - considerar otimizações para grafos muito grandes")

    if "O(n log n)" in all_complexities or "O(n²)" in all_complexities:
        print("  • Otimizações podem beneficiar grafos completos e densos")

    if overall_avg_growth > 2.0:
        print("  • Considerar cache de resultados intermediários para casos repetitivos")

    print("  • O algoritmo é robusto para a maioria dos casos de uso práticos")

    return time_complexity_results

def performance_test():
    """Testes de performance com grafos maiores - ATÉ 500 VÉRTICES"""
    print("\n=== TESTES DE PERFORMANCE (ATÉ 500 VÉRTICES) ===\n")

    performance_tests = [
        ("K_10", build_complete_graph(10)),
        ("K_20", build_complete_graph(20)),
        ("K_30", build_complete_graph(30)),
        ("C_50", build_cycle_graph(50)),
        ("C_100", build_cycle_graph(100)),
        ("C_200", build_cycle_graph(200)),
        ("C_500", build_cycle_graph(500)), 
        ("P_100", build_path_graph(100)),
        ("P_500", build_path_graph(500)),   
        ("S_100", build_star_graph(100)),
        ("S_500", build_star_graph(500)),   
        ("K10,10", build_complete_bipartite_graph(10, 10)),
        ("K20,20", build_complete_bipartite_graph(20, 20)),
        ("W_50", build_wheel_graph(50)),
        ("W_100", build_wheel_graph(100)),
    ]

    print("Grafo                      | Vértices | Arestas | Tempo (s)  | Status")
    print("------------------------------------------------------------------------")

    max_time = 0
    for name, graph in performance_tests:
        start_time = time.perf_counter()
        result = are_isomorphic(graph, graph)
        end_time = time.perf_counter()
        elapsed = end_time - start_time

        max_time = max(max_time, elapsed)
        n_vertices = len(graph.vertices())
        n_edges = len(graph.arestas())

        status = "PASS" if result else "FAIL"
        color = "\033[92m" if result else "\033[91m"
        reset = "\033[0m"

        print(f"{name:25} | {n_vertices:8} | {n_edges:7} | {elapsed:8.6f} | {color}{status}{reset}")

    print(f"\nTempo máximo: {max_time:.6f}s")

    if max_time < 0.01:
        print("Performance: Excelente ✓")
    elif max_time < 0.1:
        print("Performance: Muito Boa ✓")
    elif max_time < 1.0:
        print("Performance: Boa ✓")
    else:
        print(f"Performance: Aceitável ({max_time:.3f}s)")

def scalability_test():
    """Teste de escalabilidade com grafos progressivamente maiores - ATÉ 200 VÉRTICES"""
    print("\n=== TESTE DE ESCALABILIDADE (ATÉ 200 VÉRTICES) ===")

    print("Crescimento do tempo em função do tamanho do grafo (Ciclos C_n):")
    print("Vértices | Arestas | Tempo (s)  | Fator Cresc.")
    print("-" * 50)

    sizes = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200]
    times = []
    previous_time = 0

    for n in sizes:
        graph = build_cycle_graph(n)

        start_time = time.perf_counter()
        result = are_isomorphic(graph, graph)
        end_time = time.perf_counter()
        elapsed = end_time - start_time

        times.append(elapsed)
        n_edges = len(graph.arestas())

        growth_factor = elapsed / previous_time if previous_time > 0 else 0
        previous_time = elapsed

        status = "✓" if result else "✗"
        color = "\033[92m" if result else "\033[91m"
        growth_color = "\033[92m" if growth_factor <= 2.0 else "\033[93m" if growth_factor <= 3.0 else "\033[91m"
        reset = "\033[0m"

        if n == sizes[0]:
            print(f"{color}{n:8} | {n_edges:7} | {elapsed:9.6f} |     -{reset}")
        else:
            print(f"{color}{n:8} | {n_edges:7} | {elapsed:9.6f} | {growth_color}{growth_factor:7.2f}x{reset}")

    if len(times) > 1:
        growth_factors = []
        for i in range(1, len(times)):
            growth = times[i] / times[i - 1] if times[i - 1] > 0 else float('inf')
            growth_factors.append(growth)

        avg_growth = sum(growth_factors) / len(growth_factors)

        print(f"\nFator de crescimento médio: {avg_growth:.2f}x")

        if avg_growth < 1.2:
            complexity = "sub-linear"
            rating = "Excelente ✓"
            color = "\033[92m"
        elif avg_growth < 1.5:
            complexity = "linear"
            rating = "Muito boa ✓"
            color = "\033[92m"
        elif avg_growth < 2.0:
            complexity = "sub-quadrática"
            rating = "Boa ✓"
            color = "\033[93m"
        else:
            complexity = "quadrática ou superior"
            rating = "Aceitável ~"
            color = "\033[93m"

        reset = "\033[0m"
        print(f"Complexidade observada: {complexity}")
        print(f"{color}Avaliação: {rating}{reset}")

    print(f"\n--- ESCALABILIDADE COM GRAFOS CAMINHO (ATÉ 500 VÉRTICES) ---")
    print("Vértices | Tempo (s)  | Status")
    print("-" * 40)

    path_sizes = [100, 200, 300, 400, 500]
    for n in path_sizes:
        graph = build_path_graph(n)
        start_time = time.perf_counter()
        result = are_isomorphic(graph, graph)
        end_time = time.perf_counter()
        elapsed = end_time - start_time

        status = "✓" if result else "✗"
        color = "\033[92m" if result else "\033[91m"
        reset = "\033[0m"

        print(f"{color}{n:8} | {elapsed:9.6f} | {status}{reset}")


def stress_test():
    """Teste de estresse com múltiplas execuções - INCLUINDO 500 VÉRTICES"""
    print("\n=== TESTE DE ESTRESSE MULTINÍVEL (ATÉ 500 VÉRTICES) ===")

    print("\n🔬 TESTE COM GRAFOS PEQUENOS (K15):")
    test_graph_small = build_complete_graph(15)
    iterations_small = 20

    start_time_small = time.perf_counter()
    correct_results_small = 0

    for i in range(iterations_small):
        result = are_isomorphic(test_graph_small, test_graph_small)
        if result:
            correct_results_small += 1

    end_time_small = time.perf_counter()
    total_time_small = end_time_small - start_time_small
    avg_time_small = total_time_small / iterations_small
    success_rate_small = (correct_results_small / iterations_small) * 100

    print(f"  Resultados: {correct_results_small}/{iterations_small} corretos ({success_rate_small:.1f}%)")
    print(f"  Tempo total: {total_time_small:.4f}s")
    print(f"  Tempo médio: {avg_time_small:.6f}s")

    print("\n🔬 TESTE COM GRAFOS MÉDIOS (K50):")
    test_graph_medium = build_complete_graph(50)
    iterations_medium = 5

    start_time_medium = time.perf_counter()
    correct_results_medium = 0

    for i in range(iterations_medium):
        result = are_isomorphic(test_graph_medium, test_graph_medium)
        if result:
            correct_results_medium += 1

    end_time_medium = time.perf_counter()
    total_time_medium = end_time_medium - start_time_medium
    avg_time_medium = total_time_medium / iterations_medium
    success_rate_medium = (correct_results_medium / iterations_medium) * 100

    print(f"  Resultados: {correct_results_medium}/{iterations_medium} corretos ({success_rate_medium:.1f}%)")
    print(f"  Tempo total: {total_time_medium:.4f}s")
    print(f"  Tempo médio: {avg_time_medium:.6f}s")

    print("\n🔬 TESTE COM GRAFOS GRANDES (K100):")
    test_graph_large = build_complete_graph(100)
    iterations_large = 3

    start_time_large = time.perf_counter()
    correct_results_large = 0

    for i in range(iterations_large):
        result = are_isomorphic(test_graph_large, test_graph_large)
        if result:
            correct_results_large += 1

    end_time_large = time.perf_counter()
    total_time_large = end_time_large - start_time_large
    avg_time_large = total_time_large / iterations_large
    success_rate_large = (correct_results_large / iterations_large) * 100

    print(f"  Resultados: {correct_results_large}/{iterations_large} corretos ({success_rate_large:.1f}%)")
    print(f"  Tempo total: {total_time_large:.4f}s")
    print(f"  Tempo médio: {avg_time_large:.6f}s")

    print("\n🔬 TESTE COM GRAFOS MUITO GRANDES (Ciclo C200):")
    test_graph_huge = build_cycle_graph(200)
    iterations_huge = 3

    start_time_huge = time.perf_counter()
    correct_results_huge = 0

    for i in range(iterations_huge):
        result = are_isomorphic(test_graph_huge, test_graph_huge)
        if result:
            correct_results_huge += 1

    end_time_huge = time.perf_counter()
    total_time_huge = end_time_huge - start_time_huge
    avg_time_huge = total_time_huge / iterations_huge
    success_rate_huge = (correct_results_huge / iterations_huge) * 100

    print(f"  Resultados: {correct_results_huge}/{iterations_huge} corretos ({success_rate_huge:.1f}%)")
    print(f"  Tempo total: {total_time_huge:.4f}s")
    print(f"  Tempo médio: {avg_time_huge:.6f}s")

    print("\n🔬 TESTE COM GRAFOS DE 500 VÉRTICES (Ciclo C500):")
    test_graph_500 = build_cycle_graph(500)
    iterations_500 = 2

    start_time_500 = time.perf_counter()
    correct_results_500 = 0

    for i in range(iterations_500):
        result = are_isomorphic(test_graph_500, test_graph_500)
        if result:
            correct_results_500 += 1

    end_time_500 = time.perf_counter()
    total_time_500 = end_time_500 - start_time_500
    avg_time_500 = total_time_500 / iterations_500
    success_rate_500 = (correct_results_500 / iterations_500) * 100

    print(f"  Resultados: {correct_results_500}/{iterations_500} corretos ({success_rate_500:.1f}%)")
    print(f"  Tempo total: {total_time_500:.4f}s")
    print(f"  Tempo médio: {avg_time_500:.6f}s")

    print(f"\n📊 RESUMO DO TESTE DE ESTRESSE MULTINÍVEL (ATÉ 500 VÉRTICES):")
    print(f"  • Grafos pequenos (K15): {success_rate_small:.1f}% de sucesso")
    print(f"  • Grafos médios (K50): {success_rate_medium:.1f}% de sucesso")
    print(f"  • Grafos grandes (K100): {success_rate_large:.1f}% de sucesso")
    print(f"  • Grafos muito grandes (C200): {success_rate_huge:.1f}% de sucesso")
    print(f"  • Grafos de 500 vértices (C500): {success_rate_500:.1f}% de sucesso")

    total_tests = (iterations_small + iterations_medium + iterations_large +
                   iterations_huge + iterations_500)
    total_correct = (correct_results_small + correct_results_medium + correct_results_large +
                    correct_results_huge + correct_results_500)
    overall_success_rate = (total_correct / total_tests) * 100
    total_time = (total_time_small + total_time_medium + total_time_large +
                  total_time_huge + total_time_500)

    print(f"\n🎯 ESTATÍSTICAS GERAIS:")
    print(f"  • Total de verificações: {total_tests}")
    print(f"  • Total corretas: {total_correct}")
    print(f"  • Taxa de sucesso geral: {overall_success_rate:.1f}%")
    print(f"  • Tempo total de execução: {total_time:.4f}s")

    if overall_success_rate == 100:
        print("\n\033[92m✓ TESTE DE ESTRESSE MULTINÍVEL: PASSOU!\033[0m")
        print("  O algoritmo manteve 100% de precisão mesmo com grafos de 500 vértices")
        return True
    else:
        print(f"\n\033[91m✗ TESTE DE ESTRESSE MULTINÍVEL: FALHOU!\033[0m")
        print(f"  {total_tests - total_correct} verificação(ões) incorreta(s)")
        return False

def stress_test_very_large():
    """Teste de estresse específico para grafos muito grandes - ATÉ 500 VÉRTICES"""
    print("\n=== TESTE DE ESTRESSE - GRAFOS GRANDES (ATÉ 500 VÉRTICES) ===")

    test_cases = [
        ("Ciclo C100", build_cycle_graph(100), 5),
        ("Ciclo C200", build_cycle_graph(200), 3),
        ("Ciclo C500", build_cycle_graph(500), 2),  
        ("Path P100", build_path_graph(100), 3),
        ("Path P500", build_path_graph(500), 2),    
        ("Estrela S100", build_star_graph(100), 5),
        ("Estrela S500", build_star_graph(500), 2),
        ("Grafo Bipartido K50,50", build_complete_bipartite_graph(50, 50), 2),
    ]

    print("Testando grafos grandes (até 500 vértices)...")
    print("Grafo                | Vértices | Arestas | Iterações | Tempo Médio (s) | Status")
    print("--------------------------------------------------------------------------------")

    all_passed = True
    total_time = 0
    total_tests = 0
    total_correct = 0

    for name, graph, iterations in test_cases:
        n_vertices = len(graph.vertices())
        n_edges = len(graph.arestas())

        start_time = time.perf_counter()
        correct = 0

        for i in range(iterations):
            result = are_isomorphic(graph, graph)
            if result:
                correct += 1

        end_time = time.perf_counter()
        elapsed = end_time - start_time
        avg_time = elapsed / iterations
        total_time += elapsed
        total_tests += iterations
        total_correct += correct

        status = "PASS" if correct == iterations else "FAIL"
        color = "\033[92m" if correct == iterations else "\033[91m"
        reset = "\033[0m"

        print(f"{name:20} | {n_vertices:8} | {n_edges:7} | {iterations:9} | {avg_time:15.6f} | {color}{status}{reset}")

        if correct != iterations:
            all_passed = False

    success_rate = (total_correct / total_tests) * 100

    print(f"\n📊 RESUMO - GRAFOS GRANDES (ATÉ 500 VÉRTICES):")
    print(f"  • Total de verificações: {total_tests}")
    print(f"  • Verificações corretas: {total_correct}")
    print(f"  • Taxa de sucesso: {success_rate:.1f}%")
    print(f"  • Tempo total: {total_time:.4f}s")
    print(f"  • Tempo médio por verificação: {total_time / total_tests:.6f}s")

    if all_passed:
        print("\n\033[92m✓ TESTE COM GRAFOS ATÉ 500 VÉRTICES: PASSOU!\033[0m")
    else:
        print("\n\033[91m✗ TESTE COM GRAFOS ATÉ 500 VÉRTICES: FALHOU!\033[0m")

    return all_passed

def test_complex_graphs():
    """Testes com grafos complexos de diferentes tipos"""
    print("\n=== TESTES COM GRAFOS COMPLEXOS ===\n")

    complex_tests = []

    print("1. TESTANDO GRAFOS MOLECULARES:")
    methane = build_methane()
    ethane = build_ethane()
    butane = build_butane()
    benzene = build_benzene()
    caffeine = criar_molecula_complexa()

    complex_tests.extend([
        ("Metano vs Metano", methane, methane, True),
        ("Etano vs Etano", ethane, ethane, True),
        ("Butano vs Butano", butane, butane, True),
        ("Benzeno vs Benzeno", benzene, benzene, True),
        ("Cafeína vs Cafeína", caffeine, caffeine, True),
        ("Metano vs Etano", methane, ethane, False),
        ("Etano vs Butano", ethane, butane, False),
        ("Benzeno vs Cafeína", benzene, caffeine, False),
    ])

    print("2. TESTANDO GRAFOS TEÓRICOS:")
    cycle_6 = build_cycle_graph(6)
    cycle_8 = build_cycle_graph(8)
    path_5 = build_path_graph(5)
    path_7 = build_path_graph(7)
    complete_5 = build_complete_graph(5)
    complete_6 = build_complete_graph(6)
    bipartite_33 = build_complete_bipartite_graph(3, 3)
    bipartite_24 = build_complete_bipartite_graph(2, 4)
    petersen = build_petersen_graph()
    petersen2 = criar_grafo_petersen()
    wheel_6 = build_wheel_graph(6)
    wheel_8 = build_wheel_graph(8)
    star_5 = build_star_graph(5)
    star_7 = build_star_graph(7)
    cubic = criar_grafo_cubico()

    complex_tests.extend([
        ("C6 vs C6", cycle_6, cycle_6, True),
        ("C6 vs C8", cycle_6, cycle_8, False),
        ("C6 vs P6", cycle_6, build_path_graph(6), False),
        ("P5 vs P5", path_5, path_5, True),
        ("P5 vs P7", path_5, path_7, False),
        ("K5 vs K5", complete_5, complete_5, True),
        ("K5 vs K6", complete_5, complete_6, False),
        ("K3,3 vs K3,3", bipartite_33, bipartite_33, True),
        ("K3,3 vs K2,4", bipartite_33, bipartite_24, False),
        ("Petersen vs Petersen", petersen, petersen, True),
        ("Petersen1 vs Petersen2", petersen, petersen2, True),
        ("Petersen vs C10", petersen, build_cycle_graph(10), False),
        ("W6 vs W6", wheel_6, wheel_6, True),
        ("W6 vs W8", wheel_6, wheel_8, False),
        ("S5 vs S5", star_5, star_5, True),
        ("S5 vs S7", star_5, star_7, False),
        ("Cúbico vs Cúbico", cubic, cubic, True),
        ("Cúbico vs Petersen", cubic, petersen, False),
    ])

    print("3. TESTANDO GRAFOS DE PROTEÍNAS:")
    alpha_helix = criar_proteina_alfa_helice_complexa(12)
    alpha_helix_same = criar_proteina_alfa_helice_complexa(12)
    alpha_helix_diff = criar_proteina_alfa_helice_complexa(16)
    alpha_modified = criar_proteina_alfa_helice_modificada(12)
    beta_sheet_anti = criar_proteina_beta_folha_complexa(3, 6)
    beta_sheet_anti_same = criar_proteina_beta_folha_complexa(3, 6)
    beta_sheet_anti_diff = criar_proteina_beta_folha_complexa(4, 6)
    beta_sheet_para = criar_proteina_beta_folha_paralela(3, 6)

    complex_tests.extend([
        ("α-hélice 12 vs α-hélice 12", alpha_helix, alpha_helix_same, True),
        ("α-hélice 12 vs α-hélice 16", alpha_helix, alpha_helix_diff, False),
        ("α-hélice vs α-hélice modificada", alpha_helix, alpha_modified, False),
        ("β-folha anti 3x6 vs β-folha anti 3x6", beta_sheet_anti, beta_sheet_anti_same, True),
        ("β-folha anti 3x6 vs β-folha anti 4x6", beta_sheet_anti, beta_sheet_anti_diff, False),
        ("β-folha anti vs β-folha para", beta_sheet_anti, beta_sheet_para, False),
        ("α-hélice vs β-folha anti", alpha_helix, beta_sheet_anti, False),
    ])

    print("4. EXECUTANDO TODOS OS TESTES COMPLEXOS:")
    print("Teste                                       | Esperado | Obtido   | Tempo (s)  | Status")
    print("------------------------------------------------------------------------")

    passed = 0
    failed = 0
    total_time = 0

    for desc, g1, g2, expected in complex_tests:
        start_time = time.perf_counter()
        result = are_isomorphic(g1, g2)
        end_time = time.perf_counter()
        elapsed = end_time - start_time
        total_time += elapsed

        status = result == expected
        if status:
            passed += 1
            color = "\033[92m"
            symbol = "PASS"
        else:
            failed += 1
            color = "\033[91m"
            symbol = "FAIL"

        reset = "\033[0m"
        print(f"{desc:43} | {str(expected):8} | {str(result):8} | {elapsed:8.6f} | {color}{symbol}{reset}")

    total_tests = passed + failed
    success_rate = (passed / total_tests) * 100

    print(f"\nRESUMO DOS TESTES COMPLEXOS:")
    print(f"Total de testes: {total_tests}")
    print(f"Passaram: {passed} ({success_rate:.1f}%)")
    print(f"Falharam: {failed}")
    print(f"Tempo total: {total_time:.4f}s")
    print(f"Tempo médio por teste: {total_time / total_tests:.6f}s")

    if failed == 0:
        print("\033[92m✓ TODOS OS TESTES COMPLEXOS PASSARAM!\033[0m")
        return True
    else:
        print(f"\033[91m✗ {failed} TESTE(S) FALHARAM!\033[0m")
        return False


def main():
    """Função principal com todos os testes"""
    print("=" * 80)
    print("TESTES COMPREENSIVOS PARA ALGORITMO DE ISOMORFISMO - GRAFOS COMPLEXOS")
    print("=" * 80)
    print("=== AVALIAÇÃO COMPLETA DO ALGORITMO COM ANÁLISE DE COMPLEXIDADE ===\n")

    all_passed = True

    print("EXECUTANDO TESTES BÁSICOS...")
    if not test_basic_isomorphism():
        all_passed = False

    print("\nEXECUTANDO TESTES COM GRAFOS ISOMÓRFICOS...")
    if not test_isomorphic_graphs():
        all_passed = False

    print("\nEXECUTANDO TESTES COM GRAFOS NÃO ISOMÓRFICOS...")
    if not test_non_isomorphic_graphs():
        all_passed = False

    print("\nEXECUTANDO TESTES COM GRAFOS REGULARES...")
    if not test_regular_graphs():
        all_passed = False

    print("\nEXECUTANDO TESTES COM CASOS ESPECIAIS...")
    if not test_special_cases():
        all_passed = False

    print("\n" + "=" * 80)
    print("INICIANDO ANÁLISE DE COMPLEXIDADE ASSINTÓTICA")
    print("=" * 80)
    complexity_results = run_comprehensive_complexity_analysis()

    print("\nEXECUTANDO TESTES DE PERFORMANCE...")
    performance_test()

    print("\nEXECUTANDO TESTE DE ESCALABILIDADE...")
    scalability_test()

    print("\nEXECUTANDO TESTE DE ESTRESSE...")
    if not stress_test():
        all_passed = False

    print("\nEXECUTANDO TESTE DE ESTRESSE COM GRAFOS MUITO GRANDES...")
    if not stress_test_very_large():
        all_passed = False

    print("\nEXECUTANDO TESTES COM GRAFOS COMPLEXOS...")
    if not test_complex_graphs():
        all_passed = False

    # =========================================================================
    # NOVA SEÇÃO: ANÁLISE DE CPU
    # =========================================================================
    print("\n" + "=" * 80)
    print("INICIANDO ANÁLISE DE CPU COM MÉTRICAS CORRIGIDAS")
    print("=" * 80)
    cpu_results = run_cpu_analysis()

    print("\n" + "=" * 80)
    print("RELATÓRIO FINAL - ALGORITMO DE ISOMORFISMO")
    print("=" * 80)

    if all_passed:
        print("\033[92m✓ TODOS OS TESTES PASSARAM!\033[0m")
        print("O algoritmo demonstrou correção e robustez em todos os cenários testados.")
    else:
        print("\033[91m✗ ALGUNS TESTES FALHARAM!\033[0m")
        print("O algoritmo precisa de ajustes para lidar com certos casos.")

    print(f"\n📊 RESUMO DA ANÁLISE DE COMPLEXIDADE:")
    complexities = [result['complexity'] for result in complexity_results.values()]
    avg_exponent = sum(result['exponent'] for result in complexity_results.values()) / len(complexity_results)
    avg_growth = sum(result['avg_growth'] for result in complexity_results.values()) / len(complexity_results)

    print(f"  • Complexidade predominante: {max(set(complexities), key=complexities.count)}")
    print(f"  • Expoente médio: {avg_exponent:.2f}")
    print(f"  • Fator de crescimento médio: {avg_growth:.2f}x")
    print(f"  • Tipos de grafos analisados: {len(complexity_results)}")

    # =========================================================================
    # SEÇÃO CORRIGIDA: RESUMO DE CPU
    # =========================================================================
    print(f"\n⚡ RESUMO DA ANÁLISE DE CPU:")
    print(f"  • Tempo total de CPU (usuário): {cpu_results['total_cpu_user']:.4f}s")
    print(f"  • Tempo total de CPU (sistema): {cpu_results['total_cpu_system']:.4f}s")

    total_cpu_time = cpu_results['total_cpu_user'] + cpu_results['total_cpu_system']
    if cpu_results['total_time'] > 0:
        overall_efficiency = (total_cpu_time / cpu_results['total_time']) * 100
    else:
        overall_efficiency = 0

    print(f"  • Eficiência geral: {overall_efficiency:.1f}%")
    print(f"  • Uso médio de CPU: {cpu_results['cpu_metrics']['avg_cpu_percent']:.1f}%")

    cpu_efficiency = overall_efficiency / 100.0

    if cpu_efficiency > 0.7:
        cpu_status = "EXCELENTE"
        cpu_color = "\033[92m"
    elif cpu_efficiency > 0.4:
        cpu_status = "BOA"
        cpu_color = "\033[93m"
    else:
        cpu_status = "MODERADA"
        cpu_color = "\033[91m"

    reset = "\033[0m"
    print(f"  • Eficiência de CPU: {cpu_color}{cpu_status}{reset}")

    if avg_exponent <= 1.0 and avg_growth <= 1.5:
        assessment = "EXCELENTE - Pronto para aplicações em grande escala"
        color = "\033[92m"
    elif avg_exponent <= 1.5 and avg_growth <= 2.0:
        assessment = "MUITO BOA - Adequado para a maioria dos casos práticos"
        color = "\033[96m"
    else:
        assessment = "BOA - Recomendado para grafos de tamanho moderado"
        color = "\033[93m"

    reset = "\033[0m"
    print(f"  • Avaliação: {color}{assessment}{reset}")

    print(f"\n🔮 PREVISÕES DE PERFORMANCE:")
    base_time = 0.001

    if avg_exponent <= 1.0:
        for vertices in [100, 500, 1000]:
            estimated_time = base_time * (vertices / 10) ** avg_exponent
            print(f"  • {vertices:4} vértices: ~{estimated_time:.4f}s")
        print(f"  • Performance mantida mesmo para grafos muito grandes")
    elif avg_exponent <= 1.5:
        for vertices in [100, 500]:
            estimated_time = base_time * (vertices / 10) ** avg_exponent
            print(f"  • {vertices:4} vértices: ~{estimated_time:.4f}s")
        print(f"  • Para 1000+ vértices: ~{base_time * (1000 / 10) ** avg_exponent:.4f}s (ainda viável)")
    else:
        for vertices in [100, 500]:
            estimated_time = base_time * (vertices / 10) ** avg_exponent
            print(f"  • {vertices:4} vértices: ~{estimated_time:.4f}s")
        print(f"  ⚠️  1000+ vértices: considerar otimizações ou heurísticas")

    print(f"\n💡 RECOMENDAÇÕES ESPECÍFICAS BASEADAS NAS MÉTRICAS DE CPU:")

    if cpu_results['cpu_metrics']['avg_cpu_percent'] < 50:
        print("  • Baixo uso de CPU - algoritmo pode beneficiar de paralelização")
    elif cpu_results['cpu_metrics']['avg_cpu_percent'] > 90:
        print("  • Alto uso de CPU - algoritmo já está bem otimizado")

    if cpu_results['total_cpu_system'] > cpu_results['total_cpu_user'] * 2:
        print("  • Alto tempo de sistema - verificar operações de I/O")

    if cpu_efficiency < 0.4:
        print("  • Baixa eficiência de CPU - investigar possíveis gargalos")

    print(f"\n🔬 APLICAÇÕES RECOMENDADAS:")
    if avg_exponent <= 1.0:
        print("  ✅ Química computacional (moléculas complexas)")
        print("  ✅ Bioinformática (proteínas e estruturas grandes)")
        print("  ✅ Análise de redes sociais (grafos grandes)")
        print("  ✅ Aplicações em tempo real")
        print("  ✅ Processamento de grafos em streaming")
    elif avg_exponent <= 1.5:
        print("  ✅ Química computacional (moléculas médias e complexas)")
        print("  ✅ Bioinformática (proteínas e estruturas secundárias)")
        print("  ✅ Teoria dos grafos (até ~1000 vértices)")
        print("  ✅ Análise de redes (grafos moderados a grandes)")
    else:
        print("  ✅ Química computacional (moléculas pequenas e médias)")
        print("  ✅ Bioinformática (peptídeos e proteínas pequenas)")
        print("  ✅ Teoria dos grafos (até ~500 vértices)")
        print("  ⚠️  Para grafos maiores: considerar otimizações")

    print(f"\n📈 MÉTRICAS FINAIS DE DESEMPENHO:")
    print(f"  • Correção algorítmica: {'✓ PASSOU' if all_passed else '✗ FALHOU'}")
    print(f"  • Complexidade: {avg_exponent:.2f} (expoente médio)")
    print(f"  • Eficiência de CPU: {cpu_efficiency:.1%}")
    print(f"  • Uso de recursos: {cpu_status}")

    print("\n" + "=" * 80)

    return {
        'all_passed': all_passed,
        'complexity_results': complexity_results,
        'cpu_results': cpu_results,
        'avg_exponent': avg_exponent,
        'avg_growth': avg_growth,
        'cpu_efficiency': cpu_efficiency
    }


if __name__ == "__main__":
    main()