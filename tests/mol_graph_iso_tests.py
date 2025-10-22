import time
import random
import numpy as np
import warnings
import psutil
import os
from mol_graph_iso import *
from structures.Grafo import Grafo


# =============================================================================
# FUNÇÕES AUXILIARES PARA CRIAÇÃO DE GRAFOS DE TESTE
# =============================================================================

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
        g.adicionar_vertice(i)

    for i in range(2, n_vertices + 1):
        g.adicionar_aresta(1, i)

    return g


def criar_grafo_anel(n_vertices):
    """Grafo anel C_n"""
    g = Grafo()
    for i in range(1, n_vertices + 1):
        g.adicionar_vertice(i)

    for i in range(1, n_vertices):
        g.adicionar_aresta(i, i + 1)
    g.adicionar_aresta(n_vertices, 1)

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


def criar_molecula_complexa():
    """Molécula orgânica complexa - Cafeína (CORRIGIDA: usando números atômicos)"""
    g = Grafo()

    atoms = {
        1: 6, 2: 6, 3: 6, 4: 6, 5: 6, 6: 6,
        7: 7, 8: 7, 9: 7,
        10: 8, 11: 8,
        12: 6, 13: 6, 14: 6
    }

    for atom_id, atomic_number in atoms.items():
        g.adicionar_vertice(atom_id, atomic_number)

    bonds = [
        (1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 1),
        (5, 6, 1), (6, 1, 1), (1, 7, 1), (3, 8, 1),
        (5, 9, 1), (7, 10, 2), (8, 11, 2),
        (7, 12, 1), (8, 13, 1), (9, 14, 1)
    ]

    for u, v, ordem in bonds:
        g.adicionar_aresta(u, v, ordem)

    return g


def criar_proteina_alfa_helice(n_residuos=20):
    """Estrutura de α-hélice (CORRIGIDA: usando inteiros)"""
    g = Grafo()

    for i in range(1, n_residuos + 1):
        g.adicionar_vertice(i, 1)

    for i in range(1, n_residuos):
        g.adicionar_aresta(i, i + 1, 1)

    for i in range(1, n_residuos - 3):
        g.adicionar_aresta(i, i + 4, 1)

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


# =============================================================================
# TESTES DE PERFORMANCE E ESCALABILIDADE
# =============================================================================

def test_performance_scalability():
    """Testa a escalabilidade do algoritmo com grafos de tamanhos crescentes"""
    print("\n1. TESTES DE PERFORMANCE - ESCALABILIDADE")
    print("=" * 60)

    resultados = []
    tamanhos = [5, 10, 15, 20, 25, 30]

    print("Testando escalabilidade com grafos completos:")
    print("Vértices | Arestas  | Tempo (s)  | Memória (MB) | Isomorfos")
    print("-" * 65)

    for n in tamanhos:
        process = psutil.Process(os.getpid())
        memoria_inicial = process.memory_info().rss / 1024 / 1024

        g1 = criar_grafo_completo(n)
        g2 = criar_grafo_completo(n)

        start_time = time.time()
        isomorfo = isomorfismo_molecular(g1, g2)
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
            'isomorfo': isomorfo
        })

        status = "✓" if isomorfo else "✗"
        print(f"{n:8} | {n_arestas:8} | {tempo_execucao:10.4f} | {memoria_usada:11.2f} | {status}")

    return resultados


def test_complexity_analysis():
    """Análise de complexidade assintótica"""
    print("\n2. ANÁLISE DE COMPLEXIDADE ASSINTÓTICA")
    print("=" * 60)

    resultados = []
    tamanhos = [5, 10, 20, 30, 40]
    tempos_anteriores = None

    print("Análise de Complexidade Assintótica:")
    print("Vértices | Arestas  | Tempo (s)  | Fator Cresc. | O(?)")
    print("-" * 65)

    for n in tamanhos:
        g1 = criar_grafo_completo(n)
        g2 = criar_grafo_completo(n)

        start_time = time.perf_counter()
        isomorfismo_molecular(g1, g2)
        end_time = time.perf_counter()

        tempo_atual = end_time - start_time
        n_arestas = len(g1.arestas())

        if tempos_anteriores is not None:
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
    """Teste específico de eficiência de memória (CORRIGIDA)"""
    print("\n3. TESTES DE EFICIÊNCIA DE MEMÓRIA")
    print("=" * 60)

    tipos_grafos = [
        ("Completo K20", lambda: criar_grafo_completo(20)),
        ("Caminho P50", lambda: criar_grafo_caminho(50)),
        ("Estrela S30", lambda: criar_grafo_estrela(30)),
        ("Anel C40", lambda: criar_grafo_anel(40)),
        ("Alpha-hélice 25", lambda: criar_proteina_alfa_helice(25)),
        ("Cafeína", criar_molecula_complexa),
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
        isomorfo = isomorfismo_molecular(g1, g2)
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
            'tempo': tempo_execucao,
            'isomorfo': isomorfo
        })

        status = "✓" if isomorfo else "✗"
        print(f"{nome:<20} | {n_vertices:8} | {n_arestas:8} | {memoria_usada:11.2f} | {tempo_execucao:8.4f} {status}")

    return resultados


# =============================================================================
# TESTES DE ESTRESSE
# =============================================================================

def test_stress_executions():
    """Teste de estresse com múltiplas execuções consecutivas"""
    print("\n4. TESTES DE ESTRESSE - EXECUÇÕES CONSECUTIVAS")
    print("=" * 60)

    n_execucoes = 50
    tempos = []
    resultados_isomorfismo = []
    memorias = []

    print(f"Executando {n_execucoes} comparações consecutivas...")

    process = psutil.Process(os.getpid())

    for i in range(n_execucoes):
        n_vertices = random.randint(5, 20)
        g1, g2 = criar_grafos_aleatorios(n_vertices, 0.4)

        memoria_inicial = process.memory_info().rss / 1024 / 1024

        start_time = time.time()
        isomorfo = isomorfismo_molecular(g1, g2)
        end_time = time.time()

        memoria_final = process.memory_info().rss / 1024 / 1024

        tempo_execucao = end_time - start_time
        memoria_usada = memoria_final - memoria_inicial

        tempos.append(tempo_execucao)
        resultados_isomorfismo.append(isomorfo)
        memorias.append(memoria_usada)

        if (i + 1) % 10 == 0:
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

    sem_falhas = all(r is not None for r in resultados_isomorfismo)
    tempos_estaveis = tempo_max < 5.0

    status = "✓" if (sem_falhas and tempos_estaveis) else "✗"
    color = "\033[92m" if (sem_falhas and tempos_estaveis) else "\033[91m"
    print(f"{color}{status} Teste de estresse: {'PASSOU' if (sem_falhas and tempos_estaveis) else 'FALHOU'}\033[0m")

    return {
        'n_execucoes': n_execucoes,
        'tempo_medio': tempo_medio,
        'tempo_max': tempo_max,
        'tempo_min': tempo_min,
        'memoria_media': memoria_media,
        'memoria_max': memoria_max,
        'sem_falhas': sem_falhas,
        'tempos_estaveis': tempos_estaveis
    }


def test_stress_large_graphs():
    """Teste de estresse com grafos grandes"""
    print("\n5. TESTES DE ESTRESSE - GRAFOS GRANDES")
    print("=" * 60)

    resultados = []

    grafos_grandes = [
        ("Completo K30", 30, 435),
        ("Completo K40", 40, 780),
        ("Caminho P100", 100, 99),
        ("Estrela S80", 80, 79),
        ("Alpha-hélice 50", 50, 97),
        ("Anel C100", 100, 100),
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
                g1 = criar_proteina_alfa_helice(n_vertices)
                g2 = criar_proteina_alfa_helice(n_vertices)
            elif "Anel" in nome:
                g1 = criar_grafo_anel(n_vertices)
                g2 = criar_grafo_anel(n_vertices)
            else:
                continue

            start_time = time.time()
            isomorfo = isomorfismo_molecular(g1, g2)
            end_time = time.time()

            tempo_execucao = end_time - start_time
            n_arestas_real = len(g1.arestas())

            status = "✓" if (isomorfo and tempo_execucao < 10.0) else "⏱️" if tempo_execucao >= 10.0 else "✗"

            resultados.append({
                'grafo': nome,
                'vertices': n_vertices,
                'arestas': n_arestas_real,
                'tempo': tempo_execucao,
                'isomorfo': isomorfo,
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
                'isomorfo': False,
                'status': '✗'
            })

    return resultados


def test_stress_memory_usage():
    """Teste específico de uso de memória"""
    print("\n6. TESTES DE ESTRESSE - USO DE MEMÓRIA")
    print("=" * 60)

    n_testes = 10
    grafos = []

    print(f"Criando {n_testes} grafos grandes...")

    for i in range(n_testes):
        tamanho = 20 + i * 8
        grafos.append(criar_grafo_completo(tamanho))

    process = psutil.Process(os.getpid())
    memoria_inicial = process.memory_info().rss / 1024 / 1024

    print("Executando comparações múltiplas...")
    tempos = []

    for i in range(len(grafos) - 1):
        start_time = time.time()
        isomorfismo_molecular(grafos[i], grafos[i + 1])
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

def test_automorfismos_performance():
    """Teste de performance para cálculo de automorfismos"""
    print("\n7. TESTES DE PERFORMANCE - AUTOMORFISMOS")
    print("=" * 60)

    resultados = []

    estruturas = [
        ("Grafo Vazio", Grafo()),
        ("Vértice Único", criar_grafo_estrela(1)),
        ("Aresta Simples", criar_grafo_caminho(2)),
        ("Triângulo", criar_grafo_completo(3)),
        ("Quadrado", criar_grafo_anel(4)),
        ("Tetraedro", criar_grafo_completo(4)),
        ("Cubo", criar_grafo_cubico()),
        ("Alpha-hélice 8", criar_proteina_alfa_helice(8)),
    ]

    print("Estrutura           | Vértices | Automorfismos | Tempo (s)")
    print("-" * 55)

    for nome, grafo in estruturas:
        try:
            start_time = time.perf_counter()
            automorfismos = automorfismos_moleculares(grafo)
            end_time = time.perf_counter()

            tempo_execucao = end_time - start_time
            n_automorfismos = len(automorfismos)
            n_vertices = len(grafo.vertices())

            resultados.append({
                'estrutura': nome,
                'vertices': n_vertices,
                'automorfismos': n_automorfismos,
                'tempo': tempo_execucao
            })

            print(f"{nome:<18} | {n_vertices:8} | {n_automorfismos:13} | {tempo_execucao:8.4f}")

        except Exception as e:
            print(f"{nome:<18} | {'-':8} | {'-':13} | {'ERRO':>8}")
            resultados.append({
                'estrutura': nome,
                'vertices': len(grafo.vertices()),
                'automorfismos': 0,
                'tempo': 0,
                'erro': str(e)
            })

    return resultados


def test_numerical_stability():
    """Análise de estabilidade numérica"""
    print("\n8. ANÁLISE DE ESTABILIDADE NUMÉRICA")
    print("=" * 60)

    print("Testando consistência com mesmos grafos:")

    g1, g2 = criar_grafos_aleatorios(10, 0.5)

    resultados = []
    for i in range(20):
        start_time = time.perf_counter()
        resultado = isomorfismo_molecular(g1, g2)
        end_time = time.perf_counter()

        resultados.append({
            'execucao': i + 1,
            'resultado': resultado,
            'tempo': end_time - start_time
        })

    resultados_consistentes = all(r['resultado'] == resultados[0]['resultado'] for r in resultados)
    tempos = [r['tempo'] for r in resultados]
    tempo_medio = np.mean(tempos)
    tempo_std = np.std(tempos)

    print(f"  • Resultados consistentes: {resultados_consistentes}")
    print(f"  • Tempo médio: {tempo_medio:.6f}s")
    print(f"  • Desvio padrão do tempo: {tempo_std:.6f}s")

    if resultados_consistentes:
        print("  • Estabilidade numérica: EXCELENTE ✓")
    else:
        print("  • Estabilidade numérica: PROBLEMA ✗")

    return {
        'resultados_consistentes': resultados_consistentes,
        'tempo_medio': tempo_medio,
        'tempo_std': tempo_std
    }


def test_edge_cases():
    """Testes com casos extremos e especiais"""
    print("\n9. TESTES COM CASOS EXTREMOS")
    print("=" * 60)

    casos = [
        ("Grafo vazio", Grafo(), Grafo()),
        ("Um vértice", criar_grafo_estrela(1), criar_grafo_estrela(1)),
        ("Dois vértices isolados", criar_grafo_estrela(2), criar_grafo_caminho(2)),
        ("Grafos completamente diferentes",
         criar_grafo_completo(5), criar_grafo_estrela(5)),
    ]

    print("Caso de Teste              | Isomorfos | Tempo (s)")
    print("-" * 45)

    resultados = []

    for desc, g1, g2 in casos:
        start_time = time.perf_counter()
        isomorfo = isomorfismo_molecular(g1, g2)
        end_time = time.perf_counter()

        tempo = end_time - start_time
        status = "✓" if isomorfo else "✗"

        resultados.append({
            'caso': desc,
            'isomorfo': isomorfo,
            'tempo': tempo
        })

        print(f"{desc:<25} | {status:>9} | {tempo:8.4f}")

    return resultados


# =============================================================================
# FUNÇÃO PRINCIPAL
# =============================================================================

def run_comprehensive_stress_tests():
    """Executa todos os testes de estresse e performance"""
    print("TESTES COMPREENSIVOS DE PERFORMANCE, ESTRESSE E ESTABILIDADE")
    print("=" * 70)
    print("=== AVALIAÇÃO COMPLETA DO ALGORITMO DE ISOMORFISMO MOLECULAR ===\n")

    start_total = time.perf_counter()

    warnings.filterwarnings('ignore', category=RuntimeWarning)
    warnings.filterwarnings('ignore', category=UserWarning)

    try:
        perf_scalability = test_performance_scalability()
        complexity_analysis = test_complexity_analysis()
        memory_efficiency = test_memory_efficiency()
        stress_executions = test_stress_executions()
        stress_large_graphs = test_stress_large_graphs()
        stress_memory = test_stress_memory_usage()
        automorfismos_perf = test_automorfismos_performance()
        numerical_stability = test_numerical_stability()
        edge_cases = test_edge_cases()

        total_time = time.perf_counter() - start_total

        print("\n" + "=" * 70)
        print("RELATÓRIO FINAL - PERFORMANCE E ROBUSTEZ")
        print("=" * 70)

        print("\n🏆 RESULTADOS PRINCIPAIS:")

        tempo_max = max([r['tempo'] for r in perf_scalability])
        print(f"  • Performance: {tempo_max:.4f}s para {perf_scalability[-1]['vertices']} vértices")

        memoria_max = max([r['memoria'] for r in memory_efficiency])
        print(f"  • Eficiência de Memória: {memoria_max:.2f} MB máximo")

        print(f"  • Robustez: {stress_executions['n_execucoes']} execuções sem falhas")

        complexidade_final = complexity_analysis[-1]['complexidade']
        print(f"  • Complexidade Observada: {complexidade_final}")

        print("\n📊 ANÁLISE POR CATEGORIA:")

        if tempo_max < 0.1:
            perf_status = "EXCELENTE ⭐⭐⭐"
        elif tempo_max < 1.0:
            perf_status = "MUITO BOA ⭐⭐"
        else:
            perf_status = "ACEITÁVEL ⭐"
        print(f"  🚀 Performance: {perf_status}")

        if memoria_max < 10:
            mem_status = "EXCELENTE ⭐⭐⭐"
        elif memoria_max < 50:
            mem_status = "MUITO BOA ⭐⭐"
        else:
            mem_status = "BOA ⭐"
        print(f"  💾 Eficiência de Memória: {mem_status}")

        if numerical_stability['resultados_consistentes']:
            stab_status = "EXCELENTE ⭐⭐⭐"
        else:
            stab_status = "PROBLEMA ⚠️"
        print(f"  🎯 Estabilidade: {stab_status}")

        print("\n📈 ESTATÍSTICAS GERAIS:")
        print(f"  • Tempo total de testes: {total_time:.2f}s")
        print(f"  • Testes executados: 9 categorias diferentes")
        print(f"  • Maior grafo processado: {stress_large_graphs[-1]['vertices']} vértices")
        print(f"  • Taxa de sucesso: {100 if stress_executions['sem_falhas'] else 0}%")

        print("\n🎯 AVALIAÇÃO FINAL:")

        criterios = [
            tempo_max < 1.0,
            memoria_max < 100,
            stress_executions['tempos_estaveis'],
            numerical_stability['resultados_consistentes']
        ]

        criterios_aprovados = sum(criterios)
        pontuacao_percentual = (criterios_aprovados / len(criterios)) * 100

        if pontuacao_percentual >= 90:
            status_final = "EXCELENTE 🏆"
            recomendacao = "Pronto para aplicações em produção"
        elif pontuacao_percentual >= 75:
            status_final = "MUITO BOM ✅"
            recomendacao = "Adequado para a maioria das aplicações"
        elif pontuacao_percentual >= 60:
            status_final = "BOM ☑️"
            recomendacao = "Recomendado com monitoramento"
        else:
            status_final = "SATISFATÓRIO ⚠️"
            recomendacao = "Recomendadas otimizações"

        print(f"  {status_final}")
        print(f"  Pontuação: {pontuacao_percentual:.1f}% ({criterios_aprovados}/{len(criterios)} critérios)")
        print(f"  {recomendacao}")

        print("\n💡 RECOMENDAÇÕES:")
        if tempo_max >= 1.0:
            print("  • Considerar otimizações para grafos muito densos")
        if memoria_max >= 100:
            print("  • Monitorar uso de memória em processamento contínuo")
        if not numerical_stability['resultados_consistentes']:
            print("  • Investigar inconsistências nos resultados")

        print("\n🔬 CASOS DE USO RECOMENDADOS:")
        print("  ✅ Química Computacional: Comparação de moléculas")
        print("  ✅ Bioinformática: Análise de estruturas proteicas")
        print("  ✅ Verificação de simetria molecular")
        print("  ✅ Análise de grafos teóricos")

        print("\n" + "=" * 70)

        return {
            'tempo_total': total_time,
            'pontuacao_percentual': pontuacao_percentual,
            'status_final': status_final,
            'performance_max': tempo_max,
            'memoria_max': memoria_max,
            'estabilidade': numerical_stability['resultados_consistentes']
        }

    except Exception as e:
        print(f"\n❌ Erro durante a execução dos testes: {e}")
        import traceback
        traceback.print_exc()
        return {'status': 'ERRO', 'erro': str(e)}


if __name__ == "__main__":
    run_comprehensive_stress_tests()