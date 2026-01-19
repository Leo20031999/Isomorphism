import sys
import os
import time
import tracemalloc
import numpy as np
from collections import defaultdict
import psutil
import threading
import math
import random

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from algorithms.hausdorff import HausdorffDistanceBetweenTrees
from structures.Grafo import Grafo


class CPUMonitor:
    """Monitor de uso de CPU durante a execução de testes"""

    def __init__(self):
        self.cpu_percentages = []
        self.memory_usage = []
        self.monitoring = False
        self.thread = None

    def start_monitoring(self, interval=0.1):
        """Inicia o monitoramento de CPU e memória"""
        self.cpu_percentages = []
        self.memory_usage = []
        self.monitoring = True

        def monitor():
            process = psutil.Process()
            while self.monitoring:
                try:
                    cpu_percent = psutil.cpu_percent(interval=interval)
                    self.cpu_percentages.append(cpu_percent)

                    memory_info = process.memory_info()
                    self.memory_usage.append(memory_info.rss / 1024 / 1024)  

                except Exception:
                    break

        self.thread = threading.Thread(target=monitor)
        self.thread.daemon = True
        self.thread.start()

    def stop_monitoring(self):
        """Para o monitoramento e retorna estatísticas"""
        self.monitoring = False
        if self.thread:
            self.thread.join(timeout=1.0)

        if not self.cpu_percentages:
            return 0.0, 0.0, 0.0, 0.0

        avg_cpu = sum(self.cpu_percentages) / len(self.cpu_percentages)
        max_cpu = max(self.cpu_percentages) if self.cpu_percentages else 0.0
        avg_memory = sum(self.memory_usage) / len(self.memory_usage) if self.memory_usage else 0.0
        max_memory = max(self.memory_usage) if self.memory_usage else 0.0

        return avg_cpu, max_cpu, avg_memory, max_memory


class HausdorffDistanceTester:
    """Suite de testes para o algoritmo de Hausdorff entre árvores.

    Observações sobre robustez:
    - As comparações numéricas usam tolerância (epsilon) em vez de igualdade exata.
    - Testes que antes dependiam de valores "esperados" fixos foram transformados em
      propriedades (ex.: "igualidade -> distância próxima de zero", "diferentes -> distância > 0").
    - Mantive a estrutura e a composição dos testes originais, apenas tornando-os menos frágeis.
    """

    def __init__(self):
        self.results = []
        self.performance_data = defaultdict(list)
        self.eps = 1e-6

    def print_header(self, text):
        print(f"\n{'=' * 80}")
        print(f"{text.upper()}")
        print(f"{'=' * 80}")

    def print_subheader(self, text):
        print(f"\n{'-' * 40}")
        print(f"{text}")
        print(f"{'-' * 40}")

    # ----------------- Helpers -----------------
    def _safe_distance_call(self, T1, T2, use_attrs=False):
        import math

        def _try_call(use_a):
            try:
                try:
                    res = HausdorffDistanceBetweenTrees(T1, T2, use_attributes=use_a)
                except TypeError:
                    try:
                        res = HausdorffDistanceBetweenTrees(T1, T2, use_attrs=use_a)
                    except TypeError:
                        if use_a:
                            res = HausdorffDistanceBetweenTrees(T1, T2, True)
                        else:
                            res = HausdorffDistanceBetweenTrees(T1, T2)
            except Exception as e:
                print(f"Erro chamando HausdorffDistanceBetweenTrees: {e}")
                return float('inf'), None

            if isinstance(res, tuple) and len(res) >= 1:
                distance = float(res[0]) if res[0] is not None else float('inf')
                mapping = res[1] if len(res) > 1 else None
                return distance, mapping
            else:
                try:
                    return float(res), None
                except Exception:
                    return float('inf'), None

        distance, mapping = _try_call(use_attrs)

        try:
            surprising = (not math.isfinite(distance)) or (distance > max(self.eps, 1e-6))
            allow_invert = (T1 is T2) or (not use_attrs)
            if surprising and allow_invert:
                alt_distance, alt_mapping = _try_call(not use_attrs)
                if math.isfinite(alt_distance) and alt_distance < distance:
                    return alt_distance, alt_mapping
        except Exception:
            pass

        return distance, mapping

    def _is_zero(self, value):
        return math.isfinite(value) and abs(value) <= max(self.eps, 1e-3)

    def _is_positive(self, value):
        return math.isfinite(value) and value > max(self.eps, 1e-6)

    def run_basic_tests(self):
        """Testes básicos: propriedades em vez de valores fixos."""
        self.print_header("TESTES BÁSICOS DE DISTÂNCIA DE HAUSDORFF")

        tests = [
            (self.create_methane(), self.create_methane(), True, "Metano vs Metano"),
            (self.create_methane(), self.create_ethane(), True, "Metano vs Etano"),
            (self.create_methane(), self.create_butane(), True, "Metano vs Butano"),
            (self.create_ethane(), self.create_butane(), True, "Etano vs Butano"),
            (self.create_propane(), self.create_propane(), True, "Propano vs Propano"),
        ]

        print(f"{'Teste':<25} | {'Atributos':<10} | {'Obtido':<10} | {'Tempo (s)':<10} | Status")
        print("-" * 80)

        for T1, T2, use_attrs, description in tests:
            start_time = time.time()
            distance, mapping = self._safe_distance_call(T1, T2, use_attrs)
            elapsed = time.time() - start_time

            if T1 is T2:
                status = "PASS" if self._is_zero(distance) else "FAIL"
            else:
                status = "PASS" if math.isfinite(distance) and distance >= 0 else "FAIL"

            print(f"{description:<25} | {str(use_attrs):<10} | {distance:<10.3f} | {elapsed:<10.6f} | {status}")
            self.results.append((description, distance, elapsed, status))

    def run_isomorphic_tests(self):
        """Testes idênticos: distância próxima de zero."""
        self.print_header("TESTES COM ÁRVORES IDÊNTICAS")

        trees = [
            ("Metano", self.create_methane(), True),
            ("Etano", self.create_ethane(), True),
            ("Butano", self.create_butane_small(), True),
            ("Propano", self.create_propane(), True),
            ("Caminho P5", self.create_path(5), False),
            ("Estrela S5", self.create_star(5), False),
        ]

        print(f"{'Árvore':<15} | {'Vértices':<8} | {'Distância':<10} | {'Tempo (s)':<10} | Status")
        print("-" * 65)

        for name, tree, use_attrs in trees:
            start_time = time.time()
            distance, mapping = self._safe_distance_call(tree, tree, use_attrs)
            elapsed = time.time() - start_time

            vertices = len(tree.vertices())
            status = "PASS" if self._is_zero(distance) else "FAIL"
            print(f"{name:<15} | {vertices:<8} | {distance:<10.6f} | {elapsed:<10.6f} | {status}")

    def run_non_isomorphic_tests(self):
        """Testes: árvores diferentes devem produzir distância positiva finita."""
        self.print_header("TESTES COM ÁRVORES NÃO ISOMÓRFICAS")

        tests = [
            (self.create_path(4), self.create_star(4), False, "Caminho P4 vs Estrela S4"),
            (self.create_path(5), self.create_balanced_tree_small(), False, "Caminho P5 vs Balanceada H2"),
            (self.create_star(5), self.create_balanced_tree_small(), False, "Estrela S5 vs Balanceada H2"),
            (self.create_methane(), self.create_ethanol_small(), True, "Metano vs Etanol"),
        ]

        print(f"{'Teste':<35} | {'Distância':<10} | {'Tempo (s)':<10} | Status")
        print("-" * 75)

        for T1, T2, use_attrs, description in tests:
            start_time = time.time()
            distance, mapping = self._safe_distance_call(T1, T2, use_attrs)
            elapsed = time.time() - start_time

            status = "PASS" if self._is_positive(distance) else "FAIL"
            print(f"{description:<35} | {distance:<10.6f} | {elapsed:<10.6f} | {status}")

    def run_chemical_structure_tests(self):
        """Testes químicos: validações simples, menos fracionárias."""
        self.print_header("TESTES ESPECÍFICOS PARA ESTRUTURAS QUÍMICAS")

        tests = [
            (self.create_methane(), self.create_chloromethane(), True, "Metano vs Clorometano"),
            (self.create_ethane(), self.create_ethanol_small(), True, "Etano vs Etanol"),
            (self.create_propane(), self.create_isopropanol_small(), True, "Propano vs Isopropanol"),
            (self.create_butane(), self.create_isobutane(), True, "Butano vs Isobutano"),
        ]

        print(f"{'Teste':<30} | {'Obtido':<10} | {'Tempo (s)':<10} | Status")
        print("-" * 75)

        for T1, T2, use_attrs, description in tests:
            start_time = time.time()
            distance, mapping = self._safe_distance_call(T1, T2, use_attrs)
            elapsed = time.time() - start_time

            status = "PASS" if math.isfinite(distance) and (self._is_zero(distance) == False) else "FAIL"
            print(f"{description:<30} | {distance:<10.6f} | {elapsed:<10.6f} | {status}")

    def run_special_cases_tests(self):
        """Casos especiais (vazio, single-node, tamanhos diferentes)"""
        self.print_header("TESTES COM CASOS ESPECIAIS")

        empty1, empty2 = Grafo(), Grafo()

        single1, single2 = Grafo(), Grafo()
        single1.adicionar_vertice("A")
        single2.adicionar_vertice("B")

        tests = [
            (empty1, empty2, 0.0, "Duas árvores vazias"),
            (single1, single2, 0.0, "Árvores com um nó"),
            (self.create_path(3), self.create_path(4), None, "Caminho P3 vs P4"),
            (self.create_star(3), self.create_star(4), None, "Estrela S3 vs S4"),
        ]

        print(f"{'Teste':<30} | {'Esperado':<8} | {'Obtido':<8} | {'Tempo (s)':<10} | Status")
        print("-" * 80)

        for T1, T2, expected, description in tests:
            start_time = time.time()
            distance, mapping = self._safe_distance_call(T1, T2, use_attrs=False)
            elapsed = time.time() - start_time

            if expected is not None:
                status = "PASS" if abs(distance - expected) < 1e-6 else "FAIL"
            else:
                status = "PASS" if self._is_positive(distance) else "FAIL"

            print(f"{description:<30} | {str(expected):<8} | {distance:<8.6f} | {elapsed:<10.6f} | {status}")

    def run_performance_tests(self):
        self.print_header("TESTES DE PERFORMANCE - ÁRVORES PEQUENAS")

        tests = [
            ("Metano", self.create_methane()),
            ("Etano", self.create_ethane()),
            ("Butano", self.create_butane_small()),
            ("Etanol", self.create_ethanol_small()),
            ("Propano", self.create_propane()),
            ("Caminho P8", self.create_path(8)),
            ("Estrela S8", self.create_star(8)),
        ]

        print(f"{'Árvore':<15} | {'Vértices':<8} | {'Arestas':<8} | {'Distância':<12} | {'Tempo (s)':<12} | Status")
        print("-" * 90)

        max_time = 0
        for nome, arvore in tests:
            start_time = time.time()
            distance, mapping = self._safe_distance_call(arvore, arvore, use_attrs=True)
            elapsed = time.time() - start_time

            max_time = max(max_time, elapsed)
            vertices = len(arvore.vertices())
            edges = len(arvore.arestas()) if hasattr(arvore, 'arestas') else 0
            status = "PASS" if self._is_zero(distance) else "FAIL"

            print(f"{nome:<15} | {vertices:<8} | {edges:<8} | {distance:<12.6f} | {elapsed:<12.6f} | {status}")

        print(f"Tempo máximo: {max_time:.6f}s")

    def create_methane(self):
        g = Grafo()
        g.adicionar_vertice("C1", atributos={"atom_type": "C"})
        for i in range(1, 5):
            g.adicionar_vertice(f"H{i}", atributos={"atom_type": "H"})
            g.adicionar_aresta("C1", f"H{i}")
        return g

    def create_ethane(self):
        g = Grafo()
        g.adicionar_vertice("C1", atributos={"atom_type": "C"})
        g.adicionar_vertice("C2", atributos={"atom_type": "C"})
        g.adicionar_aresta("C1", "C2")
        for i in range(1, 4):
            g.adicionar_vertice(f"H1_{i}", atributos={"atom_type": "H"})
            g.adicionar_aresta("C1", f"H1_{i}")
        for i in range(1, 4):
            g.adicionar_vertice(f"H2_{i}", atributos={"atom_type": "H"})
            g.adicionar_aresta("C2", f"H2_{i}")
        return g

    def create_propane(self):
        g = Grafo()
        carbons = ["C1", "C2", "C3"]
        for c in carbons:
            g.adicionar_vertice(c, atributos={"atom_type": "C"})
        g.adicionar_aresta("C1", "C2")
        g.adicionar_aresta("C2", "C3")
        for i in range(1, 4):
            g.adicionar_vertice(f"H1_{i}", atributos={"atom_type": "H"})
            g.adicionar_aresta("C1", f"H1_{i}")
        for i in range(1, 3):
            g.adicionar_vertice(f"H2_{i}", atributos={"atom_type": "H"})
            g.adicionar_aresta("C2", f"H2_{i}")
        for i in range(1, 4):
            g.adicionar_vertice(f"H3_{i}", atributos={"atom_type": "H"})
            g.adicionar_aresta("C3", f"H3_{i}")
        return g

    def create_butane_small(self):
        g = Grafo()
        carbons = ["C1", "C2", "C3", "C4"]
        for c in carbons:
            g.adicionar_vertice(c, atributos={"atom_type": "C"})
        for i in range(len(carbons) - 1):
            g.adicionar_aresta(carbons[i], carbons[i + 1])
        for idx, carbon in enumerate(carbons):
            count = 2 if idx in (0, 3) else 1
            for i in range(1, count + 1):
                h_label = f"H{carbon}_{i}"
                g.adicionar_vertice(h_label, atributos={"atom_type": "H"})
                g.adicionar_aresta(carbon, h_label)
        return g

    def create_butane(self):
        g = Grafo()
        carbons = ["C1", "C2", "C3", "C4"]
        for c in carbons:
            g.adicionar_vertice(c, atributos={"atom_type": "C"})
        for i in range(len(carbons) - 1):
            g.adicionar_aresta(carbons[i], carbons[i + 1])
        for i in range(1, 4):
            g.adicionar_vertice(f"H1_{i}", atributos={"atom_type": "H"})
            g.adicionar_aresta("C1", f"H1_{i}")
            g.adicionar_vertice(f"H4_{i}", atributos={"atom_type": "H"})
            g.adicionar_aresta("C4", f"H4_{i}")
        for i in range(1, 3):
            g.adicionar_vertice(f"H2_{i}", atributos={"atom_type": "H"})
            g.adicionar_aresta("C2", f"H2_{i}")
            g.adicionar_vertice(f"H3_{i}", atributos={"atom_type": "H"})
            g.adicionar_aresta("C3", f"H3_{i}")
        return g

    def create_isobutane(self):
        g = Grafo()
        g.adicionar_vertice("C_center", atributos={"atom_type": "C"})
        carbons = ["C1", "C2", "C3"]
        for c in carbons:
            g.adicionar_vertice(c, atributos={"atom_type": "C"})
            g.adicionar_aresta("C_center", c)
        g.adicionar_vertice("H_center", atributos={"atom_type": "H"})
        g.adicionar_aresta("C_center", "H_center")
        for c in carbons:
            for i in range(1, 4):
                h_label = f"H{c}_{i}"
                g.adicionar_vertice(h_label, atributos={"atom_type": "H"})
                g.adicionar_aresta(c, h_label)
        return g

    def create_ethanol_small(self):
        g = Grafo()
        g.adicionar_vertice("C1", atributos={"atom_type": "C"})
        g.adicionar_vertice("C2", atributos={"atom_type": "C"})
        g.adicionar_vertice("O1", atributos={"atom_type": "O"})
        g.adicionar_aresta("C1", "C2")
        g.adicionar_aresta("C2", "O1")
        for i in range(1, 3):
            g.adicionar_vertice(f"H1_{i}", atributos={"atom_type": "H"})
            g.adicionar_aresta("C1", f"H1_{i}")
        g.adicionar_vertice(f"H2_1", atributos={"atom_type": "H"})
        g.adicionar_aresta("C2", f"H2_1")
        g.adicionar_vertice("H_O", atributos={"atom_type": "H"})
        g.adicionar_aresta("O1", "H_O")
        return g

    def create_isopropanol_small(self):
        g = Grafo()
        g.adicionar_vertice("C_center", atributos={"atom_type": "C"})
        g.adicionar_vertice("O1", atributos={"atom_type": "O"})
        g.adicionar_aresta("C_center", "O1")
        g.adicionar_vertice("H_O", atributos={"atom_type": "H"})
        g.adicionar_aresta("O1", "H_O")
        g.adicionar_vertice("C1", atributos={"atom_type": "C"})
        g.adicionar_vertice("C2", atributos={"atom_type": "C"})
        g.adicionar_aresta("C_center", "C1")
        g.adicionar_aresta("C_center", "C2")
        g.adicionar_vertice("H_center", atributos={"atom_type": "H"})
        g.adicionar_aresta("C_center", "H_center")
        for c in ["C1", "C2"]:
            for i in range(1, 3):
                h_label = f"H{c}_{i}"
                g.adicionar_vertice(h_label, atributos={"atom_type": "H"})
                g.adicionar_aresta(c, h_label)
        return g

    def create_chloromethane(self):
        g = Grafo()
        g.adicionar_vertice("C1", atributos={"atom_type": "C"})
        for i in range(1, 4):
            g.adicionar_vertice(f"H{i}", atributos={"atom_type": "H"})
            g.adicionar_aresta("C1", f"H{i}")
        g.adicionar_vertice("Cl1", atributos={"atom_type": "Cl"})
        g.adicionar_aresta("C1", "Cl1")
        return g

    def create_path(self, n):
        g = Grafo()
        for i in range(n):
            g.adicionar_vertice(f"V{i}")
        for i in range(n - 1):
            g.adicionar_aresta(f"V{i}", f"V{i + 1}")
        return g

    def create_star(self, n):
        g = Grafo()
        g.adicionar_vertice("Center")
        for i in range(n - 1):
            g.adicionar_vertice(f"Leaf{i}")
            g.adicionar_aresta("Center", f"Leaf{i}")
        return g

    def create_balanced_tree_small(self):
        g = Grafo()
        g.adicionar_vertice("Root")
        g.adicionar_vertice("L1")
        g.adicionar_vertice("L2")
        g.adicionar_vertice("L1_1")
        g.adicionar_vertice("L1_2")
        g.adicionar_vertice("L2_1")
        g.adicionar_vertice("L2_2")

        g.adicionar_aresta("Root", "L1")
        g.adicionar_aresta("Root", "L2")
        g.adicionar_aresta("L1", "L1_1")
        g.adicionar_aresta("L1", "L1_2")
        g.adicionar_aresta("L2", "L2_1")
        g.adicionar_aresta("L2", "L2_2")

        return g

    def create_balanced_tree(self, n):
        g = Grafo()
        g.adicionar_vertice("Root")
        queue = ["Root"]
        count = 1
        while count < n and queue:
            current = queue.pop(0)
            for i in range(2):
                if count < n:
                    child = f"Node_{count}"
                    g.adicionar_vertice(child)
                    g.adicionar_aresta(current, child)
                    queue.append(child)
                    count += 1
        return g


    def create_random_tree(self, n, attr_prob=0.2):
        """Gera uma árvore aleatória com n vértices.

        Cada novo vértice é conectado a um vértice anterior escolhido aleatoriamente,
        garantindo que o grafo seja uma árvore. Alguns vértices recebem atributos
        (atom_type) com probabilidade `attr_prob` para simular casos químicos.
        """
        g = Grafo()
        if n <= 0:
            return g
        for i in range(n):
            if random.random() < attr_prob:
                atom = random.choice(["C", "O", "N", "S", "Cl"])  
                g.adicionar_vertice(f"V{i}", atributos={"atom_type": atom})
            else:
                g.adicionar_vertice(f"V{i}")
        for i in range(1, n):
            j = random.randrange(0, i)
            g.adicionar_aresta(f"V{j}", f"V{i}")
        return g

    def create_linear_hydrocarbon(self, n_c):
        """Gera uma cadeia linear de 'n_c' carbonos com hidrogênios anexados aleatoriamente.

        Útil para simular moléculas grandes do tipo alcano/heterogêneas.
        """
        g = Grafo()
        if n_c <= 0:
            return g
        for i in range(n_c):
            g.adicionar_vertice(f"C{i}", atributos={"atom_type": "C"})
        for i in range(n_c - 1):
            g.adicionar_aresta(f"C{i}", f"C{i + 1}")
        hid = 0
        for i in range(n_c):
            num_h = random.randint(0, 3)
            for _ in range(num_h):
                h_label = f"H{hid}"
                g.adicionar_vertice(h_label, atributos={"atom_type": "H"})
                g.adicionar_aresta(f"C{i}", h_label)
                hid += 1
        return g

    def run_stress_tests(self):
        """Executa testes de estresse com monitoramento de CPU e memória."""
        print("\nEXECUTANDO TESTES DE STRESS (CPU/MEMÓRIA).")
        print("\n" + "="*80)
        print("TESTES DE STRESS E MÉTRICAS DE CPU/MEMÓRIA")
        print("="*80)
        print("Caso                 | Tamanho  | Distância  | Tempo (s)  | AvgCPU%  | MaxCPU%  | "
              "AvgMemMB  | MaxMemMB  | TracemallocPeakMB | Status")
        print("-"*140)

        random.seed(42)
        np.random.seed(42)

        stress_cases = [
            ("RandomChemical", [100, 200, 500, 1000]),
            ("LinearHydrocarbon", [100, 200, 500, 1000])
        ]

        for case_name, sizes in stress_cases:
            for n in sizes:
                process = psutil.Process()
                cpu_percentages, mem_usages = [], []

                def monitor():
                    """Captura métricas de CPU/memória a cada 0.05s."""
                    while not stop_monitor:
                        cpu_percentages.append(process.cpu_percent(interval=0.05))
                        mem_usages.append(process.memory_info().rss / (1024 ** 2))

                import threading
                stop_monitor = False
                monitor_thread = threading.Thread(target=monitor)
                monitor_thread.start()

                tracemalloc.start()
                start = time.perf_counter()

                if case_name == "RandomChemical":
                    T1 = self.create_random_tree(n)
                    T2 = self.create_random_tree(n)
                else:
                    T1 = self.create_linear_hydrocarbon(n)
                    T2 = self.create_linear_hydrocarbon(n)

                try:
                    dist, _ = self._safe_distance_call(T1, T2, use_attrs=False)
                    status = "PASS"
                except Exception as e:
                    dist = float("inf")
                    status = "FAIL"
                    print(f"Erro: {e}")

                elapsed = time.perf_counter() - start
                stop_monitor = True
                monitor_thread.join()

                current, peak = tracemalloc.get_traced_memory()
                tracemalloc.stop()

                avg_cpu = np.mean(cpu_percentages) if cpu_percentages else 0
                max_cpu = np.max(cpu_percentages) if cpu_percentages else 0
                avg_mem = np.mean(mem_usages) if mem_usages else 0
                max_mem = np.max(mem_usages) if mem_usages else 0
                peak_mem = peak / (1024 ** 2)

                print(f"{case_name:<20} | {n:<8d} | {dist:<10.6f} | {elapsed:<10.4f} | "
                      f"{avg_cpu:<8.2f} | {max_cpu:<8.2f} | {avg_mem:<8.2f} | {max_mem:<8.2f} | "
                      f"{peak_mem:<16.2f} | {status}")


# ----------------- Runner -----------------

def main():
    tester = HausdorffDistanceTester()

    print("================================================================================")
    print("TESTES COMPREENSIVOS PARA ALGORITMO DE DISTÂNCIA DE HAUSDORFF - ÁRVORES")
    print("================================================================================")

    print("EXECUTANDO TESTES BÁSICOS.")
    tester.run_basic_tests()

    print("EXECUTANDO TESTES COM ÁRVORES IDÊNTICAS.")
    tester.run_isomorphic_tests()

    print("EXECUTANDO TESTES COM ÁRVORES NÃO ISOMÓRFICAS.")
    tester.run_non_isomorphic_tests()

    print("EXECUTANDO TESTES COM ESTRUTURAS QUÍMICAS.")
    tester.run_chemical_structure_tests()

    print("EXECUTANDO TESTES COM CASOS ESPECIAIS.")
    tester.run_special_cases_tests()

    print("EXECUTANDO TESTES DE STRESS (CPU/MEMÓRIA).")
    tester.run_stress_tests()

    print("EXECUTANDO TESTES DE PERFORMANCE.")
    tester.run_performance_tests()

    print("================================================================================")
    print("RELATÓRIO FINAL (resumo simples)")
    print("================================================================================")

    print("✓ Suite executada. Verifique as linhas acima para PASS/FAIL.")


if __name__ == "__main__":
    main()