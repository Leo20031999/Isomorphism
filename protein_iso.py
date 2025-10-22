# A new Technique in Protein Structure Quantitative Identification
# (GUO et al., 2022)
import numpy as np
from collections import Counter
from scipy.linalg import svd
import warnings


class ProteinGraphDistance:

    def __init__(self, use_labels=False, sensitivity=1.0):
        self.use_labels = use_labels
        self.sensitivity = max(0.1, min(2.0, sensitivity))

    def build_vertex_adjacency_matrix(self, grafo):
        """Matriz de adjacência de vértices"""
        vertices = sorted(grafo.vertices())
        n = len(vertices)

        if n == 0:
            return np.zeros((0, 0), dtype=float), {}

        matrix = np.zeros((n, n), dtype=float)
        vertex_to_index = {v: i for i, v in enumerate(vertices)}

        for i, u in enumerate(vertices):
            for j, v in enumerate(vertices):
                if v in grafo.vizinhanca(u):
                    matrix[i][j] = 1.0
        return matrix, vertex_to_index

    def build_edge_adjacency_matrix(self, grafo):
        """Matriz de adjacência de arestas"""
        arestas = sorted(grafo.arestas())
        m = len(arestas)

        if m == 0:
            return np.zeros((0, 0), dtype=float), {}

        matrix = np.zeros((m, m), dtype=float)
        edge_to_index = {}

        for idx, (u, v) in enumerate(arestas):
            edge_to_index[(u, v)] = idx
            edge_to_index[(v, u)] = idx

        for i, (u1, v1) in enumerate(arestas):
            for j, (u2, v2) in enumerate(arestas):
                if i == j:
                    continue
                if u1 == u2 or u1 == v2 or v1 == u2 or v1 == v2:
                    matrix[i][j] = 1.0
        return matrix, edge_to_index

    def permutation_distance(self, seq1, seq2):
        """Distância de permutação com melhor granularidade"""
        if len(seq1) != len(seq2):
            return 1.0

        n = len(seq1)
        if n == 0:
            return 0.0

        if seq1 == seq2:
            return 0.0

        try:
            num_features = min(n, 8)

            power_sums1 = []
            power_sums2 = []

            for power in range(1, num_features + 1):
                sum1 = sum(float(x) ** power for x in seq1)
                sum2 = sum(float(x) ** power for x in seq2)
                power_sums1.append(sum1)
                power_sums2.append(sum2)

            if n > 1:
                features1 = [
                    np.mean(seq1), np.std(seq1),
                    np.median(seq1), max(seq1), min(seq1),
                    np.percentile(seq1, 25), np.percentile(seq1, 75)
                ]
                features2 = [
                    np.mean(seq2), np.std(seq2),
                    np.median(seq2), max(seq2), min(seq2),
                    np.percentile(seq2, 25), np.percentile(seq2, 75)
                ]

                power_sums1.extend(features1[:min(5, len(features1))])
                power_sums2.extend(features2[:min(5, len(features2))])

            max_vals = [max(abs(a), abs(b), 1e-10) for a, b in zip(power_sums1, power_sums2)]
            features1_norm = [a / max_val for a, max_val in zip(power_sums1, max_vals)]
            features2_norm = [b / max_val for b, max_val in zip(power_sums2, max_vals)]

            squared_diff = sum((a - b) ** 2 for a, b in zip(features1_norm, features2_norm))
            distance = np.sqrt(squared_diff / len(features1_norm))

            adjusted_distance = min(distance * self.sensitivity, 1.0)
            return adjusted_distance

        except:
            return 1.0

    def calculate_row_sums(self, matrix):
        """Somas das linhas"""
        if matrix.size == 0:
            return []
        try:
            return [float(sum(row)) for row in matrix]
        except:
            return []

    def get_eigenvalues(self, matrix):
        """Autovalores"""
        if matrix.size == 0:
            return []
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                eigenvalues = np.linalg.eigvals(matrix.astype(complex))
                real_eigenvalues = sorted([round(x.real, 6) for x in eigenvalues])
                return real_eigenvalues
        except:
            return []

    def convert_to_equinumerous_sequence(self, sequence):
        """Conversão equinumerosa"""
        if len(sequence) == 0:
            return []

        try:
            sequence = [round(float(x), 4) for x in sequence]
            counter = Counter(sequence)

            value_mapping = {}
            current_label = 21

            for value in sorted(set(sequence)):
                count = counter[value]
                if count == 1:
                    value_mapping[value] = [float(value)]
                else:
                    labels = [float(current_label + i) for i in range(count)]
                    value_mapping[value] = labels
                    current_label += count

            result = []
            usage_count = {value: 0 for value in set(sequence)}

            for value in sequence:
                result.append(value_mapping[value][usage_count[value]])
                usage_count[value] += 1

            return result

        except:
            return sequence

    def equinumerous_distance(self, seq1, seq2):
        """Distância equinumerosa"""
        if len(seq1) == 0 and len(seq2) == 0:
            return 0.0

        try:
            if seq1 == seq2:
                return 0.0

            equi_seq1 = self.convert_to_equinumerous_sequence(seq1)
            equi_seq2 = self.convert_to_equinumerous_sequence(seq2)

            if equi_seq1 == equi_seq2:
                return 0.0

            return self.permutation_distance(equi_seq1, equi_seq2)
        except:
            return 1.0

    def compute_svd_features(self, matrix):
        """Análise SVD"""
        if matrix.size == 0:
            return np.array([]), np.array([])

        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                matrix_float = matrix.astype(np.float64)
                U, s, Vt = svd(matrix_float, full_matrices=False)

                U = np.nan_to_num(U)
                Vt = np.nan_to_num(Vt)

                return U.flatten(), Vt.flatten()
        except:
            n = min(matrix.shape) if matrix.shape[0] > 0 else 1
            return np.ones(n), np.ones(n)

    def svd_system_distance(self, features1, features2):
        """Distância SVD"""
        if len(features1) == 0 and len(features2) == 0:
            return 0.0

        try:
            f1 = [float(x) for x in features1]
            f2 = [float(x) for x in features2]

            if all(abs(a - b) < 1e-10 for a, b in zip(f1, f2)):
                return 0.0

            return self.equinumerous_distance(f1, f2)
        except:
            return 1.0

    def quantitative_distance(self, grafo1, grafo2, verbose=True):
        """
        Algoritmo principal com melhor granularidade
        """
        if verbose:
            print("=== ALGORITMO DO ARTIGO - INICIADO ===")

        try:
            n1, m1 = len(grafo1.vertices()), len(grafo1.arestas())
            n2, m2 = len(grafo2.vertices()), len(grafo2.arestas())

            if grafo1.vertices() == grafo2.vertices() and grafo1.arestas() == grafo2.arestas():
                if verbose:
                    print("Grafos idênticos - Distância: 0.0")
                return 0.0

            if n1 != n2 or m1 != m2:
                size_penalty = min(abs(n1 - n2) / max(n1, n2, 1) + abs(m1 - m2) / max(m1, m2, 1), 1.0)
                if size_penalty > 0.9:
                    return 1.0


            if n1 == 0 and n2 == 0:
                return 0.0

            vertex_adj1, _ = self.build_vertex_adjacency_matrix(grafo1)
            vertex_adj2, _ = self.build_vertex_adjacency_matrix(grafo2)
            edge_adj1, _ = self.build_edge_adjacency_matrix(grafo1)
            edge_adj2, _ = self.build_edge_adjacency_matrix(grafo2)

            vertex_row_sums1 = self.calculate_row_sums(vertex_adj1)
            vertex_row_sums2 = self.calculate_row_sums(vertex_adj2)
            s1 = self.permutation_distance(vertex_row_sums1, vertex_row_sums2)

            edge_row_sums1 = self.calculate_row_sums(edge_adj1)
            edge_row_sums2 = self.calculate_row_sums(edge_adj2)
            s2 = self.permutation_distance(edge_row_sums1, edge_row_sums2)

            limiar = max(1e-6, 0.1 / max(n1, n2, 1))

            if verbose:
                print(f"Distancia vertices (s1): {s1:.6f}")
                print(f"Distancia arestas (s2): {s2:.6f}")
                print(f"Limiar adaptativo: {limiar:.6f}")

            if s1 > limiar or s2 > limiar:
                distance = 0.5 * s1 + 0.5 * s2
                final_distance = min(float(distance), 1.0)
                if verbose:
                    print(f"Parada na Etapa 3 - Distancia final: {final_distance:.6f}")
                return final_distance

            vertex_eigen1 = self.get_eigenvalues(vertex_adj1)
            vertex_eigen2 = self.get_eigenvalues(vertex_adj2)
            s3 = self.equinumerous_distance(vertex_eigen1, vertex_eigen2)

            edge_eigen1 = self.get_eigenvalues(edge_adj1)
            edge_eigen2 = self.get_eigenvalues(edge_adj2)
            s4 = self.equinumerous_distance(edge_eigen1, edge_eigen2)

            if verbose:
                print(f"Distancia autovalores vertices (s3): {s3:.6f}")
                print(f"Distancia autovalores arestas (s4): {s4:.6f}")

            if s3 > limiar or s4 > limiar:
                distance = 0.25 * s1 + 0.25 * s2 + 0.25 * s3 + 0.25 * s4
                final_distance = min(float(distance), 1.0)
                if verbose:
                    print(f"Parada na Etapa 6 - Distancia final: {final_distance:.6f}")
                return final_distance

            left_v1, right_v1 = self.compute_svd_features(vertex_adj1)
            left_v2, right_v2 = self.compute_svd_features(vertex_adj2)
            left_e1, right_e1 = self.compute_svd_features(edge_adj1)
            left_e2, right_e2 = self.compute_svd_features(edge_adj2)

            s5 = self.svd_system_distance(left_v1, left_v2)
            s6 = self.svd_system_distance(right_v1, right_v2)
            s7 = self.svd_system_distance(left_e1, left_e2)
            s8 = self.svd_system_distance(right_e1, right_e2)

            if verbose:
                print(f"Distancia SVD vertices: {s5:.6f}, {s6:.6f}")
                print(f"Distancia SVD arestas: {s7:.6f}, {s8:.6f}")

            if s5 > limiar or s6 > limiar or s7 > limiar or s8 > limiar:
                distance = (0.125 * s1 + 0.125 * s2 + 0.125 * s3 + 0.125 * s4 +
                            0.125 * s5 + 0.125 * s6 + 0.125 * s7 + 0.125 * s8)
                final_distance = min(float(distance), 1.0)
                if verbose:
                    print(f"Parada na Etapa 13 - Distancia final: {final_distance:.6f}")
                return final_distance

            if verbose:
                print("Grafos são ISOMORFOS!")
            return 0.0

        except Exception as e:
            if verbose:
                print(f"Erro durante cálculo: {e}")
            return 1.0
