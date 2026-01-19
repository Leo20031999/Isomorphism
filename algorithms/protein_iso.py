# A new Technique in Protein Structure Quantitative Identification
# (GUO et al., 2022)
import numpy as np
from collections import Counter
from scipy.linalg import svd
from scipy.sparse import csr_matrix
import warnings
from functools import lru_cache

class ProteinGraphDistance:

    def __init__(self, use_labels=True, sensitivity=1.0, use_sparse=False):
        self.use_labels = use_labels
        self.sensitivity = max(0.1, min(2.0, sensitivity))
        self.use_sparse = use_sparse

        self.letter_weights = {
            'A': 1.0, 'B': 2.0, 'C': 3.0, 'D': 4.0, 'E': 5.0, 'F': 6.0, 'G': 7.0,
            'H': 8.0, 'I': 9.0, 'J': 10.0, 'K': 11.0, 'L': 12.0, 'M': 13.0, 'N': 14.0,
            'O': 15.0, 'P': 16.0, 'Q': 17.0, 'R': 18.0, 'S': 19.0, 'T': 20.0,
            'U': 21.0, 'V': 22.0, 'W': 23.0, 'X': 24.0, 'Y': 25.0, 'Z': 26.0
        }

        self.element_weights = {
            'H': 0.3, 'C': 1.0, 'N': 1.4, 'O': 1.6, 'F': 1.8,
            'P': 2.0, 'S': 1.7, 'Cl': 2.2, 'Br': 2.8, 'I': 3.2,
            'AA': 1.2, 'CA': 1.0, 'CB': 1.1, 'CD': 1.1, 'CG': 1.1,
            'CZ': 1.1, 'N': 1.4, 'O': 1.6, 'OG': 1.6, 'OD1': 1.6,
            'OD2': 1.6, 'OE1': 1.6, 'OE2': 1.6, 'OH': 1.6, 'NZ': 1.4,
            'SD': 1.8, 'SG': 1.7, 'NE': 1.4, 'NH1': 1.4, 'NH2': 1.4,
            'ND1': 1.4, 'ND2': 1.4, 'NE2': 1.4, 'OG1': 1.6, 'CE': 1.1
        }

        self.bond_weights = {
            'single': 1.0, 'double': 1.8, 'triple': 2.5,
            'aromatic': 1.4, 'amide': 1.6, 'ionic': 2.0, 'hydrogen': 0.8,
            'covalent': 1.0, 'hydrogen_anti': 0.9, 'hydrogen_para': 0.9,
            'peptide': 1.2, 'disulfide': 1.7, 'backbone': 1.1,
            'sidechain': 1.0, 'salt_bridge': 1.9
        }

    def build_vertex_adjacency_matrix(self, grafo):
        """Matriz de adjacência de vértices - CORRIGIDA PARA USAR RÓTULOS"""
        vertices = sorted(grafo.vertices())
        n = len(vertices)

        if n == 0:
            if self.use_sparse:
                return csr_matrix((0, 0), dtype=float), {}
            return np.zeros((0, 0), dtype=float), {}

        if self.use_sparse:
            vertex_to_index = {v: i for i, v in enumerate(vertices)}
            row, col, data = [], [], []

            for i, u in enumerate(vertices):
                u_label = grafo.get_rotulo_vertice(u) if self.use_labels else None
                neighbors = grafo.vizinhanca(u)

                for v in neighbors:
                    j = vertex_to_index[v]
                    row.append(i)
                    col.append(j)

                    if self.use_labels and u_label:
                        u_weight = self.element_weights.get(u_label, 1.0)
                        v_label = grafo.get_rotulo_vertice(v)
                        v_weight = self.element_weights.get(v_label, 1.0) if v_label else 1.0

                        weight = (u_weight + v_weight) / 2.0

                        if u_label == v_label:
                            weight *= 1.3 
                        else:
                            weight *= 0.7  

                        data.append(weight)
                    else:
                        data.append(1.0)

            matrix = csr_matrix((data, (row, col)), shape=(n, n))
        else:
            matrix = np.zeros((n, n), dtype=float)
            vertex_to_index = {v: i for i, v in enumerate(vertices)}

            for i, u in enumerate(vertices):
                u_label = grafo.get_rotulo_vertice(u) if self.use_labels else None
                neighbors = grafo.vizinhanca(u)

                for v in neighbors:
                    j = vertex_to_index[v]

                    if self.use_labels and u_label:
                
                        u_weight = self.element_weights.get(u_label, 1.0)
                        v_label = grafo.get_rotulo_vertice(v)
                        v_weight = self.element_weights.get(v_label, 1.0) if v_label else 1.0

                        weight = (u_weight + v_weight) / 2.0

                        if u_label == v_label:
                            weight *= 1.3 
                        else:
                            weight *= 0.7 

                        matrix[i][j] = weight
                    else:
                        matrix[i][j] = 1.0

        return matrix, vertex_to_index

    def build_edge_adjacency_matrix(self, grafo):
        """Matriz de adjacência de arestas - CORRIGIDA PARA USAR RÓTULOS"""
        arestas = sorted(grafo.arestas())
        m = len(arestas)

        if m == 0:
            if self.use_sparse:
                return csr_matrix((0, 0), dtype=float), {}
            return np.zeros((0, 0), dtype=float), {}

        if self.use_sparse:
            edge_to_index = {}
            row, col, data = [], [], []

            for idx, (u, v) in enumerate(arestas):
                edge_to_index[(u, v)] = idx
                edge_to_index[(v, u)] = idx

            edge_info = []
            for u, v in arestas:
                edge_label = grafo.get_rotulo_aresta(u, v) if self.use_labels else None
                edge_weight = self.bond_weights.get(edge_label, 1.0) if edge_label else 1.0
                edge_info.append((edge_label, edge_weight))

            for i, (u1, v1) in enumerate(arestas):
                edge1_label, edge1_weight = edge_info[i]

                for j, (u2, v2) in enumerate(arestas):
                    if i == j:
                        continue

                    if u1 == u2 or u1 == v2 or v1 == u2 or v1 == v2:
                        row.append(i)
                        col.append(j)

                        if self.use_labels:
                            edge2_label, edge2_weight = edge_info[j]
                            weight = (edge1_weight + edge2_weight) / 2.0

                            if edge1_label and edge2_label:
                                if edge1_label == edge2_label:
                                    weight *= 1.2  
                                else:
                                    weight *= 0.8  
                            data.append(weight)
                        else:
                            data.append(1.0)

            matrix = csr_matrix((data, (row, col)), shape=(m, m))
        else:
            matrix = np.zeros((m, m), dtype=float)
            edge_to_index = {}

            for idx, (u, v) in enumerate(arestas):
                edge_to_index[(u, v)] = idx
                edge_to_index[(v, u)] = idx

            edge_info = []
            for u, v in arestas:
                edge_label = grafo.get_rotulo_aresta(u, v) if self.use_labels else None
                edge_weight = self.bond_weights.get(edge_label, 1.0) if edge_label else 1.0
                edge_info.append((edge_label, edge_weight))

            for i, (u1, v1) in enumerate(arestas):
                edge1_label, edge1_weight = edge_info[i]

                for j, (u2, v2) in enumerate(arestas):
                    if i == j:
                        continue

                    if u1 == u2 or u1 == v2 or v1 == u2 or v1 == v2:
                        if self.use_labels:
                            edge2_label, edge2_weight = edge_info[j]
                            weight = (edge1_weight + edge2_weight) / 2.0

                            if edge1_label and edge2_label:
                                if edge1_label == edge2_label:
                                    weight *= 1.2 
                                else:
                                    weight *= 0.8 

                            matrix[i][j] = weight
                        else:
                            matrix[i][j] = 1.0

        return matrix, edge_to_index

    def permutation_distance(self, seq1, seq2):
        seq1 = np.asarray(seq1, dtype=float)
        seq2 = np.asarray(seq2, dtype=float)

        if seq1.shape != seq2.shape:
            return 1.0

        n = len(seq1)
        if n == 0:
            return 0.0

        if np.allclose(seq1, seq2, atol=1e-10):
            return 0.0

        max_val = max(np.max(np.abs(seq1)), np.max(np.abs(seq2)), 1e-10)
        seq1_norm = seq1 / max_val
        seq2_norm = seq2 / max_val

        features1 = []
        features2 = []

        num_powers = min(n, 8)
        for power in range(1, num_powers + 1):
            features1.append(np.sum(seq1_norm ** power))
            features2.append(np.sum(seq2_norm ** power))

        if n > 1:
            stats1 = [
                np.mean(seq1_norm),
                np.std(seq1_norm),
                np.median(seq1_norm),
                np.max(seq1_norm),
                np.min(seq1_norm)
            ]
            stats2 = [
                np.mean(seq2_norm),
                np.std(seq2_norm),
                np.median(seq2_norm),
                np.max(seq2_norm),
                np.min(seq2_norm)
            ]
            features1.extend(stats1)
            features2.extend(stats2)

        features1 = np.array(features1)
        features2 = np.array(features2)

        numerator = np.sum((features1 - features2) ** 2)
        denominator = np.sum(features1 ** 2 + features2 ** 2)

        if denominator < 1e-12:
            return 0.0

        distance = np.sqrt(numerator / denominator)
        adjusted_distance = min(distance * self.sensitivity, 1.0)

        return adjusted_distance

    def calculate_row_sums(self, matrix):
        """Somas das linhas - OTIMIZADA PARA MATRIZES ESPARSAS"""
        if hasattr(matrix, 'shape') and 0 in matrix.shape:
            return []

        if hasattr(matrix, 'toarray'):
            return np.array(matrix.sum(axis=1)).flatten().tolist()
        else:
            return np.sum(matrix, axis=1).tolist()

    @lru_cache(maxsize=100)
    def _get_eigenvalues_cached(self, matrix_tuple):
        """Cache para autovalores baseado na matriz (tupla hashable)"""
        matrix = np.array(matrix_tuple)
        return self._compute_eigenvalues(matrix)

    def _compute_eigenvalues(self, matrix):
        """Cálculo interno de autovalores"""
        if matrix.size == 0:
            return []
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                eigenvalues = np.linalg.eigvalsh(matrix.astype(float))
                eigenvalues = [round(x, 6) for x in eigenvalues if abs(x) > 1e-10]
                return sorted(eigenvalues)
        except np.linalg.LinAlgError:
            return []

    def get_eigenvalues(self, matrix):
        """Autovalores - COM CACHE E TRATAMENTO ROBUSTO"""
        if hasattr(matrix, 'shape') and 0 in matrix.shape:
            return []

        if hasattr(matrix, 'toarray'):
            matrix = matrix.toarray()

        if matrix.size <= 10000:  
            matrix_tuple = tuple(matrix.flatten().tolist())
            return self._get_eigenvalues_cached(matrix_tuple)
        else:
            return self._compute_eigenvalues(matrix)

    def convert_to_equinumerous_sequence(self, sequence):
        if len(sequence) == 0:
            return []

        sequence = [round(float(x), 6) for x in sequence]
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

    def equinumerous_distance(self, seq1, seq2):
        """Distância equinumerosa - OTIMIZADA"""
        if len(seq1) == 0 and len(seq2) == 0:
            return 0.0

        if seq1 == seq2:
            return 0.0

        if len(seq1) != len(seq2):
            return 1.0

        equi_seq1 = self.convert_to_equinumerous_sequence(seq1)
        equi_seq2 = self.convert_to_equinumerous_sequence(seq2)

        return self.permutation_distance(equi_seq1, equi_seq2)

    def vector_space_similarity(self, U1, U2):
        """Calcula similaridade entre espaços vetoriais"""
        if U1.size == 0 or U2.size == 0:
            return 0.0

        min_cols = min(U1.shape[1], U2.shape[1])
        U1 = U1[:, :min_cols]
        U2 = U2[:, :min_cols]

        try:
            similarity = np.linalg.norm(U1.T @ U2, ord='fro')
            norm1 = np.linalg.norm(U1, ord='fro')
            norm2 = np.linalg.norm(U2, ord='fro')

            if norm1 * norm2 < 1e-12:
                return 0.0

            return similarity / (norm1 * norm2)
        except:
            return 0.0

    def compute_svd_features(self, matrix):
        """Análise SVD - OTIMIZADA COM VALORES SINGULARES"""
        if hasattr(matrix, 'shape') and 0 in matrix.shape:
            return np.array([]), np.array([]), np.array([])

        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")

                if hasattr(matrix, 'toarray'):
                    matrix = matrix.toarray()

                U, s, Vt = svd(matrix.astype(np.float64), full_matrices=False, lapack_driver='gesvd')

                U = np.nan_to_num(U, nan=0.0, posinf=1.0, neginf=-1.0)
                Vt = np.nan_to_num(Vt, nan=0.0, posinf=1.0, neginf=-1.0)
                s = np.nan_to_num(s, nan=0.0, posinf=1.0, neginf=-1.0)

                return U, s, Vt
        except Exception as e:
            n = min(matrix.shape) if matrix.shape[0] > 0 else 1
            return np.eye(n), np.ones(n), np.eye(n)

    def svd_system_distance(self, U1, s1, Vt1, U2, s2, Vt2):
        s_distance = self.equinumerous_distance(s1.tolist(), s2.tolist())

        left_similarity = self.vector_space_similarity(U1, U2)
        right_similarity = self.vector_space_similarity(Vt1.T, Vt2.T)

        left_distance = 1.0 - left_similarity
        right_distance = 1.0 - right_similarity

        total_distance = (s_distance + left_distance + right_distance) / 3.0
        return min(total_distance, 1.0)

    def quantitative_distance(self, grafo1, grafo2, verbose=True):
        """
        Algoritmo principal - CORRIGIDO PARA VERIFICAÇÃO DE RÓTULOS
        """
        if verbose:
            print(f"=== ALGORITMO DO ARTIGO - INICIADO (use_labels={self.use_labels}) ===")

        try:
            n1, m1 = len(grafo1.vertices()), len(grafo1.arestas())
            n2, m2 = len(grafo2.vertices()), len(grafo2.arestas())

            if n1 != n2 or m1 != m2:
                if verbose:
                    print(f"Grafos de tamanhos diferentes - Distância: 1.0")
                return 1.0

            if n1 == 0 and n2 == 0:
                return 0.0

            vertices1 = sorted(grafo1.vertices())
            vertices2 = sorted(grafo2.vertices())
            arestas1 = sorted(grafo1.arestas())
            arestas2 = sorted(grafo2.arestas())

            if vertices1 == vertices2 and arestas1 == arestas2:
                if not self.use_labels:
                    if verbose:
                        print("Grafos estruturalmente idênticos (sem considerar rótulos) - Distância: 0.0")
                    return 0.0
                else:
                    vertex_labels_match = True
                    for v in vertices1:
                        if grafo1.get_rotulo_vertice(v) != grafo2.get_rotulo_vertice(v):
                            vertex_labels_match = False
                            break

                    edge_labels_match = True
                    if vertex_labels_match:
                        for u, v in arestas1:
                            if grafo1.get_rotulo_aresta(u, v) != grafo2.get_rotulo_aresta(u, v):
                                edge_labels_match = False
                                break

                    if vertex_labels_match and edge_labels_match:
                        if verbose:
                            print("Grafos completamente idênticos (incluindo rótulos) - Distância: 0.0")
                        return 0.0
                    else:
                        if verbose:
                            print(
                                "Grafos estruturalmente idênticos mas com rótulos diferentes - Calculando distância...")
            else:
                if verbose:
                    print("Grafos estruturalmente diferentes - Calculando distância...")

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

            if verbose:
                print(f"Distancia vertices (s1): {s1:.6f}")
                print(f"Distancia arestas (s2): {s2:.6f}")

            if s1 != 0 or s2 != 0:
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

            if s3 != 0 or s4 != 0:
                distance = 0.25 * s1 + 0.25 * s2 + 0.25 * s3 + 0.25 * s4
                final_distance = min(float(distance), 1.0)
                if verbose:
                    print(f"Parada na Etapa 6 - Distancia final: {final_distance:.6f}")
                return final_distance

            U_v1, s_v1, Vt_v1 = self.compute_svd_features(vertex_adj1)
            U_v2, s_v2, Vt_v2 = self.compute_svd_features(vertex_adj2)
            U_e1, s_e1, Vt_e1 = self.compute_svd_features(edge_adj1)
            U_e2, s_e2, Vt_e2 = self.compute_svd_features(edge_adj2)

            s5 = self.svd_system_distance(U_v1, s_v1, Vt_v1, U_v2, s_v2, Vt_v2)
            s6 = self.svd_system_distance(U_e1, s_e1, Vt_e1, U_e2, s_e2, Vt_e2)

            if verbose:
                print(f"Distancia SVD vertices (s5): {s5:.6f}")
                print(f"Distancia SVD arestas (s6): {s6:.6f}")

            if s5 != 0 or s6 != 0:
                distance = (0.25 * s1 + 0.25 * s2 + 0.25 * s3 + 0.25 * s4 +
                            0.25 * s5 + 0.25 * s6)
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

    def profile_computation(self, grafo1, grafo2):
        """Método para profiling do cálculo de distância"""
        import time

        start_time = time.time()

        original_sparse = self.use_sparse
        self.use_sparse = False

        result = self.quantitative_distance(grafo1, grafo2, verbose=False)

        self.use_sparse = original_sparse

        end_time = time.time()

        print(f"Tempo de computação: {end_time - start_time:.4f} segundos")
        print(f"Resultado: {result:.6f}")

        return result