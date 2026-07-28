# A Polynomial-time Maximum Common Subgraph Algorithm for Outerplanar Graphs and its Application to Chemoinformatics
# SCHIETGAT; RAMON; BRUYNOOGHE, 2013

import time
from collections import deque, OrderedDict
from typing import Dict, List, Tuple, Optional, Any, Union
import numpy as np
from scipy.optimize import linear_sum_assignment
import logging
import sys
import hashlib

from structures.Grafo import Grafo

# Configurações padrão
DEFAULT_LABEL_WEIGHTS = {
    'H': 0.3, 'C': 1.0, 'N': 1.4, 'O': 1.6, 'F': 1.8,
    'P': 2.0, 'S': 1.7, 'Cl': 2.2, 'Br': 2.8, 'I': 3.2,
    'single': 1.0, 'double': 1.8, 'triple': 2.5,
    'aromatic': 1.4, 'amide': 1.6, 'ionic': 2.0, 'hydrogen': 0.8
}

MAX_VERTICES = 10000
TIMEOUT_DURATION = 30


class SmartCache:
    """Cache inteligente com TTL e limpeza automática"""

    def __init__(self, max_size=1000, ttl=300):
        self.cache = OrderedDict()
        self.max_size = max_size
        self.ttl = ttl
        self.access_times = {}

    def get(self, key):
        if key in self.cache:
            if time.time() - self.access_times.get(key, 0) > self.ttl:
                del self.cache[key]
                del self.access_times[key]
                return None
            self.cache.move_to_end(key)
            return self.cache[key]
        return None

    def set(self, key, value):
        if key in self.cache:
            self.cache.move_to_end(key)
        else:
            if len(self.cache) >= self.max_size:
                oldest_key = next(iter(self.cache))
                del self.cache[oldest_key]
                del self.access_times[oldest_key]
        self.cache[key] = value
        self.access_times[key] = time.time()

    def clear_expired(self):
        current_time = time.time()
        expired_keys = [k for k, t in self.access_times.items()
                        if current_time - t > self.ttl]
        for key in expired_keys:
            if key in self.cache:
                del self.cache[key]
            if key in self.access_times:
                del self.access_times[key]

    def __len__(self):
        return len(self.cache)


class GraphPart:
    """Estrutura para representar partes do grafo"""

    def __init__(self, graph: Grafo, root: Any, is_basic: bool, part_type: str,
                 orientation: Optional[str] = None, split_edge: Optional[Tuple] = None,
                 parent_root: Optional[Any] = None):
        self.graph = graph
        self.root = root
        self.is_basic = is_basic
        self.type = part_type
        self.orientation = orientation
        self.split_edge = split_edge
        self.parent_root = parent_root

    def __str__(self):
        return f"GraphPart(root={self.root}, type={self.type}, basic={self.is_basic}, orientation={self.orientation})"

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        if not isinstance(other, GraphPart):
            return False
        return (self.root == other.root and self.type == other.type and
                self.is_basic == other.is_basic and self.orientation == other.orientation)

    def __hash__(self):
        return hash((self.root, self.type, self.is_basic, self.orientation))


class RMCSRequest:
    """Estrutura para representar uma requisição de cálculo RMCS"""

    def __init__(self, P: GraphPart, Q: GraphPart, callback=None):
        self.P = P
        self.Q = Q
        self.callback = callback
        self.result = None
        self.processed = False


class OuterplanarMCS:
    """Implementação completamente iterativa do algoritmo MCS para grafos outerplanares"""

    def __init__(self, label_weights: Optional[Dict[str, float]] = None,
                 max_vertices: int = MAX_VERTICES, timeout: int = TIMEOUT_DURATION,
                 verbose: bool = False, cache_size: int = 2000):

        self.memo = SmartCache(cache_size)
        self.decomposition_cache = SmartCache(cache_size // 2)
        self.label_weights = label_weights if label_weights is not None else DEFAULT_LABEL_WEIGHTS
        self.max_vertices = max_vertices
        self.timeout = timeout
        self.start_time = None
        self.mcs_best_candidate = (Grafo(), 0.0)
        self._timeout_occurred = False
        self.verbose = verbose
        self.max_root_candidates = 12

        # Estatísticas
        self.stats = {
            'calls_rmcs': 0,
            'cache_hits': 0,
            'decompositions': 0,
            'matching_calls': 0,
            'timeout_occurred': False,
            'iterations': 0
        }

        # Sistema de processamento iterativo
        self.request_queue = deque()
        self.results_cache = {}
        self.dependencies = {}
        self.ready_requests = deque()

        # Configurar logging
        self.logger = self._setup_logger()

    def _setup_logger(self):
        """Configura sistema de logging"""
        logger = logging.getLogger('OuterplanarMCS')
        if self.verbose:
            logger.setLevel(logging.INFO)
            if not logger.handlers:
                handler = logging.StreamHandler()
                formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
                handler.setFormatter(formatter)
                logger.addHandler(handler)
        else:
            logger.setLevel(logging.WARNING)
        return logger

    def _validate_graph_methods(self, graph: Grafo) -> bool:
        """Valida se o grafo possui todos os métodos necessários"""
        required_methods = [
            'vertices', 'arestas', 'get_rotulo_vertice', 'get_rotulo_aresta',
            'vizinhanca', 'grau', 'encontrar_pontes', 'encontrar_blocos',
            'copy', 'adicionar_vertice', 'adicionar_aresta', 'existe_vertice'
        ]

        for method in required_methods:
            if not hasattr(graph, method) or not callable(getattr(graph, method)):
                self.logger.warning(f"Método '{method}' não encontrado na classe Grafo")
                return False
        return True

    def _graph_signature(self, graph: Grafo) -> str:
        """Cria uma assinatura única baseada nas características do grafo"""
        try:
            vertices = sorted(graph.vertices())
            edges = sorted(graph.arestas())

            vertex_data = []
            for v in vertices:
                try:
                    label = graph.get_rotulo_vertice(v) or "None"
                    degree = graph.grau(v)
                    vertex_data.append(f"{v}:{label}:{degree}")
                except:
                    vertex_data.append(f"{v}:Unknown")

            edge_data = []
            for u, v in edges:
                try:
                    label = graph.get_rotulo_aresta(u, v) or "None"
                    edge_data.append(f"{u}-{v}:{label}")
                except:
                    edge_data.append(f"{u}-{v}:Unknown")

            signature_str = f"V[{','.join(vertex_data)}]E[{','.join(edge_data)}]"
            return hashlib.md5(signature_str.encode()).hexdigest()
        except Exception as e:
            self.logger.warning(f"Erro ao criar assinatura do grafo: {e}")
            return f"graph_{id(graph)}"

    def _find_promising_pairs(self, G: Grafo, H: Grafo, candidates_G: List, candidates_H: List) -> List[
        Tuple]:
        """Encontra pares promissores com critérios mais inteligentes"""
        promising_pairs = []

        # Primeiro: pares com mesma label e alta conectividade
        for r in candidates_G[:15]:
            for s in candidates_H[:15]:
                try:
                    if G.get_rotulo_vertice(r) == H.get_rotulo_vertice(s):
                        score = G.grau(r) + H.grau(s)
                        promising_pairs.append((r, s, score))
                except:
                    continue

        # Se não encontrou suficientes, relaxar critério de labels
        if len(promising_pairs) < 5:
            for r in candidates_G[:10]:
                for s in candidates_H[:10]:
                    try:
                        # Mesmo com labels diferentes, se grau similar
                        if abs(G.grau(r) - H.grau(s)) <= 2:
                            score = min(G.grau(r), H.grau(s))
                            promising_pairs.append((r, s, score))
                    except:
                        continue

        # Ordenar por score
        promising_pairs.sort(key=lambda x: x[2], reverse=True)
        return promising_pairs

    def compute_mcs(self, G: Grafo, H: Grafo) -> Tuple[Grafo, float]:
        """Algoritmo principal completamente iterativo - VERSÃO OTIMIZADA"""
        self.start_time = time.time()
        self._timeout_occurred = False
        self.memo.clear_expired()
        self.decomposition_cache.clear_expired()
        self.mcs_best_candidate = (Grafo(), 0.0)

        # Reiniciar estatísticas
        self.stats = {
            'calls_rmcs': 0, 'cache_hits': 0, 'decompositions': 0,
            'matching_calls': 0, 'timeout_occurred': False,
            'iterations': 0, 'start_time': self.start_time
        }

        self.request_queue.clear()
        self.results_cache.clear()
        self.dependencies.clear()
        self.ready_requests.clear()

        self.logger.info(f"Computando MCS entre grafos com {len(G.vertices())} e {len(H.vertices())} vértices")

        # Validação de métodos
        if not self._validate_graph_methods(G) or not self._validate_graph_methods(H):
            self.logger.error("Grafos não possuem todos os métodos necessários")
            return self._compute_enhanced_fallback(G, H)

        # Casos especiais otimizados
        if len(G.vertices()) == 0 or len(H.vertices()) == 0:
            return (Grafo(), 0.0)

        if len(G.vertices()) == 1 and len(H.vertices()) == 1:
            return self._handle_single_vertex_case(G, H)

        if self._are_graphs_identical(G, H):
            return (G.copy(), self._calculate_graph_size(G))

        # Verificação de tamanho com fallback mais inteligente
        total_vertices = len(G.vertices()) + len(H.vertices())
        if total_vertices > 50:  # Limite mais conservador para grafos grandes
            self.logger.info("Grafos muito grandes, usando fallback otimizado")
            return self._compute_enhanced_fallback(G, H)

        # Verificação de conectividade melhorada
        is_connected_G = self._is_connected(G)
        is_connected_H = self._is_connected(H)

        if not is_connected_G or not is_connected_H:
            self.logger.info("Grafos desconexos, usando aproximação para componentes conectados")
            return self._compute_connected_components_mcs(G, H)

        root_candidates_G = self._find_root_candidates(G)
        root_candidates_H = self._find_root_candidates(H)

        self.logger.info(f"Encontrados {len(root_candidates_G)} e {len(root_candidates_H)} candidatos a raiz")

        # Gerar pares promissores com critérios mais flexíveis
        promising_pairs = self._find_promising_pairs(G, H, root_candidates_G, root_candidates_H)

        if not promising_pairs:
            self.logger.info("Nenhum par promissor encontrado, usando fallback avançado")
            return self._compute_enhanced_fallback(G, H)

        # Processar pares promissores com timeout individual
        best_results = []
        for i, (r, s, score) in enumerate(promising_pairs[:8]):  # Limitar a 8 melhores pares
            if self._check_timeout():
                break

            # Timeout individual por par (máximo 5 segundos por par)
            pair_start_time = time.time()

            try:
                parts_G = self._decompose_graph_fast(G, r)
                parts_H = self._decompose_graph_fast(H, s)

                root_parts_G = [p for p in parts_G if p.root == r and p.type == 'BPS']
                root_parts_H = [p for p in parts_H if p.root == s and p.type == 'BPS']

                for root_part_G in root_parts_G[:2]:  # Limitar partes raiz
                    for root_part_H in root_parts_H[:2]:
                        if time.time() - pair_start_time > 5: 
                            break

                        if self._quick_similarity_check([root_part_G], [root_part_H]):
                            result = self._rmcs_with_bbp_fast(root_part_G, root_part_H)
                            if result[1] > self.mcs_best_candidate[1]:
                                self.mcs_best_candidate = result
                            best_results.append(result)
            except Exception as e:
                self.logger.debug(f"Erro ao processar par ({r}, {s}): {e}")
                continue

        # Fallback inteligente se necessário
        if self.mcs_best_candidate[1] == 0.0:
            return self._compute_enhanced_fallback(G, H)

        self.logger.info(f"MCS encontrado com tamanho: {self.mcs_best_candidate[1]}")
        return self.mcs_best_candidate

    def _compute_connected_components_mcs(self, G: Grafo, H: Grafo) -> Tuple[Grafo, float]:
        """MCS para componentes conectados"""
        try:
            components_G = self._find_connected_components(G)
            components_H = self._find_connected_components(H)

            best_mcs = Grafo()
            best_size = 0.0

            # Ordenar componentes por tamanho (maiores primeiro)
            components_G.sort(key=len, reverse=True)
            components_H.sort(key=len, reverse=True)

            # Comparar apenas os componentes maiores para performance
            max_components_to_check = min(3, len(components_G), len(components_H))

            for i in range(max_components_to_check):
                if self._check_timeout():
                    break

                comp_G = components_G[i]
                comp_H = components_H[i]

                subgraph_G = self._create_subgraph(G, comp_G)
                subgraph_H = self._create_subgraph(H, comp_H)

                # Usar fallback para componentes muito grandes
                if len(comp_G) > 15 or len(comp_H) > 15:
                    mcs_result, size = self._compute_enhanced_fallback(subgraph_G, subgraph_H)
                else:
                    mcs_result, size = self.compute_mcs(subgraph_G, subgraph_H)

                if size > best_size:
                    best_mcs = mcs_result
                    best_size = size

            return (best_mcs, best_size)
        except Exception as e:
            self.logger.warning(f"Erro no MCS para componentes conectados: {e}")
            return self._compute_enhanced_fallback(G, H)

    def _decompose_graph_fast(self, graph: Grafo, root: Any) -> List[GraphPart]:
        """Decomposição rápida do grafo com limites rigorosos"""
        cache_key = (self._graph_signature(graph), root, "fast")
        cached_result = self.decomposition_cache.get(cache_key)
        if cached_result is not None:
            return cached_result

        is_basic = self._is_basic_root(graph, root)
        root_part = GraphPart(graph, root, is_basic, 'BPS')

        parts = [root_part]

        # Decomposição limitada para performance
        if len(graph.vertices()) <= 10:
            queue = deque([root_part])
            visited_parts = set([str(root_part)])

            while queue and len(parts) < 20:  # Limite rigoroso
                current = queue.popleft()
                elementary_parts = self._get_elementary_parts(current)

                for part in elementary_parts[:3]:  # Limitar partes elementares
                    part_str = str(part)
                    if part_str not in visited_parts:
                        parts.append(part)
                        visited_parts.add(part_str)
                        if len(parts) < 20:
                            queue.append(part)

        self.decomposition_cache.set(cache_key, parts)
        return parts

    def _compute_connected_components_mcs(self, G: Grafo, H: Grafo) -> Tuple[Grafo, float]:
        """Calcula MCS para grafos desconexos considerando componentes conectados"""
        try:
            # Encontrar componentes conectados
            components_G = self._find_connected_components(G)
            components_H = self._find_connected_components(H)

            best_mcs = Grafo()
            best_size = 0.0

            # Comparar todos os pares de componentes
            for comp_G in components_G:
                for comp_H in components_H:
                    if self._check_timeout():
                        break

                    # Criar subgrafos para os componentes
                    subgraph_G = self._create_subgraph(G, comp_G)
                    subgraph_H = self._create_subgraph(H, comp_H)

                    # Calcular MCS para os subgrafos
                    mcs_result, size = self.compute_mcs(subgraph_G, subgraph_H)

                    if size > best_size:
                        best_mcs = mcs_result
                        best_size = size

            return (best_mcs, best_size)
        except Exception as e:
            self.logger.warning(f"Erro no MCS para componentes conectados: {e}")
            return self._compute_approximate_mcs(G, H)

    def _find_connected_components(self, graph: Grafo) -> List[List[Any]]:
        """Encontra componentes conectados do grafo"""
        visited = set()
        components = []

        for vertex in graph.vertices():
            if vertex not in visited:
                component = []
                queue = deque([vertex])

                while queue:
                    current = queue.popleft()
                    if current not in visited:
                        visited.add(current)
                        component.append(current)
                        try:
                            neighbors = graph.vizinhanca(current)
                            for neighbor in neighbors:
                                if neighbor not in visited:
                                    queue.append(neighbor)
                        except:
                            continue

                if component:
                    components.append(component)

        return components

    def _create_subgraph(self, graph: Grafo, vertices: List[Any]) -> Grafo:
        """Cria um subgrafo contendo apenas os vértices especificados"""
        subgraph = Grafo()
        vertex_set = set(vertices)

        # Adicionar vértices
        for v in vertices:
            try:
                label = graph.get_rotulo_vertice(v)
                subgraph.adicionar_vertice(v, label)
            except:
                continue

        # Adicionar arestas
        for u, v in graph.arestas():
            if u in vertex_set and v in vertex_set:
                try:
                    label = graph.get_rotulo_aresta(u, v)
                    subgraph.adicionar_aresta(u, v, label)
                except:
                    continue

        return subgraph

    def _handle_single_vertex_case(self, G: Grafo, H: Grafo) -> Tuple[Grafo, float]:
        """Lida com o caso de grafos com um único vértice"""
        if len(G.vertices()) == 1 and len(H.vertices()) == 1:
            v_g = G.vertices()[0]
            v_h = H.vertices()[0]

            try:
                if G.get_rotulo_vertice(v_g) == H.get_rotulo_vertice(v_h):
                    graph = Grafo()
                    graph.adicionar_vertice(v_g, G.get_rotulo_vertice(v_g))
                    size = self._calculate_graph_size(graph)
                    return (graph, size)
            except:
                pass

        return (Grafo(), 0.0)

    def _create_cache_key(self, P: GraphPart, Q: GraphPart) -> tuple:
        """Cria chave única para cache de forma robusta e eficiente"""
        try:
            # Estratégia mais simples e direta para evitar recursão
            # Usar identificadores únicos dos grafos + características essenciais
            graph_id_P = f"G{id(P.graph)}_{hash(str(P.graph.vertices()))}"
            graph_id_Q = f"H{id(Q.graph)}_{hash(str(Q.graph.vertices()))}"

            # Características essenciais das partes
            key_parts = (
                graph_id_P,
                str(P.root)[:50],  # Limitar tamanho para evitar chaves muito longas
                str(P.is_basic),
                str(P.type)[:10],
                str(P.orientation)[:15] if P.orientation else "None",
                graph_id_Q,
                str(Q.root)[:50],
                str(Q.is_basic),
                str(Q.type)[:10],
                str(Q.orientation)[:15] if Q.orientation else "None"
            )

            # Converter para tupla hashable
            cache_key = tuple(key_parts)

            # Verificação de segurança
            try:
                hash(cache_key)
            except:
                # Fallback extremamente simples
                return (f"P_{id(P)}_{id(P.graph)}", f"Q_{id(Q)}_{id(Q.graph)}")

            return cache_key

        except Exception as e:
            # Fallback que nunca falha
            self.logger.debug(f"Cache key fallback usado: {e}")
            return (f"P_{id(P)}_{id(P.graph)}", f"Q_{id(Q)}_{id(Q.graph)}")

    def _rmcs_with_bbp_iterative(self, P: GraphPart, Q: GraphPart) -> Tuple[Grafo, float]:
        """RMCS completamente iterativo com tratamento robusto de erros - CORREÇÃO APLICADA"""
        self.stats['calls_rmcs'] += 1

        if self._check_timeout():
            return (Grafo(), 0.0)

        # Verificação de cache com tratamento robusto
        cache_key = self._create_cache_key(P, Q)

        try:
            cached_result = self.memo.get(cache_key)
            if cached_result is not None:
                self.stats['cache_hits'] += 1
                return cached_result
        except Exception as cache_error:
            self.logger.debug(f"Erro no cache get: {cache_error}")
            # Continuar sem cache em caso de erro

        try:
            # CASOS BASE - Processamento direto para evitar complexidade desnecessária

            # Caso 1: Rótulos incompatíveis
            try:
                label_P = P.graph.get_rotulo_vertice(P.root)
                label_Q = Q.graph.get_rotulo_vertice(Q.root)
                if label_P != label_Q:
                    result = (Grafo(), 0.0)
                    self.memo.set(cache_key, result)
                    return result
            except:
                pass

            # Caso 2: Grafos muito pequenos
            if (len(P.graph.vertices()) <= 2 or len(Q.graph.vertices()) <= 2 or
                    P.is_basic and Q.is_basic):
                return self._process_basic_direct_case(P, Q)

            # SISTEMA ITERATIVO SIMPLIFICADO
            request_stack = [(P, Q)]
            results_map = {}

            max_iterations = min(5000, len(P.graph.vertices()) * len(Q.graph.vertices()) * 5)
            iteration = 0

            while request_stack and iteration < max_iterations:
                if self._check_timeout():
                    break

                current_P, current_Q = request_stack.pop()
                current_key = self._create_cache_key(current_P, current_Q)

                # Pular se já processado
                if current_key in results_map:
                    continue

                # Verificar se é caso folha
                if self._is_leaf_case(current_P, current_Q):
                    result = self._process_leaf_case(current_P, current_Q)
                    results_map[current_key] = result
                    continue

                # Obter subtarefas
                subtasks = self._get_subtasks_safe(current_P, current_Q)

                if not subtasks:
                    # Caso sem subtarefas
                    result = self._process_leaf_case(current_P, current_Q)
                    results_map[current_key] = result
                else:
                    # Verificar se todas as subtarefas estão resolvidas
                    all_ready = True
                    pending_subtasks = []

                    for subtask_P, subtask_Q in subtasks:
                        subtask_key = self._create_cache_key(subtask_P, subtask_Q)

                        if subtask_key in results_map:
                            continue
                        elif self._can_process_immediately(subtask_P, subtask_Q):
                            # Processar imediatamente
                            subtask_result = self._process_immediate_case(subtask_P, subtask_Q)
                            results_map[subtask_key] = subtask_result
                        else:
                            all_ready = False
                            pending_subtasks.append((subtask_P, subtask_Q))

                    if all_ready:
                        # Todas as subtarefas prontas - processar
                        result = self._process_with_subtasks(current_P, current_Q, [
                            (sub_P, sub_Q, results_map[self._create_cache_key(sub_P, sub_Q)])
                            for sub_P, sub_Q in subtasks
                        ])
                        results_map[current_key] = result
                    else:
                        # Re-adicionar caso atual e adicionar subtarefas pendentes
                        request_stack.append((current_P, current_Q))
                        for sub_P, sub_Q in pending_subtasks:
                            if len(request_stack) < 1000:  # Limite de segurança
                                request_stack.append((sub_P, sub_Q))

                iteration += 1
                self.stats['iterations'] += 1

            # Retornar resultado ou fallback
            final_result = results_map.get(cache_key, (Grafo(), 0.0))

            # Atualizar cache apenas se não houve timeout
            if not self._check_timeout():
                try:
                    self.memo.set(cache_key, final_result)
                except:
                    pass  # Ignorar erros de cache

            return final_result

        except Exception as e:
            self.logger.warning(f"Erro em RMCS iterativo, usando fallback: {e}")
            # Fallback robusto
            try:
                return self._compute_approximate_mcs(P.graph, Q.graph)
            except Exception as fallback_error:
                self.logger.error(f"Erro no fallback também: {fallback_error}")
                return (Grafo(), 0.0)

    def _rmcs_with_bbp_fast(self, P: GraphPart, Q: GraphPart) -> Tuple[Grafo, float]:
        """Versão rápida do RMCS com limites rigorosos"""
        # Timeout rápido para esta chamada
        fast_timeout = 2.0  # 2 segundos máximo
        start_time = time.time()

        def fast_timeout_check():
            return time.time() - start_time > fast_timeout

        if fast_timeout_check():
            return (Grafo(), 0.0)

        # Usar abordagem simplificada para casos pequenos
        if (len(P.graph.vertices()) <= 5 and len(Q.graph.vertices()) <= 5):
            return self._process_small_case_directly(P, Q)

        # Limitar profundidade da recursão iterativa
        return self._rmcs_with_bbp_iterative_limited(P, Q, max_depth=3)

    def _get_subtasks_safe(self, P: GraphPart, Q: GraphPart) -> List[Tuple[GraphPart, GraphPart]]:
        """Obtém subtarefas com tratamento robusto de erros"""
        try:
            return self._get_subtasks(P, Q)
        except Exception as e:
            self.logger.debug(f"Erro ao obter subtarefas: {e}")
            return []

    def _is_leaf_case(self, P: GraphPart, Q: GraphPart) -> bool:
        """Verifica se é um caso folha que pode ser processado diretamente"""
        try:
            # Casos muito simples
            if (len(P.graph.vertices()) <= 1 or len(Q.graph.vertices()) <= 1 or
                    P.graph.grau(P.root) == 0 or Q.graph.grau(Q.root) == 0):
                return True

            # Tipos fundamentais incompatíveis
            if P.type != Q.type and not self._are_types_compatible(P.type, Q.type):
                return True

            return False
        except:
            return True  # Em caso de erro, tratar como folha

    def _process_basic_direct_case(self, P: GraphPart, Q: GraphPart) -> Tuple[Grafo, float]:
        """Processa casos básicos diretamente sem decomposição complexa"""
        try:
            # Verificar compatibilidade do vértice raiz
            if P.graph.get_rotulo_vertice(P.root) == Q.graph.get_rotulo_vertice(Q.root):
                base_graph, base_size = self._create_single_vertex_graph(P)

                # Tentar adicionar arestas compatíveis
                best_graph = base_graph
                best_size = base_size

                # Verificar vizinhança imediata
                p_neighbors = P.graph.vizinhanca(P.root)
                q_neighbors = Q.graph.vizinhanca(Q.root)

                for p_n in p_neighbors:
                    for q_n in q_neighbors:
                        try:
                            if (self._are_vertices_compatible(P.graph, p_n, Q.graph, q_n) and
                                    self._are_edges_compatible(P.graph, P.root, p_n, Q.graph, Q.root, q_n)):

                                # Criar grafo com aresta
                                candidate = self._create_edge_graph(P, P.root, p_n,
                                                                    P.graph.get_rotulo_aresta(P.root, p_n))
                                candidate_size = self._calculate_graph_size(candidate)

                                if candidate_size > best_size:
                                    best_graph = candidate
                                    best_size = candidate_size
                        except:
                            continue

                return (best_graph, best_size)
            else:
                return (Grafo(), 0.0)
        except Exception as e:
            self.logger.debug(f"Erro no processamento básico: {e}")
            return (Grafo(), 0.0)

    def _can_process_immediately(self, P: GraphPart, Q: GraphPart) -> bool:
        """Verifica se o caso pode ser processado imediatamente"""
        try:
            # Caso 1: rótulos incompatíveis
            label_P = P.graph.get_rotulo_vertice(P.root)
            label_Q = Q.graph.get_rotulo_vertice(Q.root)
            if label_P != label_Q:
                return True

            # Caso 2: tipos fundamentalmente incompatíveis
            if P.type != Q.type and not self._are_types_compatible(P.type, Q.type):
                return True

            # Caso 3: grafos muito pequenos
            if len(P.graph.vertices()) <= 1 or len(Q.graph.vertices()) <= 1:
                return True

            return False
        except:
            return True  # Em caso de erro, processar imediatamente

    def _are_types_compatible(self, type1: str, type2: str) -> bool:
        """Verifica se tipos de partes são compatíveis"""
        compatible_pairs = {('BPS', 'BPS'), ('BSS', 'BSS')}
        return (type1, type2) in compatible_pairs

    def _process_immediate_case(self, P: GraphPart, Q: GraphPart) -> Tuple[Grafo, float]:
        """Processa casos que podem ser resolvidos imediatamente"""
        try:
            # Se rótulos são compatíveis, retorna vértice único
            if P.graph.get_rotulo_vertice(P.root) == Q.graph.get_rotulo_vertice(Q.root):
                return self._create_single_vertex_graph(P)
        except:
            pass

        # Para tipos incompatíveis ou rótulos diferentes, retorna grafo vazio
        return (Grafo(), 0.0)

    def _get_subtasks(self, P: GraphPart, Q: GraphPart) -> List[Tuple[GraphPart, GraphPart]]:
        """Obtém subtarefas para processamento iterativo"""
        subtasks = []

        try:
            # Caso composto BPS
            if P.type == 'BPS' and Q.type == 'BPS' and not P.is_basic and not Q.is_basic:
                ep_P = self._get_elementary_parts(P)
                ep_Q = self._get_elementary_parts(Q)

                # Adicionar todos os pares compatíveis como subtarefas
                for part_P in ep_P:
                    for part_Q in ep_Q:
                        if self._are_parts_compatible(part_P, part_Q):
                            subtasks.append((part_P, part_Q))

            # Caso BPS básico
            elif P.type == 'BPS' and Q.type == 'BPS':
                p_neighbors = P.graph.vizinhanca(P.root)
                q_neighbors = Q.graph.vizinhanca(Q.root)

                for p_n in p_neighbors:
                    for q_n in q_neighbors:
                        if (self._are_vertices_compatible(P.graph, p_n, Q.graph, q_n) and
                                self._are_edges_compatible(P.graph, P.root, p_n, Q.graph, Q.root, q_n) and
                                self._is_edge_bbp_compatible(P.graph, Q.graph, P.root, p_n, Q.root, q_n)):
                            part_P_n = GraphPart(
                                graph=P.graph,
                                root=p_n,
                                is_basic=self._is_basic_root(P.graph, p_n),
                                part_type='BPS'
                            )
                            part_Q_n = GraphPart(
                                graph=Q.graph,
                                root=q_n,
                                is_basic=self._is_basic_root(Q.graph, q_n),
                                part_type='BPS'
                            )
                            subtasks.append((part_P_n, part_Q_n))

            # Caso BSS
            elif P.type == 'BSS' and Q.type == 'BSS':
                if (P.split_edge and Q.split_edge and
                        P.orientation == Q.orientation and
                        self._are_edges_compatible(P.graph, P.split_edge[0], P.split_edge[1],
                                                   Q.graph, Q.split_edge[0], Q.split_edge[1])):
                    r_prime = P.split_edge[0] if P.split_edge[1] == P.root else P.split_edge[1]
                    s_prime = Q.split_edge[0] if Q.split_edge[1] == Q.root else Q.split_edge[1]

                    part_P_expand = GraphPart(
                        graph=P.graph,
                        root=r_prime,
                        is_basic=self._is_basic_root(P.graph, r_prime),
                        part_type='BPS',
                        orientation=P.orientation
                    )
                    part_Q_expand = GraphPart(
                        graph=Q.graph,
                        root=s_prime,
                        is_basic=self._is_basic_root(Q.graph, s_prime),
                        part_type='BPS',
                        orientation=Q.orientation
                    )
                    subtasks.append((part_P_expand, part_Q_expand))

        except Exception as e:
            self.logger.warning(f"Erro ao obter subtarefas: {e}")

        return subtasks

    def _process_leaf_case(self, P: GraphPart, Q: GraphPart) -> Tuple[Grafo, float]:
        """Processa casos folha (sem subtarefas)"""
        try:
            # Caso básico: apenas o vértice raiz se compatível
            if P.graph.get_rotulo_vertice(P.root) == Q.graph.get_rotulo_vertice(Q.root):
                return self._create_single_vertex_graph(P)
        except:
            pass

        return (Grafo(), 0.0)

    def _process_with_subtasks(self, P: GraphPart, Q: GraphPart, subtasks: List[Tuple[GraphPart, GraphPart]]) -> Tuple[
        Grafo, float]:
        """Processa caso usando resultados das subtarefas"""
        try:
            # Coletar resultados das subtarefas
            subtask_results = []
            for subtask_P, subtask_Q, result in subtasks:
                if result[1] > 0:
                    subtask_results.append((subtask_P, subtask_Q, result))

            # Caso composto BPS - usar matching
            if P.type == 'BPS' and Q.type == 'BPS' and not P.is_basic and not Q.is_basic:
                return self._process_compound_case(P, Q, subtask_results)

            # Caso BPS básico - encontrar melhor combinação
            elif P.type == 'BPS' and Q.type == 'BPS':
                return self._process_basic_case(P, Q, subtask_results)

            # Caso BSS - expandir aresta
            elif P.type == 'BSS' and Q.type == 'BSS':
                return self._process_bss_case(P, Q, subtask_results)

            return (Grafo(), 0.0)
        except Exception as e:
            self.logger.warning(f"Erro ao processar com subtarefas: {e}")
            return (Grafo(), 0.0)

    def _process_compound_case(self, P: GraphPart, Q: GraphPart, subtask_results: List) -> Tuple[Grafo, float]:
        """Processa caso composto com matching bipartido"""
        if not subtask_results:
            return self._create_single_vertex_graph(P)

        try:
            # Agrupar por partes elementares
            ep_P_dict = {part: part for part in self._get_elementary_parts(P)}
            ep_Q_dict = {part: part for part in self._get_elementary_parts(Q)}

            # Criar matriz de pesos para matching
            ep_P_list = list(ep_P_dict.values())
            ep_Q_list = list(ep_Q_dict.values())

            weight_matrix = np.zeros((len(ep_P_list), len(ep_Q_list)))

            for i, part_P in enumerate(ep_P_list):
                for j, part_Q in enumerate(ep_Q_list):
                    # Encontrar resultado correspondente
                    for subtask_P, subtask_Q, result in subtask_results:
                        if (part_P.root == subtask_P.root and part_Q.root == subtask_Q.root and
                                part_P.type == subtask_P.type and part_Q.type == subtask_Q.type):
                            weight_matrix[i][j] = result[1]
                            break

            # Aplicar matching bipartido
            try:
                row_ind, col_ind = linear_sum_assignment(-weight_matrix)
                total_weight = weight_matrix[row_ind, col_ind].sum()
            except:
                total_weight = 0.0

            # Construir grafo resultante
            base_graph, base_size = self._create_single_vertex_graph(P)
            result_graph = base_graph.copy()

            for i, j in zip(row_ind, col_ind):
                if i < len(ep_P_list) and j < len(ep_Q_list) and weight_matrix[i][j] > 0:
                    part_P, part_Q = ep_P_list[i], ep_Q_list[j]
                    for subtask_P, subtask_Q, result in subtask_results:
                        if (part_P.root == subtask_P.root and part_Q.root == subtask_Q.root):
                            result_graph = self._merge_graphs(result_graph, result[0])
                            break

            # Ajustar tamanho
            root_weight = self.label_weights.get(P.graph.get_rotulo_vertice(P.root), 1.0)
            final_size = self._calculate_graph_size(result_graph)
            if len(row_ind) > 0:
                final_size -= (len(row_ind) - 1) * root_weight

            return (result_graph, max(base_size, final_size))
        except Exception as e:
            self.logger.warning(f"Erro no processamento de caso composto: {e}")
            return self._create_single_vertex_graph(P)

    def _process_basic_case(self, P: GraphPart, Q: GraphPart, subtask_results: List) -> Tuple[Grafo, float]:
        """Processa caso BPS básico"""
        base_graph, best_size = self._create_single_vertex_graph(P)

        if not subtask_results:
            return (base_graph, best_size)

        # Encontrar melhor combinação de arestas
        best_candidate = base_graph
        best_candidate_size = best_size

        for subtask_P, subtask_Q, result in subtask_results:
            if result[1] > 0:
                # Verificar compatibilidade da aresta
                try:
                    edge_label_P = P.graph.get_rotulo_aresta(P.root, subtask_P.root)
                    edge_label_Q = Q.graph.get_rotulo_aresta(Q.root, subtask_Q.root)

                    if edge_label_P == edge_label_Q:
                        edge_graph = self._create_edge_graph(P, P.root, subtask_P.root, edge_label_P)
                        candidate = self._merge_graphs(edge_graph, result[0])
                        candidate_size = self._calculate_graph_size(candidate)

                        if candidate_size > best_candidate_size:
                            best_candidate = candidate
                            best_candidate_size = candidate_size
                except:
                    continue

        return (best_candidate, best_candidate_size)

    def _process_bss_case(self, P: GraphPart, Q: GraphPart, subtask_results: List) -> Tuple[Grafo, float]:
        """Processa caso BSS"""
        if not P.split_edge or not Q.split_edge:
            return self._create_single_vertex_graph(P)

        # Criar grafo base com aresta split
        try:
            edge_label_P = P.graph.get_rotulo_aresta(P.split_edge[0], P.split_edge[1])
            base_graph = self._create_edge_graph(P, P.split_edge[0], P.split_edge[1], edge_label_P)
            base_size = self._calculate_graph_size(base_graph)
        except:
            return self._create_single_vertex_graph(P)

        if not subtask_results:
            return (base_graph, base_size)

        # Expandir com subtarefas
        best_candidate = base_graph
        best_candidate_size = base_size

        for subtask_P, subtask_Q, result in subtask_results:
            if result[1] > 0:
                expanded_graph = self._merge_graphs(base_graph, result[0])
                expanded_size = self._calculate_graph_size(expanded_graph)

                if expanded_size > best_candidate_size:
                    best_candidate = expanded_graph
                    best_candidate_size = expanded_size

        return (best_candidate, best_candidate_size)

    def _are_edges_compatible(self, G: Grafo, u_g: Any, v_g: Any, H: Grafo, u_h: Any, v_h: Any) -> bool:
        """Verifica compatibilidade de arestas"""
        try:
            label_G = G.get_rotulo_aresta(u_g, v_g)
            label_H = H.get_rotulo_aresta(u_h, v_h)
            return label_G == label_H
        except:
            return False

    def _get_elementary_parts(self, part: GraphPart) -> List[GraphPart]:
        """Obtém partes elementares (versão iterativa)"""
        self.stats['decompositions'] += 1
        parts = []
        graph = part.graph
        root = part.root

        try:
            if part.type == 'BPS' and not part.is_basic:
                neighbors = graph.vizinhanca(root)
                for neighbor in neighbors:
                    new_part = GraphPart(
                        graph=graph,
                        root=neighbor,
                        is_basic=self._is_basic_root(graph, neighbor),
                        part_type='BPS',
                        parent_root=root
                    )
                    parts.append(new_part)

            elif part.type == 'BPS' and part.is_basic:
                if graph.grau(root) == 1:
                    neighbors = graph.vizinhanca(root)
                    if neighbors:
                        neighbor = neighbors[0]
                        new_part = GraphPart(
                            graph=graph,
                            root=neighbor,
                            is_basic=False,
                            part_type='BPS',
                            parent_root=root
                        )
                        parts.append(new_part)
                else:
                    neighbors = graph.vizinhanca(root)
                    for neighbor in neighbors:
                        bss_clock = GraphPart(
                            graph=graph,
                            root=neighbor,
                            is_basic=True,
                            part_type='BSS',
                            split_edge=(root, neighbor),
                            orientation='clockwise',
                            parent_root=root
                        )
                        parts.append(bss_clock)

                        bss_counter = GraphPart(
                            graph=graph,
                            root=neighbor,
                            is_basic=True,
                            part_type='BSS',
                            split_edge=(root, neighbor),
                            orientation='counterclockwise',
                            parent_root=root
                        )
                        parts.append(bss_counter)

            elif part.type == 'BSS':
                orientation = part.orientation if part.orientation else 'clockwise'
                parts.extend(self._decompose_bss_with_orientation(part, orientation))

        except Exception as e:
            self.logger.warning(f"Erro na decomposição de partes: {e}")

        return parts

    def _decompose_bss_with_orientation(self, bss: GraphPart, orientation: str) -> List[GraphPart]:
        """Decomposição de BSS considerando orientação"""
        parts = []
        graph = bss.graph
        root = bss.root
        parent_root = bss.parent_root

        try:
            neighbors = graph.vizinhanca(root)
            for neighbor in neighbors:
                if neighbor != parent_root:
                    new_part = GraphPart(
                        graph=graph,
                        root=neighbor,
                        is_basic=self._is_basic_root(graph, neighbor),
                        part_type='BSS',
                        split_edge=(root, neighbor),
                        orientation=orientation,
                        parent_root=root
                    )
                    parts.append(new_part)
        except Exception as e:
            self.logger.warning(f"Erro na decomposição BSS: {e}")

        return parts

    def _decompose_graph(self, graph: Grafo, root: Any) -> List[GraphPart]:
        """Decomposição completa do grafo (iterativa)"""
        cache_key = (self._graph_signature(graph), root)
        cached_result = self.decomposition_cache.get(cache_key)
        if cached_result is not None:
            self.stats['cache_hits'] += 1
            return cached_result

        is_basic = self._is_basic_root(graph, root)

        root_part = GraphPart(
            graph=graph,
            root=root,
            is_basic=is_basic,
            part_type='BPS'
        )

        parts = [root_part]
        queue = deque([root_part])
        visited_parts = set([str(root_part)])

        while queue and len(parts) < 50:
            current = queue.popleft()
            elementary_parts = self._get_elementary_parts(current)

            for part in elementary_parts:
                part_str = str(part)
                if part_str not in visited_parts:
                    parts.append(part)
                    visited_parts.add(part_str)
                    queue.append(part)

        self.decomposition_cache.set(cache_key, parts)
        return parts

    def _find_root_candidates(self, graph: Grafo) -> List:
        """Seleção otimizada de root-candidates"""
        try:
            candidates = []
            vertices = graph.vertices()

            # Para grafos pequenos, usar todos os vértices
            if len(vertices) <= 8:
                return vertices

            # Para grafos médios, usar vértices com maior grau e rótulos comuns
            scored_vertices = []
            for v in vertices:
                try:
                    degree = graph.grau(v)
                    label = graph.get_rotulo_vertice(v)
                    score = degree

                    # Priorizar rótulos comuns em química
                    if label in ['C', 'N', 'O', 'S', 'P']:
                        score += 5
                    elif label in ['F', 'Cl', 'Br', 'I']:
                        score += 3

                    scored_vertices.append((v, score))
                except:
                    continue

            scored_vertices.sort(key=lambda x: x[1], reverse=True)
            return [v for v, score in scored_vertices[:12]]  # Limitar a 12 candidatos

        except Exception as e:
            self.logger.debug(f"Erro na seleção de candidatos: {e}")
            return graph.vertices()[:6]

    def _find_root_candidates_fallback(self, graph: Grafo) -> List:
        """Fallback para seleção de root-candidates"""
        try:
            vertices = graph.vertices()

            if len(vertices) <= 6:
                return vertices

            scored_vertices = []
            for v in vertices:
                try:
                    degree = graph.grau(v)
                    label = graph.get_rotulo_vertice(v)
                    score = degree * 2
                    if label in ['C', 'N', 'O']:
                        score += 3
                    scored_vertices.append((v, score))
                except:
                    continue

            scored_vertices.sort(key=lambda x: x[1], reverse=True)
            return [v for v, score in scored_vertices[:8]]
        except:
            return graph.vertices()[:4]

    def _quick_similarity_check(self, parts_G: List[GraphPart], parts_H: List[GraphPart]) -> bool:
        """Verificação rápida de similaridade"""
        if not parts_G or not parts_H:
            return False

        labels_G = set()
        for part in parts_G:
            try:
                label = part.graph.get_rotulo_vertice(part.root)
                if label:
                    labels_G.add(label)
            except:
                continue

        labels_H = set()
        for part in parts_H:
            try:
                label = part.graph.get_rotulo_vertice(part.root)
                if label:
                    labels_H.add(label)
            except:
                continue

        common_labels = labels_G & labels_H
        return len(common_labels) > 0

    def _is_basic_root(self, graph: Grafo, root: Any) -> bool:
        """Verificação se raiz é básica"""
        try:
            degree = graph.grau(root)
            if degree <= 1:
                return True

            # Verificar se está em apenas um bloco
            if hasattr(graph, 'encontrar_blocos') and callable(graph.encontrar_blocos):
                blocks = graph.encontrar_blocos()
                incident_blocks = [block for block in blocks if root in block]
                return len(incident_blocks) == 1
            else:
                # Fallback: considerar como básico se grau <= 2
                return degree <= 2
        except:
            return True

    def _is_bridge(self, graph: Grafo, u: Any, v: Any) -> bool:
        """Verifica se uma aresta é uma bridge"""
        try:
            if hasattr(graph, 'encontrar_pontes') and callable(graph.encontrar_pontes):
                bridges = graph.encontrar_pontes()
                return (u, v) in bridges or (v, u) in bridges
            else:
                # Fallback: criar cópia e verificar conectividade
                temp_graph = graph.copy()
                # Assumindo que existe método para remover aresta
                if hasattr(temp_graph, 'remover_aresta'):
                    temp_graph.remover_aresta(u, v)
                    return not self._is_connected(temp_graph)
                else:
                    return False
        except:
            return False

    def _is_connected(self, graph: Grafo) -> bool:
        """Verifica se o grafo é conectado"""
        try:
            vertices = graph.vertices()
            if not vertices:
                return True
            visited = set()
            queue = deque([vertices[0]])
            while queue:
                current = queue.popleft()
                if current not in visited:
                    visited.add(current)
                    neighbors = graph.vizinhanca(current)
                    queue.extend([n for n in neighbors if n not in visited])
            return len(visited) == len(vertices)
        except:
            return False

    def _are_vertices_compatible(self, G: Grafo, v_g: Any, H: Grafo, v_h: Any) -> bool:
        """Verifica compatibilidade de vértices"""
        try:
            return G.get_rotulo_vertice(v_g) == H.get_rotulo_vertice(v_h)
        except:
            return False

    def _is_edge_bbp_compatible(self, G: Grafo, H: Grafo, u_g: Any, v_g: Any, u_h: Any, v_h: Any) -> bool:
        """Verifica compatibilidade BBP para arestas"""
        try:
            is_bridge_G = self._is_bridge(G, u_g, v_g)
            is_bridge_H = self._is_bridge(H, u_h, v_h)
            return is_bridge_G == is_bridge_H
        except:
            return True

    def _are_parts_compatible(self, P: GraphPart, Q: GraphPart) -> bool:
        """Verificação de compatibilidade entre partes"""
        try:
            label_P = P.graph.get_rotulo_vertice(P.root)
            label_Q = Q.graph.get_rotulo_vertice(Q.root)
            return label_P == label_Q
        except:
            return False

    def _are_graphs_identical(self, G: Grafo, H: Grafo) -> bool:
        """Verificação robusta de grafos idênticos"""
        try:
            if len(G.vertices()) != len(H.vertices()):
                return False
            if len(G.arestas()) != len(H.arestas()):
                return False

            # Verificar vértices
            g_vertices = sorted(G.vertices())
            h_vertices = sorted(H.vertices())

            for v_g, v_h in zip(g_vertices, h_vertices):
                if G.get_rotulo_vertice(v_g) != H.get_rotulo_vertice(v_h):
                    return False

            # Verificar arestas
            g_edges = sorted(G.arestas())
            h_edges = sorted(H.arestas())

            for (u_g, v_g), (u_h, v_h) in zip(g_edges, h_edges):
                if G.get_rotulo_aresta(u_g, v_g) != H.get_rotulo_aresta(u_h, v_h):
                    return False

            return True
        except:
            return False

    def _create_single_vertex_graph(self, part: GraphPart) -> Tuple[Grafo, float]:
        """Cria grafo com vértice único"""
        try:
            graph = Grafo()
            root = part.root
            label = part.graph.get_rotulo_vertice(root)
            graph.adicionar_vertice(root, label)
            size = self._calculate_graph_size(graph)
            return (graph, size)
        except Exception as e:
            self.logger.warning(f"Erro ao criar vértice único: {e}")
            return (Grafo(), 0.0)

    def _create_edge_graph(self, part: GraphPart, u: Any, v: Any, edge_label: str) -> Grafo:
        """Cria grafo com aresta"""
        try:
            graph = Grafo()
            label_u = part.graph.get_rotulo_vertice(u)
            label_v = part.graph.get_rotulo_vertice(v)
            graph.adicionar_vertice(u, label_u)
            graph.adicionar_vertice(v, label_v)
            graph.adicionar_aresta(u, v, edge_label)
            return graph
        except Exception as e:
            self.logger.warning(f"Erro ao criar grafo com aresta: {e}")
            return Grafo()

    def _merge_graphs(self, graph1: Grafo, graph2: Grafo) -> Grafo:
        """Fusão de dois grafos"""
        try:
            result = graph1.copy()
            for v in graph2.vertices():
                if not result.existe_vertice(v):
                    label = graph2.get_rotulo_vertice(v)
                    result.adicionar_vertice(v, label)
            for u, v in graph2.arestas():
                if not result.existe_aresta(u, v):
                    label = graph2.get_rotulo_aresta(u, v)
                    result.adicionar_aresta(u, v, label)
            return result
        except Exception as e:
            self.logger.warning(f"Erro ao mesclar grafos: {e}")
            return graph1

    def _calculate_graph_size(self, graph: Grafo) -> float:
        """Calcula tamanho do grafo usando pesos dos rótulos"""
        try:
            total_size = 0.0
            for v in graph.vertices():
                label = graph.get_rotulo_vertice(v)
                if label is None:
                    label = "C"
                weight = self.label_weights.get(label, 1.0)
                total_size += weight
            for u, v in graph.arestas():
                label = graph.get_rotulo_aresta(u, v)
                if label is None:
                    label = "single"
                weight = self.label_weights.get(label, 1.0)
                total_size += weight
            return total_size
        except:
            return 0.0

    def _compute_enhanced_fallback(self, G: Grafo, H: Grafo) -> Tuple[Grafo, float]:
        """Fallback aprimorado para casos difíceis"""
        try:
            # Tentar abordagem gulosa primeiro
            greedy_result = self._compute_greedy_mcs(G, H)
            if greedy_result[1] > 0:
                return greedy_result

            # Fallback para vértices únicos
            return self._compute_vertex_only_mcs(G, H)
        except Exception as e:
            self.logger.warning(f"Erro no fallback aprimorado: {e}")
            return (Grafo(), 0.0)

    def _compute_vertex_only_mcs(self, G: Grafo, H: Grafo) -> Tuple[Grafo, float]:
        """Fallback: MCS apenas com vértices (sem arestas)"""
        try:
            mcs_graph = Grafo()
            total_size = 0.0

            # Adicionar todos os vértices com labels comuns
            common_labels = set()
            for v in G.vertices():
                label_g = G.get_rotulo_vertice(v)
                for u in H.vertices():
                    if H.get_rotulo_vertice(u) == label_g:
                        common_labels.add(label_g)
                        break

            for label in common_labels:
                vertices_g = [v for v in G.vertices() if G.get_rotulo_vertice(v) == label]
                vertices_h = [v for v in H.vertices() if H.get_rotulo_vertice(v) == label]

                # Adicionar até min(len(vertices_g), len(vertices_h)) vértices
                count = min(len(vertices_g), len(vertices_h))
                for i in range(count):
                    v = vertices_g[i]
                    mcs_graph.adicionar_vertice(v, label)
                    total_size += self.label_weights.get(label, 1.0)

            return (mcs_graph, total_size)
        except Exception as e:
            self.logger.debug(f"Erro no MCS de vértices: {e}")
            return (Grafo(), 0.0)

    def _compute_greedy_mcs(self, G: Grafo, H: Grafo) -> Tuple[Grafo, float]:
        """Abordagem gulosa para MCS"""
        try:
            best_graph = Grafo()
            best_size = 0.0

            # Encontrar todos os vértices compatíveis
            compatible_vertices = []
            for v_g in G.vertices():
                for v_h in H.vertices():
                    if self._are_vertices_compatible(G, v_g, H, v_h):
                        compatible_vertices.append((v_g, v_h))

            # Ordenar por grau (mais conectados primeiro)
            compatible_vertices.sort(key=lambda pair: G.grau(pair[0]) + H.grau(pair[1]), reverse=True)

            # Tentar construir a partir dos melhores candidatos
            for v_g, v_h in compatible_vertices[:10]:  # Limitar a 10 tentativas
                if self._check_timeout():
                    break

                candidate = self._grow_mcs_from_vertex(G, H, v_g, v_h)
                candidate_size = self._calculate_graph_size(candidate)

                if candidate_size > best_size:
                    best_graph = candidate
                    best_size = candidate_size

            return (best_graph, best_size)
        except Exception as e:
            self.logger.debug(f"Erro no MCS guloso: {e}")
            return (Grafo(), 0.0)

    def _grow_mcs_from_vertex(self, G: Grafo, H: Grafo, start_g: Any, start_h: Any) -> Grafo:
        """Cresce MCS a partir de um vértice inicial"""
        try:
            mcs_graph = Grafo()
            mcs_graph.adicionar_vertice(start_g, G.get_rotulo_vertice(start_g))
            mapping = {start_g: start_h}
            queue = deque([start_g])

            while queue and not self._check_timeout():
                current_g = queue.popleft()
                current_h = mapping[current_g]

                # Explorar vizinhança
                for neighbor_g in G.vizinhanca(current_g):
                    if neighbor_g in mapping:
                        continue

                    # Encontrar vizinho compatível em H
                    found = False
                    for neighbor_h in H.vizinhanca(current_h):
                        if neighbor_h in mapping.values():
                            continue

                        if (self._are_vertices_compatible(G, neighbor_g, H, neighbor_h) and
                                self._are_edges_compatible(G, current_g, neighbor_g, H, current_h, neighbor_h)):
                            # Adicionar ao MCS
                            mcs_graph.adicionar_vertice(neighbor_g, G.get_rotulo_vertice(neighbor_g))
                            mcs_graph.adicionar_aresta(current_g, neighbor_g,
                                                       G.get_rotulo_aresta(current_g, neighbor_g))
                            mapping[neighbor_g] = neighbor_h
                            queue.append(neighbor_g)
                            found = True
                            break

                # Limitar crescimento para evitar explosão combinatorial
                if len(mcs_graph.vertices()) > 20:
                    break

            return mcs_graph
        except Exception as e:
            self.logger.debug(f"Erro no crescimento do MCS: {e}")
            return Grafo()

    def _compute_approximate_mcs(self, G: Grafo, H: Grafo) -> Tuple[Grafo, float]:
        """Computação aproximada para fallback"""
        try:
            best_graph = Grafo()
            best_size = 0.0

            # Encontrar vértices com rótulos comuns
            common_labels = set()
            for v in G.vertices():
                label_g = G.get_rotulo_vertice(v)
                for u in H.vertices():
                    if H.get_rotulo_vertice(u) == label_g:
                        common_labels.add(label_g)
                        break

            # Tentar construir MCS a partir de cada par compatível
            for label in common_labels:
                if self._check_timeout():
                    break

                g_vertices = [v for v in G.vertices() if G.get_rotulo_vertice(v) == label]
                h_vertices = [v for v in H.vertices() if H.get_rotulo_vertice(v) == label]

                for g_start in g_vertices[:3]:
                    for h_start in h_vertices[:3]:
                        candidate, size = self._build_mcs_from_pair(G, H, g_start, h_start)
                        if size > best_size:
                            best_graph = candidate
                            best_size = size

            return (best_graph, best_size)
        except Exception as e:
            self.logger.warning(f"Erro no fallback aproximado: {e}")
            return (Grafo(), 0.0)

    def _build_mcs_from_pair(self, G: Grafo, H: Grafo, g_start: Any, h_start: Any) -> Tuple[Grafo, float]:
        """Construir MCS a partir de um par de vértices compatíveis"""
        try:
            mcs_graph = Grafo()
            mcs_graph.adicionar_vertice(g_start, G.get_rotulo_vertice(g_start))
            mapping = {g_start: h_start}
            queue = deque([g_start])

            while queue and not self._check_timeout():
                g_v = queue.popleft()
                h_v = mapping[g_v]

                for g_neighbor in G.vizinhanca(g_v):
                    if g_neighbor in mapping:
                        continue

                    g_label = G.get_rotulo_vertice(g_neighbor)
                    edge_label = G.get_rotulo_aresta(g_v, g_neighbor)

                    for h_neighbor in H.vizinhanca(h_v):
                        if h_neighbor in mapping.values():
                            continue

                        if (H.get_rotulo_vertice(h_neighbor) == g_label and
                                H.get_rotulo_aresta(h_v, h_neighbor) == edge_label):
                            mcs_graph.adicionar_vertice(g_neighbor, g_label)
                            mcs_graph.adicionar_aresta(g_v, g_neighbor, edge_label)
                            mapping[g_neighbor] = h_neighbor
                            queue.append(g_neighbor)
                            break

            return (mcs_graph, self._calculate_graph_size(mcs_graph))
        except Exception as e:
            self.logger.warning(f"Erro ao construir MCS a partir de par: {e}")
            return (Grafo(), 0.0)

    def _check_timeout(self) -> bool:
        """Verifica timeout com tratamento robusto"""
        try:
            if time.time() - self.start_time > self.timeout:
                self._timeout_occurred = True
                self.stats['timeout_occurred'] = True
                self.logger.warning("Timeout ocorreu durante o processamento")
                return True
            return False
        except:
            return False

    def get_statistics(self) -> Dict:
        """Retorna estatísticas de performance"""
        stats = self.stats.copy()
        if stats.get('start_time'):
            stats['total_time'] = time.time() - stats['start_time']
        return stats


# ============================ FUNÇÃO PRINCIPAL ============================

def calcular_mcs_outerplanar(grafo_G: Grafo, grafo_H: Grafo,
                             label_weights: Optional[Dict[str, float]] = None,
                             timeout: int = TIMEOUT_DURATION,
                             return_stats: bool = False,
                             verbose: bool = False) -> Union[Tuple[Grafo, float], Tuple[Grafo, float, Dict]]:
    """
    Função principal para cálculo de MCS entre grafos outerplanares.
    """
    try:
        solver = OuterplanarMCS(
            label_weights=label_weights,
            timeout=timeout,
            verbose=verbose
        )
        result = solver.compute_mcs(grafo_G, grafo_H)

        if return_stats:
            stats = solver.get_statistics()
            return (*result, stats)
        return result

    except Exception as e:
        logging.error(f"Erro no cálculo do MCS: {e}")
        if return_stats:
            return (Grafo(), 0.0, {'error': str(e)})
        return (Grafo(), 0.0)