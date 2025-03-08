from collections import deque
import networkx as nx

class Grafo:
    def __init__(self):
        self.grafo = nx.Graph()

    def vertices(self):
        return list(self.grafo.nodes)

    def arestas(self):
        return [(min(u, v), max(u, v)) for u, v in self.grafo.edges]

    def adicionar_aresta(self, u, v, peso=None):
        if not self.grafo.has_node(u):
            self.grafo.add_node(u)
        if not self.grafo.has_node(v):
            self.grafo.add_node(v)
        if not self.grafo.has_edge(u, v):
            self.grafo.add_edge(u, v, weight=peso)

    def remover_aresta(self, u, v):
        if not self.grafo.has_edge(u, v):
            raise KeyError(f"Aresta {u}-{v} não existe.")
        self.grafo.remove_edge(u, v)

    def get_peso(self, u, v):
        if self.grafo.has_edge(u, v):
            return self.grafo[u][v].get('weight', None)
        raise KeyError(f"Aresta {u}-{v} não existe.")

    def set_peso(self, u, v, peso):
        if self.grafo.has_edge(u, v):
            self.grafo[u][v]['weight'] = peso
        else:
            raise KeyError(f"Aresta {u}-{v} não existe.")

    def vizinhanca(self, u, parent_map=None):
        """Retorna os filhos de u, excluindo o pai (se houver)."""
        if parent_map is None:
            return list(self.grafo.neighbors(u))
        parent = parent_map.get(u)
        return [v for v in self.grafo.neighbors(u) if v != parent]

    def grau(self, u):
        self._verificar_vertice(u)
        return self.grafo.degree(u)

    def is_leaf(self, v, parent_map=None):
        """Verifica se v é folha, considerando o parent_map para árvores enraizadas."""
        if parent_map is None:
            # Folha em árvore não enraizada: grau <= 1
            return self.grafo.degree(v) <= 1
        # Em árvore enraizada, folha não tem filhos além do pai
        return len(self.vizinhanca(v, parent_map)) == 0

    def altura(self, node, parent_map):
        """Calcula a altura como número de arestas no caminho mais longo."""
        if self.is_leaf(node, parent_map):
            return 0
        max_depth = 0
        stack = [(node, 0)]  # (nó, profundidade em arestas)
        while stack:
            current, depth = stack.pop()
            max_depth = max(max_depth, depth)
            for child in self.vizinhanca(current, parent_map):
                stack.append((child, depth + 1))
        return max_depth

    def sao_adj(self, u, v):
        self._verificar_vertice(u, v)
        return self.grafo.has_edge(u, v)

    def n(self, v, tipo="*"):
        self._verificar_vertice(v)
        vizinhos = self.vizinhanca(v)
        if tipo == "*":
            return vizinhos
        return [w for w in vizinhos if self.grafo[v][w].get("tipo", None) == tipo]
    
    def preorder(self, root):
        visited = set()
        order = []

        def dfs(v):
            if v in visited:
                return
            visited.add(v)
            order.append(v)
            for neighbor in self.n(v):
                if neighbor not in visited:
                    dfs(neighbor)

        dfs(root)
        return order

    def center(self):
        """
        Retorna o vértice central da árvore.
        Se houver dois centros, os centros.
        """
        return nx.center(self.grafo)