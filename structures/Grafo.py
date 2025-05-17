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
        if parent_map is None:
            return list(self.grafo.neighbors(u))
        parent = parent_map.get(u)
        return [v for v in self.grafo.neighbors(u) if v != parent]

    def _verificar_vertice(self, *vertices):
        """Verifica se um ou mais vértices existem no grafo (método interno)"""
        for v in vertices:
            if not self.grafo.has_node(v):
                raise ValueError(f"Vértice {v} não existe no grafo")
            
    def grau(self, u):
        self._verificar_vertice(u)
        return self.grafo.degree(u)

    def is_leaf(self, v, parent_map=None):
        if parent_map is None:
            return self.grafo.degree(v) <= 1
        return len(self.vizinhanca(v, parent_map)) == 0

    def altura(self, node, parent_map):
        if self.is_leaf(node, parent_map):
            return 0
        max_depth = 0
        stack = [(node, 0)]
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
    
    def preorder(self, root, parent_map):
        visited = []
        
        def dfs(node):
            if node is None:
                return
            visited.append(node)
            children = sorted(
                [child for child in self.vizinhanca(node, parent_map) if parent_map.get(child, None) == node],
                key=lambda x: (-self.altura(x, parent_map), x)
            )
            for child in children:
                dfs(child)
        
        dfs(root)
        return visited

    def center(self):
        return nx.center(self.grafo)
    
    def levelize(self, root):
        parent_map = {}
        layers = []
        queue = deque([root])
        parent_map[root] = None
        current_layer = [root]
        while current_layer:
            layers.append(current_layer)
            next_layer = []
            for u in current_layer:
                for v in self.grafo.neighbors(u):
                    if v not in parent_map:
                        parent_map[v] = u
                        next_layer.append(v)
            current_layer = next_layer
        return layers, parent_map