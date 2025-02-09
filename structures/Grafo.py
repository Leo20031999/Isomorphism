import networkx as nx

class Grafo:
    def __init__(self):
        self.grafo = nx.Graph()
        self.num_nos = 0

    def definir_n(self, n):
        if len(self.grafo.nodes) < n:
            self.grafo.add_nodes_from(range(len(self.grafo.nodes) + 1, n + 1))
        self.num_nos = n

    def vertices(self):
        return list(self.grafo.nodes)

    def arestas(self):
        return [(min(u, v), max(u, v)) for u, v in self.grafo.edges]

    def adicionar_aresta(self, u, v, peso=None):
        self._verificar_vertice(u, v)
        
        if not self.grafo.has_edge(u, v) and not self.grafo.has_edge(v, u):
            self.grafo.add_edge(u, v, weight=peso)
        else:
            print(f"Aresta entre {u}-{v} já existe ou é equivalente à aresta {v}-{u}.")


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

    def vizinhanca(self, u):
        self._verificar_vertice(u)
        return [v for v in self.grafo.neighbors(u) if v != self.num_nos + 1 and v != self.num_nos + 2]

    def grau(self, u):
        self._verificar_vertice(u)
        return self.grafo.degree(u)

    def is_leaf(self, v):
        self._verificar_vertice(v)
        return len(list(self.grafo.neighbors(v))) == 0

    def sao_adj(self, u, v):
        self._verificar_vertice(u, v)
        return self.grafo.has_edge(u, v)

    def altura(self, v=None):
        def dfs_iterativo(v, heights):
            visited = set()
            stack = [(v, 0)]
            while stack:
                node, depth = stack.pop()
                if node not in visited:
                    visited.add(node)
                    heights[node] = depth
                    for neighbor in self.vizinhanca(node):
                        if neighbor not in visited:
                            stack.append((neighbor, depth + 1))

        heights = {}
        if v is not None:
            dfs_iterativo(v, heights)
        else:
            for node in self.vertices():
                if node not in heights:
                    dfs_iterativo(node, heights)
        return heights if v is None else heights[v]

    def n(self, v, tipo="*"):
        self._verificar_vertice(v)
        vizinhos = self.vizinhanca(v)
        if tipo == "*":
            return vizinhos
        return [w for w in vizinhos if self.grafo[v][w].get("tipo", None) == tipo]

    def _verificar_vertice(self, *vertices):
        for v in vertices:
            if v not in self.grafo.nodes:
                raise ValueError(f"Vértice {v} não existe no grafo.")
    
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

