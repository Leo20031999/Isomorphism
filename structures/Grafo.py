class Grafo:
    def __init__(self):
        import networkx as nx
        self._nx = nx
        self.grafo = self._nx.Graph()
        self._rotulos_vertices = {}
        self._rotulos_arestas = {}
        self._atributos_vertices = {}
        self._atributos_arestas = {}

    def adicionar_vertice(self, v, rotulo=None, atributos=None):
        self.grafo.add_node(v)
        if rotulo is not None:
            self._rotulos_vertices[v] = rotulo
        if atributos is not None:
            self._atributos_vertices[v] = atributos
        else:
            self._atributos_vertices[v] = {}

    def adicionar_aresta(self, u, v, rotulo=None, atributos=None):
        self.grafo.add_edge(u, v)
        if rotulo is not None:
            self._rotulos_arestas[(u, v)] = rotulo
            self._rotulos_arestas[(v, u)] = rotulo
        if atributos is not None:
            self._atributos_arestas[(u, v)] = atributos
            self._atributos_arestas[(v, u)] = atributos
        else:
            self._atributos_arestas[(u, v)] = {}
            self._atributos_arestas[(v, u)] = {}

    def adicionar_multiplos_vertices(self, list_v):
        for v in list_v:
            self.adicionar_vertice(v)

    def adicionar_multiplas_arestas(self, list_e):
        for u, v in list_e:
            self.adicionar_aresta(u, v)

    def vertices(self):
        return list(self.grafo.nodes())

    def arestas(self):
        return list(self.grafo.edges())

    def get_atributos_vertice(self, v):
        return self._atributos_vertices.get(v, {})

    def get_atributos_aresta(self, u, v):
        return self._atributos_arestas.get((u, v), {})

    def grau(self, v):
        return self.grafo.degree(v)

    def get_rotulo_vertice(self, v):
        return self._rotulos_vertices.get(v, None)

    def get_rotulo_aresta(self, u, v):
        return self._rotulos_arestas.get((u, v), None)

    def set_rotulo_vertice(self, v, rotulo):
        self._rotulos_vertices[v] = rotulo

    def set_rotulo_aresta(self, u, v, rotulo):
        self._rotulos_arestas[(u, v)] = rotulo
        self._rotulos_arestas[(v, u)] = rotulo

    def center(self):
        return self._nx.center(self.grafo)

    def bfs(self, root):
        visited = []
        queue = [root]
        while queue:
            v = queue.pop(0)
            if v not in visited:
                visited.append(v)
                queue.extend([n for n in self.grafo.neighbors(v) if n not in visited])
        return visited

    def dfs(self, root):
        visited = []
        stack = [root]
        while stack:
            v = stack.pop()
            if v not in visited:
                visited.append(v)
                stack.extend([n for n in self.grafo.neighbors(v) if n not in visited])
        return visited

    def encontrar_bridges(self):
        return list(self._nx.bridges(self.grafo))

    def encontrar_blocos(self):
        try:
            return list(self._nx.biconnected_components(self._nx.Graph(self.grafo)))
        except Exception:
            return []

    def encontrar_pontes(self):
        try:
            return list(self._nx.bridges(self._nx.Graph(self.grafo)))
        except Exception:
            return []

    def eh_outerplanar(self):
        try:
            return self._nx.is_outerplanar(self.grafo)
        except Exception:
            return False

    def vizinhanca(self, v):
        return list(self.grafo.neighbors(v))

    def copy(self):
        """Cria uma cópia profunda do grafo"""
        novo_grafo = Grafo()
        for v in self.vertices():
            novo_grafo.adicionar_vertice(v, self.get_rotulo_vertice(v))
        for u, v in self.arestas():
            novo_grafo.adicionar_aresta(u, v, self.get_rotulo_aresta(u, v))
        return novo_grafo

    def existe_vertice(self, v):
        """Verifica se um vértice existe no grafo"""
        return v in self.grafo

    def is_leaf(self, v, parent_map):
        """Verifica se v é uma folha na árvore enraizada, dado o parent_map."""
        return len(self.obter_filhos(v, parent_map)) == 0

    def obter_filhos(self, v, parent_map):
        """Obtém os filhos de um vértice na árvore enraizada, dado o parent_map."""
        filhos = []
        for w in self.vizinhanca(v):
            if parent_map.get(w) == v:
                filhos.append(w)
        return filhos

    def altura(self, v, parent_map):
        """Calcula a altura da subárvore enraizada em v (versão iterativa)."""
        alturas = {}
        stack = [v]
        visitados = set()

        while stack:
            node = stack[-1]
            if node in visitados:
                stack.pop()
                filhos = self.obter_filhos(node, parent_map)
                alturas_filhos = [alturas[f] for f in filhos if f in alturas]
                alturas[node] = 1 + max(alturas_filhos) if alturas_filhos else 0
            else:
                visitados.add(node)
                filhos = self.obter_filhos(node, parent_map)
                for filho in filhos:
                    if filho not in visitados:
                        stack.append(filho)

        return alturas.get(v, 0)

    def preorder(self, root, parent_map):
        """Retorna a lista de vértices em pré-ordem na árvore enraizada."""
        ordem = []
        stack = [root]
        while stack:
            node = stack.pop()
            ordem.append(node)
            filhos = self.obter_filhos(node, parent_map)
            for filho in reversed(filhos):
                stack.append(filho)
        return ordem