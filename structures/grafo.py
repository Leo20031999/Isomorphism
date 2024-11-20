import networkx as nx

class Grafo(object):
    def __init__(self):
        self.G = nx.Graph()

    def DefinirN(self, n):
        self.n = n
        self.G.add_nodes_from(range(1, n+1))

    def V(self):
        for v in self.G.nodes:
            yield v

    def E(self, IterarSobreNo=False):
        arestas = set()
        for u, v in self.G.edges:
            arestas.add((min(u, v), max(u, v)))
        return list(arestas)

    def AdicionarAresta(self, u, v, peso=None):
        self.G.add_edge(u, v, weight=peso)

    def RemoverAresta(self, u, v):
        self.G.remove_edge(u, v)

    def Peso(self, u, v):
        if self.G.has_edge(u, v):
            return self.G[u][v].get('weight', None)
        return None

    def Vizinhanca(self, u):
        return list(self.G.neighbors(u))

    def Grau(self, u):
        return self.G.degree(u)


