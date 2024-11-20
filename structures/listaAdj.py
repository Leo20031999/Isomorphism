import networkx as nx
from structures.grafo import Grafo

class ListaAdjG(Grafo):
    def DefinirN(self, n):
        super().DefinirN(n)

    def AdicionarAresta(self, u, v, peso=0):
        if u < 1 or u > self.n or v < 1 or v > self.n:
            print(f"Erro ao tentar adicionar aresta {u}-{v}: vértices fora do intervalo [1, {self.n}]")
            return None
        
        if self.SaoAdj(u, v):
            print(f"Aresta {u} - {v} já existe")
            return None
        
        # Adicionando aresta ao grafo
        self.G.add_edge(u, v, weight=peso)
        return self.G[u][v]

    def RemoverAresta(self, u, v):
        if self.G.has_edge(u, v):
            self.G.remove_edge(u, v)

    def SaoAdj(self, u, v):
        return self.G.has_edge(u, v)

    def N(self, v, Tipo="*"):
        if v < 1 or v > self.n:
            return []
        
        neighbors = list(self.G.neighbors(v))
        if Tipo == "*":
            return neighbors
        else:
            return [w for w in neighbors if self.G[v][w].get("Tipo") == Tipo]

    def setPeso(self, u, v, peso):
        if self.G.has_edge(u, v):
            self.G[u][v]['weight'] = peso

    def getPeso(self, u, v):
        if self.G.has_edge(u, v):
            return self.G[u][v].get('weight', None)
        return None

    def is_leaf(self, v):
        return len(list(self.G.neighbors(v))) == 0
