class Grafo(object):
    def __init__(self, orientado = False):
        self.n = None
        self.m = None
        self.orientado = orientado
    def DefinirN(self, n):
        self.n = n
        self.m = 0
    def V(self):
        for i in range(1,self.n+1):
            yield i
    
    def E(self, IterarSobreNo = False):
        arestas = set()
        for u in range(1, self.n + 1):
            w = self.L[u].Prox
            while w is not None:
                arestas.add((min(u, w.Viz), max(u, w.Viz)))
                w = w.Prox
        return list(arestas)

