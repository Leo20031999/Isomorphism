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
        for v in self.V():
            for w in self.N(v, Tipo = "+" if self.orientado else "*", IterarSobreNo = IterarSobreNo):
                enumerar = True
                if not self.orientado:
                    wint = w if isinstance(w,int) else w.Viz
                    enumerar = v < wint
            if enumerar:
                yield(v,w)


