from structures.grafo import Grafo

class listaAdj(Grafo):
    class NoAresta(object):
        def __init__(self):
            self.Viz = None
            self.e = None
            self.Prox = None
            self.Peso = 0
    
    class Aresta(object):
        def __init__(self):
            self.v1, self.No1 = None, None
            self.v2, self.No2 = None, None
            self.Peso = 0 

    def DefinirN(self, n, VizinhancaDuplamenteLigada = False):
        super(listaAdj, self).DefinirN(n)
        self.L = [None]*(self.n+1)
        for i in range(1, self.n+1):
            self.L[i] = listaAdj.NoAresta()
        self.VizinhancaDuplamenteLigada = VizinhancaDuplamenteLigada

    def AdicionarAresta(self, u, v, peso = 0):
        if self.SaoAdj(u,v):
            return None
        def AdicionarLista(u, v, e, Tipo):
            No = listaAdj.NoAresta()
            No.Viz,No.e,No.Prox, No.Peso = v, e, self.L[u].Prox, peso
            self.L[u].Prox=No
            return No
        
        e = listaAdj.Aresta()
        e.v1, e.v2 = u, v
        e.No1 = AdicionarLista(u,v,e,"+")
        if not self.orientado:
            if not self.SaoAdj(v,u):
                e.No2 = AdicionarLista(v,u,e,"-")
        e.Peso = peso
        self.m = self.m+1
        return e
        
    def RemoverAresta(self, uv):
        def RemoverLista(No):
            No.Ant.Prox = No.Prox
            if No.Prox != None:
                No.Prox.Ant = No.Ant
        RemoverLista(uv.No1)
        RemoverLista(uv.No2)
    
    def SaoAdj(self, u, v):
        for w in self.N(u):
            if w==v:
                return True
        return False
        
    def N(self, v , Tipo = "*", Fechada = False, IterarSobreNo=False):
        if Fechada:
            No = listaAdj.NoAresta()
            No.Viz, No.e, No.Prox = v, None, None
            yield No if IterarSobreNo else No.Viz
        w = self.L[v].Prox
        while w != None:
            if Tipo == "*" or w.Tipo == Tipo:
                yield w if IterarSobreNo else w.Viz
            w = w.Prox
    
    def getPeso(self, u, v):
        for w in self.N(u):
            if w == v:
                return self.L[u].Prox.Peso
        return None
    
    def is_leaf(self, v):
        return self.L[v].Prox is None
    
    def preorder(self,v,visit):
        pilha = [v]
        while pilha:
            atual = pilha.pop()
            visit(atual)
            for w in reversed(list(self.N(atual))):
                pilha.append(w)
        