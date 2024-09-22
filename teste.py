from structures.listaAdj import listaAdj
from structures.grafo import Grafo

def imprimir(grafo):
    for v in grafo.V():
        adjacencias = []
        for w in grafo.N(v):
            vizinho = w.Viz if hasattr(w,'Viz') else w
            adjacencias.append(str(vizinho))
        print(f"{v}: " + ", ".join(adjacencias))

grafo1 = listaAdj(orientado=False)

grafo1.DefinirN(4)

grafo1.AdicionarAresta(1, 2)
grafo1.AdicionarAresta(1, 3)
grafo1.AdicionarAresta(2, 4)
grafo1.AdicionarAresta(3, 4)

imprimir(grafo1)