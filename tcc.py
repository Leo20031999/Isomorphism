class Graph:
    def __init__(self):
        self.adj_list = {}

    def add_edge(self, u, v):
        if u not in self.adj_list:
            self.adj_list[u] =[]
        self.adj_list[u].append(v)
        
    def print_graph(self):
        for vertex, adjacent in self.adj_list.items():
            print(f'{vertex}: {"->".join(map(str,adjacent))}')

    def is_leaf(self, u):
        if len(self.adj_list[u]) == 1:
            return True
        else:
            return False

    def center(self):
        lista = {k: list(v) for k, v in self.adj_list.items()}
        leafs = [vertex for vertex in lista if len(lista[vertex]) == 1]
        while len(lista) > 2:
            new_leafs = []
            for leaf in leafs:
                if leaf in lista:
                    neighbour = lista[leaf][0]
                    if leaf in lista[neighbour]:
                        lista[neighbour].remove(leaf)
                    if len(lista[neighbour]) == 1:
                        new_leafs.append(neighbour)
                    del lista[leaf]
            leafs = new_leafs
        return list(lista.keys())

def hausdorffDistanceBetweenTrees(grafo1, grafo2):
    hd = 99
    O = 0
    list_v = grafo1.center()
    v = list_v[0]
    

grafo = Graph()

grafo.add_edge('A', 'B')
grafo.add_edge('B', 'C')
grafo.add_edge('B', 'D')
grafo.add_edge('D', 'E')
grafo.add_edge('D', 'F')

grafo.print_graph()
print(grafo.center())
