class Graph:
    def __init__(self):
        self.adj_list = {}

    def add_edge(self, u, v):
        if u not in self.adj_list:
            self.adj_list[u] =[]
        self.adj_list[u].append(v)
        
    def print_graph(self):
        allVertices = set(self.adj_list.keys())
        for neighbours in self.adj_list.values():
            allVertices.update(neighbours)
        adjacencyList = {vertex: set() for vertex in allVertices}

        for vertex, neighbours in self.adj_list.items():
            for neighbour in neighbours:
                adjacencyList[vertex].add(neighbour)
                adjacencyList[neighbour].add(vertex)
        for vertex in sorted(adjacencyList):
            adjacent = sorted(adjacencyList[vertex])
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
    
    def calculateHeight(self, node, parent=None):
        height = 0
        for neighbor in self.adj_list.get(node, []):
            if neighbor != parent:
                height = max(height, self.calculateHeight(neighbor, node) + 1)
        return height

    def heights(self):
        heights = {}
        allVertices = set(self.adj_list.keys())
        for vertex in self.adj_list.values():
            allVertices.update(vertex)
        for vertex in allVertices:
            heights[vertex] = self.calculateHeight(vertex,None)
        return heights

def hausdorffDistanceBetweenTrees(grafo1, grafo2):
    return 0   

graph = Graph()
graph.add_edge(0, 1)
graph.add_edge(0, 2)
graph.add_edge(1, 3)
graph.add_edge(1, 4)
graph.add_edge(2, 5)
graph.add_edge(2, 6)
graph.add_edge(5, 7)
graph.add_edge(6, 8)

print("Representação do Grafo:")
graph.print_graph()

heights = graph.heights()
print("\nAlturas dos vértices:")
for vertex, height in heights.items():
    print(f"Vértice {vertex}: Altura {height}")