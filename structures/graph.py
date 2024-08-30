class Graph:
    def __init__(self):
        self.adj_list = {}

    def add_edge(self, u, v):
        if u not in self.adj_list:
            self.adj_list[u] =[]
        if v not in self.adj_list:
            self.adj_list[v] =[]
        self.adj_list[u].append(v)
        self.adj_list[v].append(u)
        
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
    
    def calculateHeight(self):
        heights = {}
        stack = []
        visited = set()
        node_to_children = {}
        node_parents = {}

        for node in self.adj_list:
            if node not in visited:
                visited.add(node)
                stack.append((node,None))
                while stack:
                    current, parent = stack.pop()
                    if current not in node_to_children:
                        node_to_children[current] = []
                        node_parents[parent] = parent
                    for neighbor in self.adj_list.get(current,[]):
                        if neighbor!=parent and neighbor not in visited:
                            visited.add(neighbor)
                            stack.append((neighbor, current))
                            node_to_children[current].append(neighbor)
        
        stack = [node for node, children in node_to_children.items() if not children]

        while stack:
            node = stack.pop()

            if node in heights:
                continue

            parent = node_parents.get(node)
            max_child_height = 0

            for child in node_to_children.get(node,[]):
                if child in heights:
                    max_child_height = max(max_child_height, heights[child])
                else:
                    stack.append(node)
                    break
            else:
                heights[node] = max_child_height + 1
            
            if parent is not None:
                stack.append(parent)
        return heights
        

    def heights(self):
        heights = self.calculateHeight()
        return heights