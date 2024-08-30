from structures import graph
class WeightedGraph(graph.Graph):
    def __init__(self):
        super().__init__()
        self.weights = {}
    
    def add_edge(self, u, v, weight):
        super().add_edge(u,v)
        self.weights[(u,v)] = weight
        self.weights[(v,u)] = weight

    def get_weight(self, u, v):
        return self.weights.get((u,v), float('inf'))

    def get_edges(self):
        return self.weights.keys()