
#Tree Isomorphism

 #Description:
#   Rooted trees (S,s) and (T,t) are isomorphic
#   iff there is a bijection between childs of s and childs of t
#   such that s.child[i] and t.child[pi[i]] is isomorphic.

#   Two trees are isomorphic iff these are isomorphic with some roots.

# Algorithm:
#   Aho-Hopcroft-Ullmann's algorithm for rooted isomorphism.
#   It simultaneously scans the vertices in the two trees from bottom to root,
#   and assigns unique id (among the level) for each subtrees.
#
#   For unrooted isomorphism, it first finds centers of tree,
#   and try all possibility of rooted tree isomorphism.
#   Since the number of centers in a tree is at most two,
#   it can be performed as the same complexity as rooted isomorphism.

# Complexity:
#   O(n log n). 

# Verified:
#   SPOJ7826.

# References:
#   A. V. Aho, J. E. Hopcroft, and J. D. Ullman (1974):
#   The Design and Analysis of Computer Algorithms.
#   Addison-Wesley.

from structures.Grafo import Grafo

def rooted_isomorphism(S, s, T, t):
    layers_S, parent_S = S.levelize(s)
    layers_T, parent_T = T.levelize(t)
    if len(layers_S) != len(layers_T):
        return False
    
    code_S = {}
    code_T = {}
    long_code_S = {u: [] for u in S.grafo.nodes}
    long_code_T = {u: [] for u in T.grafo.nodes}
    
    for h in reversed(range(len(layers_S))):
        bucket = {}
        
        # Process nodes in S's layer h
        for u in layers_S[h]:
            children = [v for v in S.grafo.neighbors(u) if parent_S.get(v) == u]
            sorted_codes = sorted([code_S[v] for v in children])
            key = tuple(sorted_codes)
            long_code_S[u] = sorted_codes
            bucket[key] = None
        
        # Process nodes in T's layer h
        for u in layers_T[h]:
            children = [v for v in T.grafo.neighbors(u) if parent_T.get(v) == u]
            sorted_codes = sorted([code_T[v] for v in children])
            key = tuple(sorted_codes)
            long_code_T[u] = sorted_codes
            bucket[key] = None
        
        # Assign IDs to each unique key
        keys = sorted(bucket.keys(), key=lambda x: (len(x), x))
        id = 0
        for key in keys:
            bucket[key] = id
            id += 1
        
        # Assign codes for S's nodes
        for u in layers_S[h]:
            key = tuple(long_code_S[u])
            code_S[u] = bucket[key]
        
        # Assign codes for T's nodes
        for u in layers_T[h]:
            key = tuple(long_code_T[u])
            code_T[u] = bucket[key]
    
    return code_S[s] == code_T[t]

def isomorphic(S, T):
    centers_S = S.center()
    centers_T = T.center()
    if len(centers_S) != len(centers_T):
        return False
    if rooted_isomorphism(S, centers_S[0], T, centers_T[0]):
        return True
    if len(centers_S) > 1:
        return rooted_isomorphism(S, centers_S[1], T, centers_T[0])
    return False

def doit():
    print("Bem vindo ao do-it")
    print("Insira a quantidade de vertices: ")
    n = int(input())
    S = Grafo()
    T = Grafo()
    # Build S
    for i in range(n-1):
        print("Insira as arestas (u,v): ")
        u, v = map(int, input().split())
        S.adicionar_aresta(u-1, v-1)
    print("Grafo S feito!")
    # Build T
    for i in range(n-1):
        print("Insira as arestas (u,v): ")
        u, v = map(int, input().split())
        T.adicionar_aresta(u-1, v-1)
    print("Grafo T feito!")
    if isomorphic(S, T):
        print("YES")
    else:
        print("NO")

def main():
    print("Insira a quantidade de casos: ")
    ncase = int(input())
    for _ in range(ncase):
        doit()

if __name__ == "__main__":
    main()