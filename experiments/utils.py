
def add(adj, a, b):
    if not a in adj.keys():
        adj[a] = []
    adj[a].append(b)
    
def read_graph(filename):
    adj = {}
    with open(filename, 'r') as f:
        for line in f:
            a, b = [int(x) for x in line.split(' ')]
            add(adj, a, b)
            add(adj, b, a)
    return adj