from utils import read_graph

def compare(adj1, adj2):
    try:
        for key in adj1.keys():
            l1, l2 = adj1[key], adj2[key]
            l1.sort()
            l2.sort()
            if l1 != l2: return False
    except KeyError:
        return False

    return True

if __name__ == "__main__":
    adj1 = read_graph('../cmake-build-debug/bin/delaunay.txt')
    adj2 = read_graph('../cmake-build-debug/bin/delaunay_old.txt')
    print(compare(adj1, adj2))