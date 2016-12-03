def trees(records):
    M = len(records)
    I = sorted(range(M), key=lambda j: (records[j].left, records[j].time))
    O = sorted(range(M), key=lambda j: (records[j].right, -records[j].time))
    N = 1 + max( max(r.node,max(r.children)) for r in records)
    pi = [-1 for j in range(N)]
    chi = [[] for j in range(N)]
    time = [0.0 for j in range(N)]
    j = 0
    k = 0
    while j < M:
        x = records[I[j]].left
        while records[O[k]].right == x:
            h = O[k]
            print("\tout:", records[h])
            chi[records[h].node] = []
            for q in records[h].children:
                pi[q] = -1
            k += 1

        while j < M and records[I[j]].left == x:
            h = I[j]
            print("\tin:", records[h])
            chi[records[h].node] = records[h].children
            if time[records[h].node]==0.0:
                time[records[h].node]=records[h].time
            elif time[records[h].node]!=records[h].time:
                raise ValueError("Inconsistent node times.")
            for q in records[h].children:
                if time[q] >= time[records[h].node]:
                    raise ValueError("Node times must be stictly increasing within a tree.")
                pi[q] = records[h].node
            j += 1
        yield pi, chi


def parent_dict(pi):
    parent_dict = {}
    for child, parent in enumerate(pi):
        parent_dict[child] = parent
    return parent_dict


if __name__ == "__main__":
    import msprime
    tsgood = msprime.simulate(3, recombination_rate=1, random_seed=111)
    recs = list(tsgood.records())
    for r in recs:
        print(r)
    for t in trees(recs):
        print(t)
