def trees(records):
    M = len(records)
    I = sorted(range(M), key=lambda j: (records[j].left, records[j].time))
    O = sorted(range(M), key=lambda j: (records[j].right, -records[j].time))
    pi = [-1 for j in range(max(r.node for r in records) + 1)]
    chi = [[] for j in range(max(r.node for r in records) + 1)]
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
            for q in records[h].children:
                pi[q] = records[h].node
            j += 1
        yield pi, chi

if __name__ == "__main__":
    import msprime
    tsgood = msprime.simulate(3, recombination_rate=1, random_seed=111)
    recs = list(tsgood.records())
    for r in recs:
        print(r)
    for t in trees(recs):
        print(t)
