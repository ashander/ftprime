# Sparse trees

A sparse tree is a sequence $\pi$ of integers, 
with $\pi_k$ denoting the parent of node $k$,
and $0$ denoting no parent (as for example the root).


# Building a sparse tree from coalescent records

Since coalescent records record parent $u$ and children $c$,
to build a tree out of a set of records
one only needs to initialize the entries of $\pi$ all to 0,
and then to set the $c_i$th entries of $\pi$ all equal to $u$ for each $i$.


# Generating trees (the tree iterator)

Suppose we have the tree at $x$ and the next record to the right
begins at $\ell>x$.
The algorithm moves from this tree to the next (at $\ell$)
by first removing all records that end between $x$ and $\ell$
(by setting the parents of each of the children in each record to zero)
and then adding all new records beginning at $\ell$.

![Algorithm T](algorithm_t.png)

Note that this assumes that all relevant entries get zeroed out
when 'removing' records.

## Python version from Jerome

```
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
```
