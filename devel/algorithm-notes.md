# Requirements of coalescence records in a tree sequence

1. Offspring must be born after their parents (and hence, no loops).
2. The set of intervals on which individual $a$ is a child must be disjoint, for every $a$.
3. The set of intervals on which individual $a$ is a parent must be disjoint, for every $a$.
4. All records with the same parent must occur at the same time.
5. The samples must be numbered 0,...,n-1, and the smallest non-sampled label must be an internal node.

The first two disallow time travel and multiple inheritance;
the third and fifth are algorithmic requirements; 
and the fourth implies that we measure branch lengths in clock time
(and hence get ultrametric trees).


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

## Questions:

1.  Can we have partially observed tips that are not samples?

    * **Yes:** if the largest ID is an internal node (or we change how memory is allocated)

2.  Can we have the same internal node be the parent of more than on coalescent record
    at different times?

    * **No.**


## Python version from Jerome

... slightly modified.

```{python}

def trees(records):
    M = len(records)
    I = sorted(range(M), key=lambda j: (records[j].left, records[j].time))
    O = sorted(range(M), key=lambda j: (records[j].right, -records[j].time))
    N = 1 + max( max(r.node,max(r.children)) for r in records)
    pi = [-1 for j in range(N)]
    chi = [[] for j in range(N)]
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

# Clock time or meioses

Suppose that:

1. We start with individual $a$ at $t=0$.
2. $a$ gives birth to $b$ at $t=1$.
3. $b$ gives birth to $c$ at $t=2$.
4. $a$ gives birth to $d$ and $e$ at $t=3$.
5. $a$ and $b$ die and we sample $c$, $d$, and $e$ at $t=4$.

What should the tree be?
If we count in **meioses**, it should be

```
                                        
2            a _                        
            / \ \                       
1          b   d e                      
          /                             
0        c                              
```

Here, nodes = individuals.

On the other hand, if we count in **clock time**, then it should be

```
                                        
4       a                               
        |\                              
3       a b                             
        | |\                            
2       a b c                           
       /|\ \ \                          
1     a d e b c                         
        | |   |                         
0       d e   c                         
```

Here, nodes = (individuals at a particular point in time)



# Forwards-time algorithms

Requirements:

1.  Coalescence records have at least two offspring.
2.  Coalescence records are output in time-order.
3.  We need to be able to re-label so samples have the first $n$ labels at the end.



## Without recombination

This algorithm will be in **clock time**.

Now:
the idea is that each branch of the tree gets a label;
so we record coalescence records then a branch splits.
As we move forwards in time, we keep track of:

-  `L` : a vector of unique labels, one for each individual currently alive

and we write out coalescent records, which are:

- `[parent, offspring, time]` : parent and offspring are labels of branch tips; 
    the `time` is the time of *birth of the offspring*,
    so we measure time in *clock time*.

### Example

Here's the example above, extended slightly:

1. We start with individual $a$ at $t=0$.
2. $a$ gives birth to $b$ at $t=1$.
3. $b$ gives birth to $c$ at $t=2$.
4. $a$ gives birth to $d$ and $e$ at $t=3$.
5. $b$ gives birth to $f$ and $a$ and $b$ die at $t=4$
6. We sample $c$, $d$, $e$, and $f$.

With lineage labels on the right, and time now moving forwards:
```
------------
time   tree        |   lineages      | current state |  records output

0       a          |       0         |  a:0          |              

------------

0       a          |       0         |               |  [ 0, (1,2), 0 ]               
        |\         |       |\        |  a:1          |                 
1       a b        |       1 2       |  b:2          |  

------------

0       a          |       0         |               |  
        |\         |       |\        |               |                 
1       a b        |       1 2       |  a:1          |  [ 2, (3,4), 1 ]               
        | |\       |       | |\      |  b:3          |                 
2       a b c      |       1 3 4     |  c:4          |  

------------
                                     |               |   
0       a          |       0         |               |  
        |\         |       |\        |               |                 
1       a b        |       1 2       |   a:5         |  
        | |\       |       | |\      |   b:3         |                 
2       a b c      |       1 3 4     |   c:4         |  [ 1, (5,6,7), 2 ]               
       /|\ \ \     |      /|\ \ \    |   d:6         |                 
3     a d e b c    |     5 6 7 3 4   |   e:7         |                 

------------
                                     |               |   
0       a          |       0         |               |  
        |\         |       |\        |               |                 
1       a b        |       1 2       |               |  
        | |\       |       | |\      |               |                 
2       a b c      |       1 3 4     |               |  
       /|\ \ \     |      /|\ \ \    |    c:4        |                 
3     a d e b c    |     5 6 7 3 4   |    d:6        |                 
      | | | |\ \   |       | |  \ \  |    e:7        |                 
4     * d e * f c  |       6 7   3 4 |    f:3        |  

```

Here we noticed that b died and was replaced by f,
so kept just one lineage;
otherwise we'd have:

```
0       a          |       0         |               |  
        |\         |       |\        |               |                 
1       a b        |       1 2       |               |  
        | |\       |       | |\      |               |                 
2       a b c      |       1 3 4     |    c:4        |  
       /|\ \ \     |      /|\ \ \    |    d:6        |                 
3     a d e b c    |     5 6 7 3 4   |    e:7        |                 
      | | | |\ \   |       | | |\ \  |    f:9        |   [ 3, (8,9), 4 ]              
4     * d e * f c  |       6 7 8 9 4 |               |  

```

... where we only really have to output the coalescence record 
so as to remember that lineage 9 is a continuation of lineage 3.
(But, either works.)


### Algorithm, no recombination


At each time step, some individuals die and some give birth.
Let $L$ denote the current vector of labels, one for each extant individual.

0.  Begin with $L$ a vector of unique labels.
1.  For each birth and/or death at time $t$, consider the parent:
    
    a.  If they die with no offspring, remove them.
    b.  If they die and leave one offspring,
        let their offspring inherit their label.
    c.  Otherwise, assign a new label to each offspring,
        and to *themselves*, if they survive.
        Let the label of the parent be $\ell_0$ and these new labels be $\ell_1, \ldots, \ell_k$,
        and output the coalescence record
        $$( \ell_0, (\ell_1, \ldots, \ell_k), t ).$$


Check: every time a new label is assigned, 
it should appear in a coalescence record (so we know who it's parent is).


<!--
## Jerome's outline of forward-time alg

edited from his email per correction

I think you're on the right track here, but I'm a bit confused as to why you're
recording 'coalescences' where there is a single child. In my mind, the
simplest way to think about this is to just imagine running the coalescent
algorithm in reverse, using the same data structures. Suppose we were doing
this for a single locus, and we started with ancestral individuals `[0, 1, 2, 3, 4, 5]`. Now, we generate some children randomly from this, and we get

`[3, 2, 1, 3, 3, 0]`

Since multiple children have chosen the same parent, we must have a
coalescence. Replace the coalesced node with new IDs, giving

`[6, 2, 1, 7, 8, 0]`

and record the coalescence (6, 7, 8) -> 3. Repeating, we get the next generation

`[1, 1, 7, 8, 0, 0]`

Because we have two distinct parents here, we have two coalescences:

`[9, 10, 7, 8, 11, 12]`

with records `(9, 10) -> 1` and `(11, 12) -> 0`.

The key thing here is that we only record _coalescences_, where actual changes
to the topology of the trees occur. If only one child in the next generation
inherits from a particular individual, then its node ID is passed on unchanged.
When coalescences occur, we allocate new nodes for the children, so that they
each become nodes in the genealogy.

A nice property of this representation is that it's easy to tell when we've
simulated far enough forward in time to be sure that there is an MRCA: if you
started with `n` individuals, and any `ID < n` remains in the population, then you
need to continue (thinking about this again, it may not be a sufficient
condition though).

Once you get to the end, you should be able to rename your present day
population as `0` to `n - 1`, update any affected records, and then use the tree
sequence directly. Calling subset on this should then trim out all the stuff
that you don't want. (Warning though: there is still a bug in subset when
dealing with weird corner cases. If you hit a 'Bad records' when running subset
it's because of this bug. Let me know if you do.)

This is nice and easy for a single locus, and it's a lot more tricky for
multiple loci, of course. However, I think the same strategy as used in
`msprime` will work, where each individual is a linked list of ancestral
segments. The algorithms listed in the paper for manipulating these linked
lists should also work here, as we're doing essentially the same thing. I think
this will probably end up being equivalent to the approach that you're
outlining, but it may be a bit simpler.

I would try to get things working for a single locus model first though, and
see if the basic logic holds.
 -->

## With recombination

We need only make sure we record from where each individual has inherited her genome.
However, coalescence records focus on the "other end" of the relationship, siblingship.

Here's a minimal example:
with `(i,j,x)->k` denoting that individual `k` inherits from `i` on `[0,x)` and from `j` on `[x,1)`:

1. Begin with an individual `a` (and another anonymous one) at `t=0`.
2. `(a,?,1.0)->b` at `t=1`
3. `(a,b,0.5)->c` at `t=2`
4. `(a,c,0.2)->d` and `(c,b,0.6)->e` at `t=3`.
5. `a` and `b` die, we sample `c,d,e` at `t=4`.


### Algorithm No. 1

As a first pass, 
we'll "deal with death last" at each time point,
so that everyone has their parent as a "sibling".
Coalescent records are
```
(left, right, node, (children), time)
```
The current state is a list of labels,
with each label corresponding to a currently alive individual.


```
t  |   trees       |             lineages                                    |    state           |   records output                                              
                   |                                                                                             
0  |        a      |                         1                               |    1:a             |                        

# (a,?,1.0)->b , t=1
1  |        a      |                         1                               |                    |                        
   |       / \     |                        / \                              |    2:a             |   ( 0.0, 1.0, 1, (2,3), 1 )    
   |      a   b    |                       2   3                             |    3:b             |                        

# (a,b,0.5)->c , t=2
 2 |        a      |             1                           1               |                    |                        
   |       / \     |            / \                         / \              |                    |                        
   |      a   b    |           2   3                       2   3             |    4:a             |   ( 0.0, 0.5, 2, (4,6), 2 )
   |     / \ / \   |          / \   \                     /   / \            |    5:b             |   ( 0.5, 1.0, 3, (5,6), 2 )
   |    a   c   b  |         4   6   5                   4   6   5           |    6:c             |                         
   |               |                                                         |                    |                         
   |               |         [0.0,0.5)                    [0.5,1.0)          |                    |                         
                                                                    
# (a,c,0.2)->d , t=3                                                
3  |       a       |       1             1                    1              |                    |                         
   |      / \      |      / \           / \                  / \             |                    |                         
   |     a   b     |     2   3         2   3                2   3            |                    |   ( 0.0, 0.2, 4, (7,8), 3 )
   |    / \ / \    |    / \   \       / \   \              /   / \           |    5:b             |                         
   |   a   c   b   |   4   6   5     4   6   5            4   6   5          |                    |   ( 0.2, 1.0, 6, (8,9), 3 )    ** see below
   |  / \ /        |  / \           /   /|               /   /|              |    7:a             |                        
   | a   d         | 7   8         7   8 9              7   8 9              |    8:d             |                        
   |               |                                                         |    9:c             |                        
   |               |                                                         |                    |                        
   |               |  [0.0,0.2)      [0.2,0.5)             [0.5,1.0)         |                    |                        

# (c,b,0.6)->e , still t=3
3  |       a       |       1             1           1              1        |                    |                        
   |      / \      |      / \           / \         / \            / \       |                    |                        
   |     a   b     |     2   3         2   3       2   3          2   3      |                    |                        
   |    / \ / \    |    / \   \       / \   \     /   / \        /   / \     |                    |   ( 0.0, 0.2, 6, (9,10), 3 )
   |   a   c   b   |   4   6   5     4   6   5   4   6   5      4   6   5    |    7:a             |   ( 0.2, 0.6, 6, (8,9,10), 3 )   ** split
   |  / \ / \ / \  |  /|   |\   \    |  /|\   \  |  /|\   \    /   /|   |\   |    8:d             |   ( 0.6, 1.0, 6, (8,9), 3 )      ** split
   | a   d   e   b | 7 8   9 10 11   7 8 9 10 11 7 8 9 10 11  7   8 9  10 11 |    9:c             |                        
   |               |                                                         |   10:e             |                                   
   |               |  [0.0,0.2)      [0.2,0.5)    [0.5,0.6)      [0.6,1.0)   |   11:b             |   ( 0.6, 1.0, 5, (10,11), 3 )
                                                                           
# a,b die; sample c,d,e at t=4
3  |       a       |       1             1           1              1        |                    |                        
   |      / \      |      / \           / \         / \            / \       |                    |  
   |     a   b     |     2   3         2   3       2   3          2   3      |                    |  
   |    / \ / \    |    / \   \       / \   \     /   / \        /   / \     |                    |  
   |   a   c   b   |   4   6   5     4   6   5   4   6   5      4   6   5    |                    |  
   |  / \ /|\ / \  |  /|   |\   \    |  /|\   \  |  /|\   \    /   /|   |\   |    8:d             |  
   | a   d | e   b | 7 8   9 10 11   7 8 9 10 11 7 8 9 10 11  7   8 9  10 11 |                    |  
   |    /  |  \    |   |   |  \       /  |  \     /  |  \        /  |   |    |    9:c             |  
   |   d   c   e   |   8   9   10    8   9   10  8   9   10     8   9  10    |   10:e             |  
   |               |                                                         |                    |  
   |               |  [0.0,0.2)      [0.2,0.5)    [0.5,0.6)      [0.6,1.0)   |                    |  
   |               |                                                         |                    |  

````


What'd we do there?
Every time there's a new individual, they get a lineage all to themselves (their whole genome).
When anyone has offspring
we create new lineages for the parents as well; 
these new lineages serve as siblings to their offspring in the new coalescent records.
Sometimes in building coalescent records for events happening in the same generation,
we modified one already in the pipeline.

Here's the order we add and remove records to create the trees
in Algorithm T:

```
  Adding                           |   Removing
  ------                           |   --------
1  ( 0.0, 1.0, 1, (2,3),    1 )    |
1  ( 0.0, 0.5, 2, (4,6),    2 )    |
1  ( 0.0, 0.2, 4, (7,8),    3 )    |
1  ( 0.0, 0.2, 6, (9,10),   3 )    |
                                   | 2  ( 0.0, 0.2, 4, (7,8),   3 )
                                   | 2  ( 0.0, 0.2, 6, (9,10),  3 )
2  ( 0.2, 0.6, 6, (8,9,10), 3 )    |
                                   | 3  ( 0.0, 0.5, 2, (4,6),   2 )
3  ( 0.5, 1.0, 3, (5,6),    2 )    |
                                   | 4  ( 0.2, 0.6, 6, (8,9,10),3 )
4  ( 0.6, 1.0, 6, (8,9),    3 )    |
4  ( 0.6, 1.0, 5, (10,11),  3 )    |
                                   | *  ( 0.6, 1.0, 6, (8,9),   3 )
                                   | *  ( 0.6, 1.0, 5, (10,11), 3 )
                                   | *  ( 0.5, 1.0, 3, (5,6),   2 )
                                   | *  ( 0.0, 1.0, 1, (2,3),   1 )    
```


**Algorithm:**

Let $L$ be the dictionary of labels (initiallized to a unique set);
and let $u$ be a list of flags denoting whether we've processed a reproduction for the individual
this generation.
At time point $t$:

0.  Set every element of $u$ to False.
1.  For each birth, with parents $i$, $j$ and offspring $k$:
    
    a.  If $u(i)$ is False, 
        assign a new label to $i$, and set $u(i)$ to True,
        and likewise for $j$.
        Let $L'(i)$ be the previous label for $i$ 
        and $L(i)$ the new label for $i$ (possibly already assigned), 
        and likewise for $j$;
        and let $L(k)$ be the label of the offspring.

    b.  Let $x_0=0 \le x_1,x_2,\ldots,x_n \le x_{n=1}=1$ be the recombination breakpoints,
        with the offspring inheriting from $i$ on $[x_{2\ell},x_{2\ell+1})$
        and from $j$ on $[x_{2\ell+1},x_{2\ell+2})$;
        and output coalescence records
        $(x_{2\ell},x_{2\ell+1},L'(i),(L(i),L(k)),t)$
        and
        $(x_{2\ell+1},x_{2\ell+2},L'(j),(L(j),L(k)),t)$
        for each $\ell$.

    c.  Remove dead individuals.

    d.  Merge-and-split coalescence records that share a parent.



## Renumbering

At the end we need to renumber,
so that the samples have labels $1,\ldots,n$
and the largest-numbered node is an internal node.



<!--
# An alternative algorithm (probably ignore this)


The first algorithm has an inefficiency
in that it records coalescent events for EVERY new offspring.
Here is an algorithm that initially keeps track of *more*,
but may be more amenable to subsequent pruning of redundant information.

The current state is a list of labels,
with each label corresponding to a segment of a currently alive individual.


```
t  |   trees       |             lineages                                    |    state           |   records output                                              
                   |                                                                                             
0  |        a      |                         1                               |    1:a,[0,1)       |                        

# (a,?,1.0)->b , t=1
1  |        a      |                         1                               |                    |                        
   |       / \     |                        / \                              |    2:a,[0,1)       |   ( 0.0, 1.0, 1, (2,3), 1 )    
   |      a   b    |                       2   3                             |    3:b,[0,1)       |                        

# (a,b,0.5)->c , t=2
 2 |        a      |             1                           1               |    2:a,[0.5,1.0) * |                        
   |       / \     |            / \                         / \              |    3:b,[0.0,0.5) * |                        
   |      a   b    |           2   3                       2   3             |    4:a,[0.0,0.5)   |   ( 0.0, 0.5, 2, (4,6), 2 )
   |     / \ / \   |          / \   \                     /   / \            |    5:b,[0.5,1.0)   |   ( 0.5, 1.0, 3, (5,6), 2 )
   |    a   c   b  |         4   6   3                   2   6   5           |    6:c,[0.0,1.0)   |                         
   |               |                                                         |                    |                         
   |               |         [0.0,0.5)                    [0.5,1.0)          |                    |                         
                                                                    
# (a,c,0.2)->d , t=3                                                
3  |       a       |       1             1                    1              |    2:a,[0.5,1.0)   |                         
   |      / \      |      / \           / \                  / \             |    3:b,[0.0,0.5)   |                         
   |     a   b     |     2   3         2   3                2   3            |    4:a,[0.2,0.5) * |   ( 0.0, 0.2, 4, (7,8), 3 )
   |    / \ / \    |    / \   \       / \   \              /   / \           |    5:b,[0.5,1.0)   |                         
   |   a   c   b   |   4   6   3     4   6   3            2   6   5          |    6:c,[0.0,0.2) * |   ( 0.2, 1.0, 6, (8,9), 3 )    ** see below
   |  / \ /        |  / \           /   / \                  / \             |    7:a,[0.0,0.2)   |                        
   | a   d         | 7   8         4   8   9                8   9            |    8:d,[0.0,1.0)   |                        
   |               |                                                         |    9:c,[0.2,1.0)   |                        
   |               |                                                         |                    |                        
   |               |  [0.0,0.2)      [0.2,0.5)             [0.5,1.0)         |                    |                        

# (c,b,0.6)->e , still t=3
3  |       a       |       1             1           1              1        |    2:a,[0.5,1.0)   |                        
   |      / \      |      / \           / \         / \            / \       |    3:b,[0.0,0.5)   |                        
   |     a   b     |     2   3         2   3       2   3          2   3      |    4:a,[0.2,0.5)   |                        
   |    / \ / \    |    / \   \       / \   \     /   / \        /   / \     |    5:b,[0.5,0.6) * |   ( 0.0, 0.2, 6, (10,11), 3 )
   |   a   c   b   |   4   6   3     4   6   3   2   6   5      2   6   5    |    7:a,[0.0,0.2)   |   ( 0.2, 0.6, 6, (8,9,11), 3 )   ** split
   |  / \ / \ / \  |  /|   |\       /   /|\         /|\            /|  / \   |    8:d,[0.0,1.0)   |   ( 0.6, 1.0, 6, (8,9), 3 )      ** split
   | a   d   e   b | 7 8  10 11    4   8 9 11      8 9 11         8 9 11 12  |    9:c,[0.2,1.0)   |                        
   |               |                                                         |   10:c,[0.0,0.2)   |                                   
   |               |  [0.0,0.2)      [0.2,0.5)    [0.5,0.6)      [0.6,1.0)   |   11:e,[0.0,1.0)   |                       
   |               |                                                         |   12:b,[0.6,1.0)   |   ( 0.6, 1.0, 5, (11,12), 3 )
                                                                           
# a,b die at t=4; sample c,d,e
4  |       a       |       1             1           1              1        |                    | 
   |      / \      |      / \           / \         / \            / \       |                    | 
   |     a   b     |     2   3         2   3       2   3          2   3      |     8:d,[0.0,1.0)  | 
   |    / \ / \    |    / \   \       / \   \     /   / \        /   / \     |     9:c,[0.2,1.0)  | 
   |   a   c   b   |   4   6   3     4   6   3   2   6   5      2   6   5    |    10:c,[0.0,0.2)  | 
   |  / \ /|\ / \  |  /|   |\       /   /|\         /|\            /|  / \   |    11:e,[0.0,1.0)  | 
   | a   d | e   b | 7 8  10 11    4   8 9 11      8 9 11         8 9 11 12  |                    | 
   |    /  |  \    |   |   |  \       /  |  \     /  |  \        /  | |      |                    | 
   |   d   c   e   |   8  10   11    8   9   11  8   9   11     8   9 11     |                    | 
   |               |                                                         |                    | 
   |               |  [0.0,0.2)      [0.2,0.5)    [0.5,0.6)      [0.6,1.0)   |                    | 
   |               |                                                         |                    | 

````


What'd we do there?
Every time there's a new individual, they get a lineage all to themselves (their whole genome).
When anyone has offspring
we create new lineages for the parents as well, 
but only the parts of their genomes inherited by the offspring;
these new lineages serve as siblings to their offspring in the new coalescent records.

Here's the order we add and remove records to create the trees
in Algorithm T:

```
   Adding                       |      Removing
1 ( 0.0, 1.0, 1, (2,3),   1 )   |   2 ( 0.0, 0.2, 4, (7,8),   3 )
1 ( 0.0, 0.5, 2, (4,6),   2 )   |   2 ( 0.0, 0.2, 6, (10,11), 3 )
1 ( 0.0, 0.2, 4, (7,8),   3 )   |   3 ( 0.0, 0.5, 2, (4,6),   2 )
1 ( 0.0, 0.2, 6, (10,11), 3 )   |   4 ( 0.2, 0.6, 6, (8,9,11),3 )
2 ( 0.2, 0.6, 6, (8,9,11),3 )   |   * ( 0.6, 1.0, 6, (8,9),   3 )
3 ( 0.5, 1.0, 3, (5,6),   2 )   |   * ( 0.6, 1.0, 5, (11,12), 3 )              
4 ( 0.6, 1.0, 6, (8,9),   3 )   |   * ( 0.5, 1.0, 3, (5,6),   2 )
4 ( 0.6, 1.0, 5, (11,12), 3 )   |   * ( 0.0, 1.0, 1, (2,3),   1 )                        
```
-->

# Algorithm M: including single-offspring events.

Here's a minimal example:
with `(i,j,x)->k` denoting that individual `k` inherits from `i` on `[0,x)` and from `j` on `[x,1)`:

1. Begin with an individual `a` (and another anonymous one) at `t=0`.
2. `(a,?,1.0)->b` at `t=1`
3. `(a,b,0.5)->c` at `t=2`
4. `(a,c,0.2)->d` and `(c,b,0.6)->e` at `t=3`.
5. `a` and `b` die, we sample `c,d,e` at `t=4`.


Coalescent records are
```
(left, right, node, (children), time)
```

The current state will be a list of coalescent records for which extant individuals are the parents;
these are output when the individual dies.


```
t  |   trees       |             lineages                                    |    state                         |   records output                                              
                   |                                                                                                           
0  |        a      |                         a                               |                                  |                        

# (a,?,1.0)->b , t=1
1  |        a      |                         a                               |                                  |                        
   |       / \     |                        / \                              |    ( 0.0, 1.0, a, (b,), 1 )      |  
   |      a   b    |                       a   b                             |                                  |                        

# (a,b,0.5)->c , t=2
 2 |        a      |             a                           a               |                                  |                        
   |       / \     |            / \                         / \              |                                  |                        
   |      a   b    |           a   b                       a   b             |                                  |   
   |     / \ / \   |          / \   \                     /   / \            |                                  |   
   |    a   c   b  |         a   c   b                   a   c   b           |    ( 0.0, 0.5, a, (b,c), 1 )     |                         
   |               |                                                         |    ( 0.5, 1.0, a, (b,),  1 )     |                         
   |               |         [0.0,0.5)                    [0.5,1.0)          |    ( 0.5, 1.0, b, (c,),  2 )     |                         
                                                                    
# (a,c,0.2)->d , t=3                                                
3  |       a       |       a             a                    a              |                                  |                         
   |      / \      |      / \           / \                  / \             |                                  |                         
   |     a   b     |     a   b         a   b                a   b            |     ( 0.0, 0.2, a, (b,c,d), 1 )  |
   |    / \ / \    |    / \   \       / \   \              /   / \           |     ( 0.2, 0.5, a, (b,c),   1 )  |
   |   a   c   b   |   a   c   b     a   c   b            a   c   b          |     ( 0.5, 1.0, a, (b,),    1 )  |  <-- this one gets split next time
   |  / \ /        |  / \           /   /|               /   /|              |     ( 0.5, 1.0, b, (c,),    2 )  |                        
   | a   d         | a   d         a   d c              a   d c              |     ( 0.2, 0.8, c, (d,),    3 )  |  <-- so does this one
   |               |                                                         |                                  |                        
   |               |                                                         |                                  |                        
   |               |  [0.0,0.2)      [0.2,0.5)             [0.5,1.0)         |                                  |                        

# (c,b,0.6)->e , still t=3
3  |       a       |       a             a           a              a        |                                  |                        
   |      / \      |      / \           / \         / \            / \       |     ( 0.0, 0.2, a, (b,c,d), 1 )  |
   |     a   b     |     a   b         a   b       a   b          a   b      |     ( 0.2, 0.5, a, (b,c),   1 )  |
   |    / \ / \    |    / \   \       / \   \     /   / \        /   / \     |     ( 0.5, 1.0, a, (b,),    1 )  |
   |   a   c   b   |   a   c   b     a   c   b   a   c   b      a   c   b    |     ( 0.5, 0.6, b, (c,),    2 )  |
   |  / \ / \ / \  |  /|   |\   \    |  /|\   \  |  /|\   \    /   /|   |\   |     ( 0.6, 1.0, b, (c,e),   2 )  |
   | a   d   e   b | a d   c e   b   a d c e   b a d c e   b  a   d c   e  b |     ( 0.0, 0.2, c, (e,),    3 )  |  <-- this is a 'new' record
   |               |                                                         |     ( 0.2, 0.6, c, (d,e),   3 )  |
   |               |  [0.0,0.2)      [0.2,0.5)    [0.5,0.6)      [0.6,1.0)   |     ( 0.6, 0.8, c, (d,),    3 )  |  <-- so is this
                                                                           
# a,b die at t=4
4  |       a       |       a             a           a              a        |                                  |                        
   |      / \      |      / \           / \         / \            / \       |                                  |  ( 0.0, 0.2, a, (b,c,d), 1 )
   |     a   b     |     a   b         a   b       a   b          a   b      |                                  |  ( 0.2, 0.5, a, (b,c),   1 )
   |    / \ / \    |    / \   \       / \   \     /   / \        /   / \     |                                  |  ( 0.5, 1.0, a, (b,),    1 )
   |   a   c   b   |   a   c   b     a   c   b   a   c   b      a   c   b    |                                  |  ( 0.5, 0.6, b, (c,),    2 )
   |  / \ /|\ / \  |  /|   |\   \    |  /|\   \  |  /|\   \    /   /|   |\   |                                  |  ( 0.6, 1.0, b, (c,e),   2 )
   | a   d | e   b | a d   c e   b   a d c e   b a d c e   b  a   d c   e b  |     ( 0.0, 0.2, c, (e,),    3 )  |
   |    /  |  \    |   |   |  \       /  |  \     /  |  \        /  |   |    |     ( 0.2, 0.6, c, (d,e),   3 )  |
   |   d   c   e   |   d   c   e     d   c   e   d   c   e      d   c   e    |     ( 0.6, 0.8, c, (d,),    3 )  |
   |               |                                                         |                                  |
   |               |  [0.0,0.2)      [0.2,0.5)    [0.5,0.6)      [0.6,1.0)   |                                  |

# sample c,d,e still at t=4
4  |       a       |       a             a           a              a        |                                  |                        
   |      / \      |      / \           / \         / \            / \       |                                  |
   |     a   b     |     a   b         a   b       a   b          a   b      |                                  |
   |    / \ / \    |    / \   \       / \   \     /   / \        /   / \     |                                  |
   |   a   c   b   |   a   c   b     a   c   b   a   c   b      a   c   b    |                                  |
   |  / \ /|\ / \  |  /|   |\   \    |  /|\   \  |  /|\   \    /   /|   |\   |                                  |
   | a   d | e   b | a d   c e   b   a d c e   b a d c e   b  a   d c   e b  |                                  |  ( 0.0, 0.2, c, (e,),    3 )
   |    /  |  \    |   |   |  \       /  |  \     /  |  \        /  |   |    |                                  |  ( 0.2, 0.6, c, (d,e),   3 )
   |   d   c   e   |   d   c   e     d   c   e   d   c   e      d   c   e    |                                  |  ( 0.6, 0.8, c, (d,),    3 )
   |               |                                                         |                                  |
   |               |  [0.0,0.2)      [0.2,0.5)    [0.5,0.6)      [0.6,1.0)   |                                  |

````

What happens there?

- Every time we produce a new offspring
    we check for any overlapping records that share a parent with that offspring,
    and merge the new 'record' with the existing one(s).

- When individuals die we output all records for which they are the parent.


# Another example situation:

With `(i,j,x)->k` denoting that individual `k` inherits from `i` on `[0,x)` and from `j` on `[x,1)`:

1. Begin with an individual `a` (and another anonymous one) at `t=0`.
2. `(a,?,1.0)->b` and `(a,?,1.0)->c` at `t=1`
3. `(b,a,0.9)->d` and `(a,c,0.1)->e` and then `a` dies at `t=2`
4. `(d,e,0.7)->f` at `t=3`
5. `(f,d,0.8)->g` and `(e,f,0.2)->h` at `t=4`.
6. `(b,g,0.6)->i` and `(g,h,0.5)->j` and `(h,c,0.4)->k` at `t=5`.
7. We sample `i` and `j`.


Here are the trees:
```
t                  |              |              |             |             |             |             |             |             |            
                                                                                                                                                  
0       --a--      |     --a--    |     --a--    |    --a--    |    --a--    |    --a--    |    --a--    |    --a--    |    --a--    |    --a--   
       /  |  \     |    /  |  \   |    /     \   |   /     \   |   /     \   |   /     \   |   /     \   |   /     \   |   /     \   |   /  |  \  
1     b   |   c    |   b   |   c  |   b       c  |  b       c  |  b       c  |  b       c  |  b       c  |  b       c  |  b       c  |  b   |   c 
      |\ / \ /|    |   |\   \  |  |   |\     /|  |  |\     /|  |  |\     /   |  |\     /   |   \     /   |   \     /   |   \     /   |     /   /  
2     | d   e |    |   | d   e |  |   | d   e |  |  | d   e |  |  | d   e    |  | d   e    |    d   e    |    d   e    |    d   e    |    d   e   
      | |\ /| |    |   |  \  | |  |   |  \  | |  |  |  \    |  |  |  \       |  |  \       |     \       |       /     |    |  /     |    |  /    
3     | | f | |    |   |   f | |  |   |   f | |  |  |   f   |  |  |   f      |  |   f      |      f      |      f      |    | f      |    | f     
      | |/ \| |    |   |  /  | |  |   |  /  | |  |  |  / \  |  |  |  / \     |  |  / \     |     / \     |     / \     |    |  \     |    |  \    
4     | g   h |    |   | g   h |  |   | g   h |  |  | g   h |  |  | g   h    |  | g   h    |    g   h    |    g   h    |    g   h    |    g   h   
      |/ \ / \|    |   |  \    |  |   |  \    |  |  |  \    |  |  |  \   \   |  |    / \   |   /   / \   |   /   / \   |   /   / \   |   /   / \  
5     i   j   k    |   i   j   k  |   i   j   k  |  i   j   k  |  i   j   k  |  i   j   k  |  i   j   k  |  i   j   k  |  i   j   k  |  i   j   k 
                                                                                                                                                  
                   |   0.0 - 0.1  |   0.1 - 0.2  |  0.2 - 0.4  |  0.4 - 0.5  |  0.5 - 0.6  |  0.6 - 0.7  |  0.7 - 0.8  |  0.8 - 0.9  |  0.9 - 1.0 

```


Here is one labeling of the lineages of i, j and k:
```
t                  |              |              |             |             |             |             |             |             |            
                                                                                                                                                  
5       --a--      |     --6--    |     --6--    |    --6--    |    --.--    |    --.--    |    --.--    |    --.--    |    --6--    |    --6--   
       /  |  \     |    /  |  \   |    /     \   |   /     \   |   /     \   |   /     \   |   /     \   |   /     \   |   /     \   |   /  |  \  
4     b   |   c    |   5   |   2  |   5       2  |  5       2  |  5       .  |  5       .  |  .       .  |  .       .  |  0       3  |  .   |   3 
      |\ / \ /|    |   |\   \  |  |   |\     /|  |  |\     /|  |  |\     /   |  |\     /   |   \     /   |   \     /   |   \     /   |     /   /  
3     | d   e |    |   | 1   . |  |   | 1   . |  |  | 1   . |  |  | 4   .    |  | 3   .    |    .   .    |    .   .    |    0   3    |    0   3   
      | |\ /| |    |   |  \  | |  |   |  \  | |  |  |  \    |  |  |  \       |  |  \       |     \       |       /     |    |  /     |    |  /    
2     | | f | |    |   |   1 | |  |   |   1 | |  |  |   1   |  |  |   4      |  |   3      |      4      |      4      |    | 3      |    | 3     
      | |/ \| |    |   |  /  | |  |   |  /  | |  |  |  / \  |  |  |  / \     |  |  / \     |     / \     |     / \     |    |  \     |    |  \    
1     | g   h |    |   | 1   . |  |   | 1   . |  |  | 1   . |  |  | 1   2    |  | .   3    |    0   3    |    0   3    |    0   3    |    0   3   
      |/ \ / \|    |   |  \    |  |   |  \    |  |  |  \    |  |  |  \   \   |  |    / \   |   /   / \   |   /   / \   |   /   / \   |   /   / \  
0     i   j   k    |   0   1   2  |   0   1   2  |  0   1   2  |  0   1   2  |  0   1   2  |  0   1   2  |  0   1   2  |  0   1   2  |  0   1   2 
                                                                                                                                                  
                   |   0.0 - 0.1  |   0.1 - 0.2  |  0.2 - 0.4  |  0.4 - 0.5  |  0.5 - 0.6  |  0.6 - 0.7  |  0.7 - 0.8  |  0.8 - 0.9  |  0.9 - 1.0 

```
and the corresponding set of records:
```
(left, right, parent, children, time)
( 0.5,   1.0,      3,    (1,2), 1.0 )
( 0.4,   0.5,      4,    (1,2), 2.0 )
( 0.6,   0.8,      4,    (0,3), 2.0 )
( 0.0,   0.4,      5,    (0,1), 4.0 )
( 0.4,   0.5,      5,    (0,4), 4.0 )
( 0.5,   0.6,      5,    (0,3), 4.0 )
( 0.0,   0.4,      6,    (2,5), 5.0 )
( 0.8,   1.0,      6,    (0,3), 5.0 )
```


Alternatively, here is the full labeling (using letters instead of numbers to reduce confusion):
```
t                  |              |              |             |             |             |             |             |             |            
                                                                                                                                                  
5       --a--      |     --a--    |     --a--    |    --a--    |    --.--    |    --.--    |    --.--    |    --.--    |    --a--    |    --a--   
       /  |  \     |    /  |  \   |    /     \   |   /     \   |   /     \   |   /     \   |   /     \   |   /     \   |   /     \   |   /  |  \  
4     b   |   c    |   b   |   c  |   b       c  |  b       c  |  b       .  |  b       .  |  .       .  |  .       .  |  b       c  |  .   |   c 
      |\ / \ /|    |   |\   \  |  |   |\     /|  |  |\     /|  |  |\     /   |  |\     /   |   \     /   |   \     /   |   \     /   |     /   /  
3     | d   e |    |   | d   . |  |   | d   . |  |  | d   . |  |  | d   .    |  | d   .    |    .   .    |    .   .    |    d   e    |    d   e   
      | |\ /| |    |   |  \  | |  |   |  \  | |  |  |  \    |  |  |  \       |  |  \       |     \       |       /     |    |  /     |    |  /    
2     | | f | |    |   |   f | |  |   |   f | |  |  |   f   |  |  |   f      |  |   f      |      f      |      f      |   d2 f      |   d2 f     
      | |/ \| |    |   |  /  | |  |   |  /  | |  |  |  / \  |  |  |  / \     |  |  / \     |     / \     |     / \     |    |  \     |    |  \    
1     | g   h |    |  b2 g   . c2 |  b2 g   . c2 | b2 g   . c2 | b2 g   h    | b2 .   h    |    g   h    |    g   h    |    g   h    |    g   h   
      |/ \ / \|    |   |  \    |  |   |  \    |  |  |  \    |  |  |  \   \   |  |    / \   |   /   / \   |   /   / \   |   /   / \   |   /   / \  
0     i   j   k    |   0   1   2  |   0   1   2  |  0   1   2  |  0   1   2  |  0   1   2  |  0   1   2  |  0   1   2  |  0   1   2  |  0   1   2 
                                                                                                                                                  
                   |   0.0 - 0.1  |   0.1 - 0.2  |  0.2 - 0.4  |  0.4 - 0.5  |  0.5 - 0.6  |  0.6 - 0.7  |  0.7 - 0.8  |  0.8 - 0.9  |  0.9 - 1.0 

```
and here are two possible sets of records;
on the left are all one-parent records;
on the right are records produced by merging within one generation overlapping ones:
```
(left, right, parent, children, time)    |  (left, right, parent, children, time)
#  (b,g,0.6)->i and (g,h,0.5)->j
#      and (h,c,0.4)->k        
( 0.0,   0.4,      h,    (k,), 1.0 )     | ( 0.0,   0.4,      h,    (k,), 1.0 )           
( 0.4,   1.0,      c2,   (k,), 1.0 )     | ( 0.4,   1.0,      c2,   (k,), 1.0 )           
( 0.0,   0.5,      g,    (j,), 1.0 )     | ( 0.0,   0.5,      g,    (j,), 1.0 )           
( 0.5,   1.0,      h,    (j,), 1.0 )     | ( 0.5,   1.0,      h,    (j,), 1.0 )           
( 0.0,   0.6,      b2,   (i,), 1.0 )     | ( 0.0,   0.6,      b2,   (i,), 1.0 )           
( 0.6,   1.0,      g,    (i,), 1.0 )     | ( 0.6,   1.0,      g,    (i,), 1.0 )           

#  (f,d,0.8)->g and (e,f,0.2)->h
( 0.0,   0.2,      e,    (h,), 2.0 )     | ( 0.0,   0.2,      e,    (h,), 2.0 )           
( 0.2,   1.0,      f,    (h,), 2.0 )     | ( 0.8,   1.0,      f,    (h,), 2.0 )           
( 0.0,   0.8,      f,    (g,), 2.0 )     | ( 0.2,   0.8,      f,   (g,h), 2.0 )           
( 0.8,   1.0,      d2,   (g,), 2.0 )     | ( 0.0,   0.2,      f,    (g,), 2.0 )           
                                         | ( 0.8,   1.0,      d2,   (g,), 2.0 )        

#  (d,e,0.7)->f                
( 0.0,   0.7,      d,    (f,), 3.0 )     | ( 0.0,   0.7,      d,  (f,d2), 3.0 )           
( 0.7,   1.0,      e,    (f,), 3.0 )     | ( 0.7,   1.0,      e,    (f,), 3.0 )           
( 0.0,   1.0,      d,   (d2,), 3.0 )     | ( 0.7,   1.0,      d,   (d2,), 3.0 )

#  (b,a,0.9)->d and (a,c,0.1)->e
( 0.0,   0.1,      a,    (e,), 4.0 )     | ( 0.0,   0.1,      a,    (e,), 4.0 )           
( 0.1,   1.0,      c,    (e,), 4.0 )     | ( 0.1,   1.0,      c,  (e,c2), 4.0 )           
( 0.0,   0.9,      b,    (d,), 4.0 )     | ( 0.0,   0.9,      b,  (d,b2), 4.0 )           
( 0.9,   1.0,      a,    (d,), 4.0 )     | ( 0.9,   1.0,      a,    (d,), 4.0 )           
( 0.0,   1.0,      b,   (b2,), 4.0 )     | ( 0.9,   1.0,      b,   (b2,), 4.0 ) 
( 0.0,   1.0,      c,   (c2,), 4.0 )     | ( 0.0,   0.1,      c,   (c2,), 4.0 ) 

#  (a,?,1.0)->b and (a,?,1.0)->c
( 0.0,   1.0,      a,    (c,), 5.0 )     | ( 0.0,   1.0,      a,   (b,c), 5.0 )           
( 0.0,   1.0,      a,    (b,), 5.0 )     |
```

How can we merge these to get closer to the smaller set of records above?
In general, we can keep track of who inherits from whom,
splitting and pruning as we go along.
For instance, since `h` has only offspring on `[0,0.4)` and `[0.5,0.8)`;
none on `[0.4,0.5)`,
and no grandchildren on `[0.8,1.0)`,
we could remove the label `h` entirely, combining:
```
                            from      --->   to
 ( 0.0,   0.4,      h,    (k,), 1.0 )  |
 ( 0.5,   1.0,      h,    (j,), 1.0 )  |                                                  
 ( 0.0,   0.2,      e,    (h,), 2.0 )  | ( 0.0,   0.2,      e,    (k,), 2.0 ) 
 ( 0.2,   0.8,      f,   (g,h), 2.0 )  | ( 0.2,   0.4,      f,   (g,k), 2.0 ) 
 ( 0.8,   1.0,      f,    (h,), 2.0 )  | ( 0.4,   0.5,      f,    (g,), 2.0 ) 
                                       | ( 0.5,   0.8,      f,   (g,j), 2.0 ) 
                                       | ( 0.8,   1.0,      f,    (j,), 2.0 ) 
```
We can check that we've still got complete records for the `[0,0.4)` segment of `k`,
and the `[0.5,1)` segment of `j`, which is what we started with.

