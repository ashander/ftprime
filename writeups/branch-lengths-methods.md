# General method

Let $A_1, ..., A_n$ be fixed sets of samples.
For a node $z$ in a tree, let $t(z)$ be the length of the branch above $z$, 
and let $X(z) = ( x_1(z), ... x_n(z) )$ be a list of integers, 
with $x_k(z)$ the number of samples in $A_k$ that subtend $z$.  
We will want to count the total length of all branches weighted by a function $f(X)$.
Then for a given tree $T$ we define $L(T)$ to be this total length:
$$
   L(T) = \sum_{z} f(X(z)) t(z)
$$

Now, let $L$ be the sum of these along the tree sequence: if $T_i$ extends for length $u_i$ along the genome, then we want
$$
   S = \sum_i u_i L(T_i) .
$$

For instance: we get sequence length multiplied by average TMRCA between $a$ and $b$ with $A_1 = {a}$ and $A_2 = {b}$ and $f(x_1,x_2) = \text{xor}(x_1,x_2)$
(mapping True to 1 and False to 0).

To compute this, we'll keep track of X(z) for each node of the current tree, and update these as we move along the tree sequence:

 0) initialize $X$
 1) $S += u_i * L(T_i)$
 2) move to the next tree and update $X.$

As for updating $X$, note that $X$ is additive: 
if node $z$ has children $y_1, ..., y_m$ then $X(z) = sum_j X(Y_j)$.
So, to update we just need to propagate changes up from any node that has changed.
Concretely, we can do this by:

  - For each record with parent $z$ and children $y_1, ..., y_m$ that is removed,

      - $L -= t(y_i) f(X(y_i))$ for each $i$ because these branches no longer exist
      - $f_\text{old} = f(X(z))$
      - $dx = X(y_1)+\cdots+X(y_m)$
      - $X(z) -= dx$ (note $X(z)$ is not 0 if $z\in A_i$ for some $i$)
      - $L += t(z) (f(X(z))-f_\text{old})$
      - for each ancestor $u$ of $z$ up the tree,
          do $f_\text{old}(u)=f(X(u))$ and $X(u) -= dx$ and $L += t(u) (f(X(u))-f_\text{old}(u))$ 

  - then, for each record with parent $z$ and children $y_1, ... y_m$ that is added,

      - $L += t(y_i) f(X(y_i))$ for each $i$ 
      - $f_\text{old} = f(X(z))$
      - $dx = X(y_1) + ... + X(y_m)$
      - $X(z) += dx$
      - $L += t(z) (f(X(z))-f_\text{old})$
      - for each ancestor $u$ of $z$ up the tree,
          do $f_\text{old}(u)=f(X(u))$ and $X(u) += dx$ and $L += t(u) (f(X(u))-f_\text{old}(u))$ 

After this we have $L(T)$ for the new tree, so we can do $S += L * u_i$.

Note that I've said this in terms of branch lengths, 
but by swapping out branch lengths for numbers of mutations, 
we'd be able to easily compute a large class of sequence statistics 
(e.g., Patterson's $f$ statistics, which involve four groups $A_1, ... A_4$).  

# Statistics

## Divergence

To compute mean divergence between two groups, $A_1$ and $A_2$, of sizes $n_1$ and $n_2$,
the weight we assign to a branch is equal to the number of paths from an element of $A_1$ to an element of $A_2$
that pass through that branch
divided by the total number of paths, $n_1 n_2$.
If below a branch there are $x_1$ from $A_1$ and $x_2$ from $A_2$
then the number of such paths are $x_1 (n_2-x_2) + (n_1-x_1) x_2$,
so we want
$$\begin{aligned}
    f(x_1,x_2) = x_1 (n_2-x_2) + (n_1-x_1) x_2 .
\end{aligned}$$

## Y statistic

By $Y(a;b,c)$ we mean the mean length of any branches from which either ($a$), or (both of $b$ and $c$),
but not $a$ and any of $b$ and $c$, inherit.
With $A_1=\{a\}$ and $A_2 = \{b,c\}$, this therefore corresponds to the function
$$\begin{aligned}
    f(x_1,x_2) = (( (x_1==1) \text{and} (x_2==0) ) \text{or} ( (x_1==0) \text{and} (x_2==2) ))
\end{aligned}$$
The generalization gives the length of the internal branches 
that separate $A_1$ from $A_2$, if any, using
$$\begin{aligned}
    f(x_1,x_2) = (( (x_1==n_1) \text{and} (x_2==0) ) \text{or} ( (x_1==0) \text{and} (x_2==n_2) ))
\end{aligned}$$


## $f_4$ statistic

The statistic $f_4(A_1,A_2,A_3,A_4)$
is the mean path length from MRCA$(a_1,a_3)$ to MRCA$(a_2,a_4)$,
minus the mean path length from MRCA$(a_1,a_4)$ to MRCA$(a_2,a_3)$.
Note that those will often be zeros.
Equivalently (see my arXiv paper), this uses the function
$$\begin{aligned}
    f(x_1,x_2,x_3,x_4)
    = \left( \frac{x_1}{n_1} - \frac{x_2}{n_2} \right)\left( \frac{x_3}{n_3} - \frac{x_4}{n_4} \right)
\end{aligned}$$


# Multiple mutations

Observation: branches are equivalent to splits (i.e., bipartitions) are equivalent to biallelic sites.

Currently, the statistic is defined by

-  a list of sets of samples $A_1,...,A_n$, and
-  a function $f()$ that takes a list of integers and returns a number;
-  then each mutation (or unit of branchlength) counts towards the statistic weighted by $f()$ of the vector of numbers of individuals in each set of samples inheriting from that mutation (or branch), $N_1,...,N_n$.

This definition applies perfectly fine to possibly recurrent mutations - the only issue is that we have to do a bit more work to determine what the $N$'s are for a particular mutation, since if it has occurred more than once then these aren't the same as the $N$'s for the branch it lies on.

What about sites with more than one allele? There is not currently a consensus about how to compute single-site statistics like divergence using e.g. triallelic sites, and there's no obvious (simple) best answer, so this definition is as good as any. Concretely, let's take the definition above after replacing "each mutation" by "each allele". So, at sites with $k$ alleles, we'd compute the vector of $N$'s for each of the $k$ alleles, apply $f()$ to each of these, and sum them.

This is slightly different than before, because with biallelic sites if we know the $N$'s for a particular mutation, we also know them for the alternate allele, and we exploited this in defining particular $f()$s. We can't just sum over derived mutations, either: define divergence to be "density of differing sites"; let $A_1$ and $A_2$ both contain exactly one sample, and let $f(x,y) = xor(x,y)$. This counts derived alleles inherited by $A_1$ or $A_2$ but not both. Biallelic sites are counted correctly, but triallelic sites where $A_1$ and $A_2$ each inherit different derived alleles would get counted twice, incorrectly. The solution is to let $f(x,y) = xor(x,y)/2$, and sum over all alleles, so that

-  at a site with ancestral state G, sampled sequences G, C, we'd get $N(G) = (1,0)$ and $N(C) = (0,1)$ so that $f(N(G)) + f(N(C)) = 1$
-  at a site with ancestral state G, sampled sequences T, C, we'd get $N(T) = (1,0)$ and $N(C) = (0,1)$ and $N(G) = (0,0)$, so that again $f(N(G)) + f(N(C)) + f(N(T))= 1$.

Any statistic that doesn't require knowing the ancestral state should be invariant under changing which allele is 'ancestral', so I think (but haven't checked) that translating previous $f()$'s for particular statistics to this scheme only requires dividing them all by two (as in the divergence example).

So, my proposal is to replace the definition of the statistic so that for branch lengths both $f(N)$ and $f(|A|-N)$ contribute. This then requires us to work out the $N$s for each allele at each polymorphic site, maybe by

-  start with the $N$s computed for each branch in the current algorithm
-  postprocess these to obtain the $N$ for alleles at that site (need to work out how to do this)


We have two proposals for which direction to move in.
Pros/cons:

- roots-to-leaves allows only comparison of subsequent mutations for inheritance,
    if mutations are in preorder
- leaf-to-roots does not change an individual's genotype once it is found; 
    roots-to-leaves must traverse the entire tree to find an individual's genotype


## Algorithm 1: from leaves to root

Suppose that for a given site we have

a. a list of mutations $(i_k,x_k)$: node and derived state, ordered by increasing time-ago of the node,
b. the sparse tree for this site
c. for each node $i$, the set of samples $D(i)$ that descend from the node

I'll assume there is at most one mutation per node; for things we're considering here only the most recent one matters.

Here is an algorithm to obtian genotypes at this site. For each $k$ beginning from $k=1$:

1. Assign everyone in $D(i_k)$ the allele $x_k$.
2. Move back up through the tree from $i_k$, checking to see if we pass through any other $i_j$ along the way (easy as they are ordered by time); 
    for any $j$ that $i_k$ inherits from $i_j$, update $D(i_j) \mapsto D(i_j) \setdiff D(i_k)$.
3. At the end, assign anyone without an allele to the ancestral state.

As for computing statistics: the key is producing, for a given subset $A$ of samples and for each allele $x$ at the site, 
how many individuals in $A$ have allele $x$. Now instead of (c) above suppose we know
(c') for each node $i$, the number of individuals in $A$ inherited from that node, denoted $N(i)$.
To produce the required numbers, say, $N(x)$, the algorithm above can be used, but updating $N(i_j) \mapsto N(i_j) - N(i_k)$ instead. 
We'll also need to check if the allele $x_k=x_j$ 
for some $j \lt k$, in which case we update $N(i_j) += N(i_k)$.


## Algorithm 2: From root to leaves

It's interesting that your algorithm for generating genotypes goes from leaf-to-root; I would have thought the other way around was the way to do it. Consider the following tree:

[!example_tree.png]

Suppose we have one site with ancestral state A, and three mutations: [(19, T), (13, C), (14, G)]. We want to compute the frequency of all alleles.

1. First, set f[T] = 6, because there are 6 leaves below 19.
2. Then, set f[C] = 2 (two leaves below 13) and f[T] = 6 - 2 = 4, because 19 is an ancestor of 13.
3. Then, set f[G] = 2.
4. Finally, set f[A] = n - sum(f).

I've ordered the nodes at the site in pre-order here, as I think it makes it easier to reason about what's happening in the current subtree (i.e., here we fully finish up with the subtree below 19 before moving to look at the mutation at 14). It seems to me that there is a bit of efficiency to be gained there because we don't have to look at the entire list of nodes but only the last when we're moving back up the tree to see if we're affecting any existing mutations.

There is a slight trickiness when computing the frequency of the alleles when we have back mutations to the ancestral state, but it's a small detail. I think this will work for calculating frequencies within subsets as well. The cost should be O((k - 1) * height of tree) if we have k mutations at a site, since in worst case we traverse up to the root k - 1 times to find ancestral mutations. This probably about as efficient as it's going to get.


## Comparison

A possible advantage of the leaf-to-root direction is that once you find a genotype, it doesn't change.  If trees might have many mutations this could be good, as we wouldn't have to actually reach the root of the tree (only up to the most recent mutation for each sample of interest).

Your point about the number of comparisons between mutations is probably better, though.  However, I think you'd have to go a bit further than just the last node.  For instance, suppose that the mutations on that tree were  [(20, T), (16,C), (13, C), (14, G)].  We'd need to keep track of the current list of mutations on the way back to the root, and compare all the way up that.  Call this list S.

1.  First, set f[T] = 8, because there are 8 leaves below 20; set S=[20].
2.  Then, since 20 is an ancestor of 16 set f[C] = 2 (two leaves below 16) and f[T] = 8 - 2 = 6; and S=[20,16].
3.  Next, check that 16 is not an ancestor of 13, but 20 is, so set f[C] = 2 + 2 = 4 (two more leaves below 16) and f[T] = 6 - 2 = 4; and S=[20,13].
4.  Then, since 14 is not descended from 20 or 14, set f[G] = 2.
5.  Finally, set f[A] = n - sum(f).
