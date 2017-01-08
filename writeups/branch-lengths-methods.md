
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

# Y statistic

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


$f_4$ statistic

The statistic $f_4(A_1,A_2,A_3,A_4)$
is the mean path length from MRCA$(a_1,a_3)$ to MRCA$(a_2,a_4)$,
minus the mean path length from MRCA$(a_1,a_4)$ to MRCA$(a_2,a_3)$.
Note that those will often be zeros.
Equivalently (see my arXiv paper), this uses the function
$$\begin{aligned}
    f(x_1,x_2,x_3,x_4)
    = \left( \frac{x_1}{n_1} - \frac{x_2}{n_2} \right)\left( \frac{x_3}{n_3} - \frac{x_4}{n_4} \right)
\end{aligned}$$

