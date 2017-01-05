
Let $A_1, ..., A_n$ be fixed sets of samples.
For a node $z$ in a tree, let $t(z)$ be the length of the branch above $z$, 
and let $X(z) = ( x_1(z), ... x_n(z) )$ be a list of integers, 
with $x_k(z)$ the number of samples in $A_k$ that subtend $z$.  
We will want to count the total length of all branches with $X$ satisfying a certain condition, 
which we encode by a function $f : X \to {True, False}$.  
Then for a given tree $T$ we define $L(T)$ to be this total length:
$$
   L(T) = \sum_{z : f(X(z))} t(z)
$$

Now, let $L$ be the sum of these along the tree sequence: if $T_i$ extends for length $u_i$ along the genome, then we want
$$
   S = \sum_i u_i L(T_i) .
$$

For instance: we get sequence length multiplied by average TMRCA between $a$ and $b$ with $A_1 = {a}$ and $A_2 = {b}$ and $f(x_1,x_2) = \text{xor}(x_1,x_2)$.

To compute this, we'll keep track of X(z) for each node of the current tree, and update these as we move along the tree sequence:

 0) initialize $X$
 1) $S += u_i * L(T_i)$
 2) move to the next tree and update $X.$

As for updating $X$, note that $X$ is additive: 
if node $z$ has children $y_1, ..., y_m$ then $X(z) = sum_j X(Y_j)$.  
So, to update we just need to propagate changes up from any node that has changed.
To avoid re-computing $f(X)$ we would also cache these and only update those that were changed: 
to do this, let $F(z) = f(X(z))$. Concretely, we can do this by:

  - For each record with parent $z$ and children $y_1, ..., y_m$ that is removed,

      - if F(X(z)) then L -= t(z)
      - if F(X(y_i)) then L -= t(y_i) for each $i$
      - $X(u) -= X(z)$ for each ancestor $u$ of $z$ up the tree, updating $L$ if $F(X(u))$ changes

  - then, for each record with parent $z$ and children $y_1, ... y_m$ that is added,

      - if $F(X(y_i))$ then $L += t(y_i)$ for each $i$ 
      - $X(z) = X(y_1) + ... + X(y_m)$
      - if $F(X(z))$ then $L += t(z)$
      - $X(u) += X(z)$ for each ancestor $u$ of $z$ up the tree, updating $L$ if $F(X(u))$ changes

After this we have $L(T)$ for the new tree, so we can do $S += L * u_i$.

Note that I've said this in terms of branch lengths, 
but by swapping out branch lengths for numbers of mutations, 
we'd be able to easily compute a large class of sequence statistics 
(e.g., Patterson's $f$ statistics, which involve four groups $A_1, ... A_4$).  
