# General method

Let $A_1, ..., A_n$ be fixed sets of samples.
For a node $z$ in a tree, let $t(z)$ be the length of the branch above $z$, 
and let $X(z) = ( (x_1(z), \bar x_1(z)),  ... (x_n(z), \bar x_n(z)) )$ be a list of **pairs** of integers, 
with $x_k(z)$ the number of samples in $A_k$ that subtend $z$,
and $\bar x_k(z)$ the number of samples in $A_k$ that do *not* subtend $z$.
(This is redundant given the sizes of $A_k$, but allows for easier generalization to mutations.)
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
(mapping True to 1 and False to 0 and discarding the complements $\bar x_1$ and $\bar x_2$).

To compute this, we'll keep track of $(x_1, \ldots, x_n)$ for each node of the current tree, 
and update these as we move along the tree sequence
(from these at the end we get $\bar x_k = |A_k| - x_k$):

 0) initialize $X$
 1) $S \mathrel{+}= u_i * L(T_i)$
 2) move to the next tree and update $X$.

As for updating $X$, note that $X$ is additive: 
if node $z$ has children $y_1, ..., y_m$ then $X(z) = \sum_j X(Y_j)$.
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

Here we describe how to compute a collection of statistics.
All the statistics we define will be functions of a collection of *disjoint* leaf sets $A_1$, \ldots, $A_n$;
and we will require them to have the following properties:

- The expected value of a statistic $T(A_1, \ldots, A_n)$
    should not depend on the sizes of the leaf sets,
    as long as all members of each leaf set are exchangeable.

This is natural if you think about, say, $A_1$ as a sample of individuals from a single "population",
so that larger samples should only give us a better estimate of whatever $T$ is estimating.

This probably means they are all unbiased estimators of some (single) descriptive quantity
in some sense.

It would also be natural to allow overlap between leaf sets and require statistics to be insensitive to degree of overlap,
but this becomes onerous (see examples below).
We can allow overlap naturally between leaf sets by setting two leaf sets equal,
and correcting for this fact;
the relationship of $f_2$ and $f_3$ to the $f_4$ is an example of this.

## Mean TMRCA

Mean time to most recent common ancestor between two individuals is the mean distance in the tree between them, 
averaged across sites, divided by two.
Between two groups, it also averages over choices of individuals from those groups. 
We assume that the sets are disjoint, but for computing "self" comparisons -- mean TMRCA of a leaf set to itself --
we correct to exclude self comparisons.

To compute mean divergence between two groups, $A_1$ and $A_2$, of sizes $n_1$ and $n_2$,
the weight we assign to a branch is equal to the number of paths from an element of $A_1$ to an element of $A_2$
that pass through that branch
divided by the total number of paths, $n_1 n_2$.
If below a branch there are $x_1$ from $A_1$ and $x_2$ from $A_2$,
and as always $\bar x_k = n_k - x_k$,
then the number of such paths are $x_1 \bar x_2 + \bar x_1 x_2$,
so that we should use
$$f((x_1, \bar x_1), (x_2, \bar x_2)) = \frac{x_1 \bar x_2 + \bar x_1 x_2}{n_1 n_2} .$$

**For self comparisons** when $A_1 = A_2$,
we want the mean pairwise TMRCA between two *distinct* samples from the set.
Not accounting for this makes the statistic depend on sample size,
as in smaller samples, self comparisons account for a larger fraction.
$$\begin{aligned}
    f((x_1,\bar x_1)) = \frac{2 x_1 \bar x_1}{n_1 (n_1 - 1)} .
\end{aligned}$$

## Y statistics

Take a set of three distinct samples $a$, $b$, $c$.
By $Y(a;b,c)$ we mean the mean length of any branches from which either ($a$), or (both of $b$ and $c$),
but not $a$ and any of $b$ and $c$, inherit.
With $A_1=\{a\}$ and $A_2 = \{b\}$ and $A_3 = \{c\}$, this therefore corresponds to the function
$$\begin{aligned}
    f((x_1,\bar x_1), (x_2, \bar x_2), (x_3, \bar x_3))
    &=
    ( x_1=1 \;\text{and}\; x_2=x_3=0 ) \;\text{or}\; ( x_1=0 \;\text{and}\; x_2=x_3=1 )
    &= 
    x_1 \bar x_2 \bar x_3
    + \bar x_1 x_2 x_3
\end{aligned}$$

Extending this to groups we would average over choices of $a$, $b$, and $c$
from those groups, to obtain $Y_3(A;B,C)$.
This at first seems simliar to $f_3$, but the $f_3$ statistic involves averaging over choices of *four distinct* individuals,
and this statistic uses only three.
The weighting function for the $Y_3$ is then
$$\begin{aligned}
    f_{Y3}((x_1,\bar x_1),(x_2,\bar x_2),(x_3,\bar x_3)) 
    &= 
        \frac{ x_1 \bar x_2 \bar x_3 + \bar x_1 x_2 x_3 }{ n_1 n_2 n_3 } 
\end{aligned}$$

We can also extend this: $Y_2(A;B) = Y_3(A;B,B)$ and $Y_1(A) = Y_3(A;A,A)$ but corrected for sample sizes as above,
we have that
$$\begin{aligned}
    f_{Y2}((x_1,\bar x_1),(x_2,\bar x_2))
    &= 
        \frac{ x_1 \bar x_2 (\bar x_2 - 1) + \bar x_1 x_2 (x_2-1) }{ n_1 n_2 (n_2-1) } 
\end{aligned}$$
and
$$\begin{aligned}
    f_{Y1}((x_1,\bar x_1),(x_2,\bar x_2),(x_3,\bar x_3)) 
    &= 
        \frac{ x_1 \bar x_1 (\bar x_1-1) + \bar x_1 x_1 (x_1-1) }{ n_1 (n_1-1) (n_1-1) } .
\end{aligned}$$


## $f_4$ statistic

The statistic $f_4(A_1,A_2,A_3,A_4)$
is the mean path length from MRCA$(a_1,a_3)$ to MRCA$(a_2,a_4)$,
minus the mean path length from MRCA$(a_1,a_4)$ to MRCA$(a_2,a_3)$.
Note that those will often be zeros.
To rewrite this,
note that if $Z_1$ is a randomly chosen allele from $A_1$,
then $\E[Z_1]=p_1=x_1/n_1$, the allele frequency in $A_1$;
if these are independent then $(p_1-p_2)(p_3-p_4) = \E[(Z_1-Z_2)(Z_3-Z_4)]$.
This makes it easy to rewrite the weighting function 
(useful later in correcting for overlapping samples):
$$\begin{aligned}
    f_4((x_1,\bar x_1),(x_2,\bar x_2),(x_3,\bar x_3),(x_4,\bar x_4))
    &= 
        \left( \frac{x_1}{n_1} - \frac{x_2}{n_2} \right)\left( \frac{x_3}{n_3} - \frac{x_4}{n_4} \right) \\
    &= 
        \frac{ ( x_1 n_2 - x_2 n_1 ) ( x_3 n_4 - x_4 n_3 ) }{ n_1 n_2 n_3 n_4 } \\
    &= 
        \frac{ ( x_1 \bar x_2 - x_2 \bar x_1 ) ( x_3 \bar x_4 - x_4 \bar x_3 ) }{ n_1 n_2 n_3 n_4 } \\
    &= 
        \frac{ x_1 \bar x_2 x_3 \bar x_4 + \bar x_1 x_2 \bar x_3 x_4  - x_1 \bar x_2 \bar x_3 x_4 - \bar x_1 x_2 x_3 \bar x_4  }{ n_1 n_2 n_3 n_4 } 
\end{aligned}$$

<!--
```r
# (a-b)(c-d) = ac + bd - ad - bc = abcd + a(1-b)c(1-d) + abcd + (1-a)b(1-c)d - abcd - a(1-b)(1-c)d - abcd (1-a)bc(1-d)
h <- function (x,n) { (x[1]/n[1] - x[2]/n[2]) * (x[3]/n[3] - x[4]/n[4]) }
f <- function (x,n) { ( ( x[1] * n[2] - x[2] * n[1] ) * ( x[3] * n[4] - x[4] * n[3] ) )/( n[1] * n[2] * n[3] * n[4] ) }
g <- function (x,n) { 
   (  x[1] * (n[2] - x[2]) * x[3] * (n[4] - x[4]) 
    + (n[1] - x[1]) * x[2] * (n[3] - x[3]) * x[4]  
    - x[1] * (n[2] - x[2]) * (n[3] - x[3]) * x[4] 
    - (n[1] - x[1]) * x[2] * x[3] * (n[4] - x[4]) 
   )/( n[1] * n[2] * n[3] * n[4] ) }
```
-->

## $f_3$ statistic

The statistic $f_3(A_1;A_2,A_3)$ is just $f_4(A_1,A_2;A_1,A_3)$ corrected for overlapping samples.
The function we use is the probability that the branch
separates $(a,b)$ from $(c,d)$, where $a$ and $b$ are *distinct* draws from $A_1$
and $c$ and $d$ are from $A_2$ and $A_3$,
minus the probability that $(a,d)$ are distinct from $(b,c)$.
This uses the weighting function
$$\begin{aligned}
    f_3((x_1,\bar x_1),(x_2,\bar x_2),(x_3,\bar x_3))
    &= \frac{ x_1 (x_1-1) (n_2 - x_2) (n_3 - x_3) 
            + (n_1 - x_1) (n_1 - x_1 - 1) x_2 x_3 
            - x_1 (n_1 - x_1) (n_2 - x_2) x_3
            - (n_1 - x_1) x_1 x_2 (n_3 - x_3)
        }{ n_1 (n_1-1) n_2 n_3 } \\
    &= \frac{n_1}{n_1-1} \left( \frac{x_1}{n_1} - \frac{x_2}{n_2} \right)\left( \frac{x_1}{n_1} - \frac{x_3}{n_3} \right)
        - \frac{1}{n_1-1} \left\{ 
            \frac{x_1}{n_1} 
            \frac{\bar x_2}{n_2} 
            \frac{\bar x_3}{n_3} 
            + 
            \frac{\bar x_1}{n_1}
            \frac{x_2}{n_2}
            \frac{x_3}{n_3}
        \right\} 
\end{aligned}$$

## $f_2$ statistic

The statistic $f_2(A_1;A_2)$ is just $f_4(A_1,A_2;A_1,A_2)$ corrected for overlapping samples.
The function we use is the probability that the branch
separates $(a,b)$ from $(c,d)$, where $a$ and $b$ are *distinct* draws from $A_1$
and $c$ and $d$ are *distinct* draws from $A_2$,
minus the probability that $(a,d)$ are distinct from $(b,c)$.
This uses the weighting function
$$\begin{aligned}
    f_2((x_1,\bar x_1),(x_2,\bar x_2))
    &= \frac{ x_1 (x_1-1) (n_2 - x_2) (n_2 - x_2 - 1) 
            + (n_1 - x_1) (n_1 - x_1 - 1) x_2 (x_2 - 1) 
            - x_1 (n_1 - x_1) (n_2 - x_2) x_2
            - (n_1 - x_1) x_1 x_2 (n_2 - x_2)
        }{ n_1 (n_1-1) n_2 (n_2-1) } 
    &= \frac{ x_1 (x_1-1) \bar x_2 (\bar x_2 - 1)
            + \bar x_1 (\bar x_1 - 1) x_2 (x_2 - 1) 
            - 2 x_1 \bar x_1 \bar x_2 x_2
        }{ n_1 (n_1-1) n_2 (n_2-1) } 
\end{aligned}$$



## Covariance

If the genotype of the $x^\text{th}$ individual at site $i$ is $G_{xi}$,
then the genetic covariance is
$$\begin{aligned}
    \Sigma_{xy} 
        &= \frac{1}{n^2} \sum_{uv} (G_{xi} - G_{ui})(G_{yi} - G_{vi}) \\
        &= \left(G_{xi} - \frac{1}{n}\sum_v G_{ui}\right)\left(G_{yi} - \frac{1}{n}\sum_v G_{vi}\right) 
\end{aligned}$$

Note that the function $(a-b)(c-d)$ 
takes the value $+1$ if $a=c \neq b=d$ 
and takes the value $-1$ if $a=d \neq b=c$,
and is zero otherwise.
Therefore, the function we need to implement takes the value $z$, where

- if $G_x=G_y$, with $f$ the frequency of the *other* allele, $z=f^2$.
- if $G_x \neq G_y$, with $f$ the frequency of one allele (doesn't matter which), $z=-f(1-f)=f^2-f=(1-f)^2-(1-f)$.

If we let $A_1=\{x\}$, $A_2=\{y\}$, and $A_3=X$ (all the samples, including $x$ and $y$),
then this is
$$\begin{aligned}
    f((x,\bar x),(y,\bar y),(z, \bar z))
    =
    -\frac{z}{n}\frac{\bar z}{n} +
    \begin{cases}
        z/n \qquad & \text{if}\; x=y=0 \\
        \bar z/n \qquad & \text{if}\; x=y=1  .
    \end{cases}
\end{aligned}$$

Note that theory tells us that the function of branch lengths that this estimates is actually different.

## Mean TMRCA with overlapping leaf sets

**However,** the above does not account for a case when $A_1$ and $A_2$ overlap.
Indeed, it may be desirable to pass in $A_1 = A_2$ to obtain the mean pairwise TMRCA between two *distinct* samples from the set.
Not accounting for this makes the statistic depend on sample size,
as in smaller samples, self comparisons account for a larger fraction.
Since if the two elements are equal divergence is zero,
we need only divide by the probability that the two chosen samples are distinct,
resulting in
$$\begin{aligned}
    f(x_1,x_2) = \frac{x_1 (n_2-x_2) + (n_1-x_1) x_2}{n_1 n_2 - n_{1 \cap 2}^2} ,
\end{aligned}$$
where $n_{1 \cap 2}$ is the number of samples in both $A_1$ and $A_2$.

# Considerations with overlapping leaf sets

Here is discussion of how to compute the statistics in a consistent way (as defined above)
when leaf sets are not disjoint.  We don't use this.

## Y statistic with overlapping leaf sets

To compute a $Y$ statistic with overlapping leaf sets, 
we want $f(x_A;x_B,x_C)$ to be equal to the probability that
the chosen edge separates $a$ from $b$ and $c$, where $a$, $b$, and $c$ are uniform random samples from $A$, $B$, and $C$ respectively,
but *chosen to be all different*, i.e., conditioned on $a \neq b \neq c \neq a$.

Equivalently: suppose we have three boxes, labeled $A$, $B$, and $C$;
in box $A$ there are $n_A$ balls, $x_A$ of which are colored red and the rest black, and on which are written some numbers; 
same for $B$ and $C$.
We imagine those balls in different boxes with the same number are the *same*, in particular, they have the same color.
Suppose we pick one ball from each box.
We want the probability that the ball from $A$ differs in color from the ones from $B$ and $C$,
conditioned on the three balls all having different numbers.
Write $p_A = x_A/n_A$ for the proportion of balls in box $A$ that are red,
and $p_{B \cap C} = x_{B \cap C}/n_{B \cap C}$
likewise for the proportion of balls in *both* $A$ and $B$ that are red.
The probability that the three chosen balls are of the required colors without regards identiy is
$$\begin{aligned}
    p_A (1-p_B) (1-p_C) + (1-p_A) p_B p_C .
\end{aligned}$$
However, this has included the case where $b=c$:
the probability that $b=c$ and the balls are of the required color is
$$\begin{aligned}
    p_A (1-p_{B \cap C})  + (1-p_A) p_{B \cap C} .
\end{aligned}$$
Overall, the probability that we get three distinct numbers is
$$\begin{aligned}
    p_d = 1 
    - \frac{n_{A \cap B}^2}{n_A n_B}
    - \frac{n_{A \cap C}^2}{n_A n_C}
    - \frac{n_{B \cap C}^2}{n_B n_C}
    + 2 \frac{n_{A \cap B \cap C}^3}{n_A n_B n_C} .
\end{aligned}$$
Putting this together, we get that
$$\begin{aligned}
    f_Y(x_1;x_2,x_3)
    &=
    \frac{1}{p_d} \left( p_A (1-p_B) (1-p_C) + (1-p_A) p_B p_C 
            - p_A (1-p_{B \cap C})  + (1-p_A) p_{B \cap C} \right) \\
    &=
    \frac{ x_A (n_B - x_B) (n_C - x_C) + (n_A - x_A) x_B x_C - x_A (n_{B \cap C} - x_{B \cap C}) + (1 - x_A)( n_{B \cap C} - x_{B \cap C}) }
        { n_A n_B n_C - n_{A \cap B}^2 n_C - n_{A \cap C}^2 n_B - n_A n_{B \cap C}^2 + 2 n_{A \cap B \cap C}^3} .
\end{aligned}$$

## $f_4$ statistic, with overlapping samples

The statistic $f_4(A,B;C,D)$ uses the function $(p_A - p_C)(p_B - p_D)$,
corrected for sample size.
This adds branches that separate $(a,b)$ from $(c,d)$
and subtracts branches that separate $(a,d)$ from $(b,c)$,
where $a$, $b$, $c$, and $d$ are *distinct* samples from $A$, $B$, $C$ and $D$ respectively.
Chasing through the inclusion-exclusions,
considering the cases 

- a b | c d
- a d | b c
- a=b | c d
- a b | c=d
- a=b | c=d
- a=d | b c
- a d | b=c
- a=d | b=c

we get a numerator for the $f_4$ function of 
$$\begin{aligned}
    u_4(x_A,x_B;x_C,x_D)
    &=
    ( x_A x_B (n_C - x_C) (n_D - x_D) + (n_A - x_A) (n_B - x_B) x_C x_D )                                  
    \\ &\qquad {} 
    - ( x_A x_C (n_B - x_B) (n_D - x_D) + (n_A - x_A) (n_C - x_C) x_B x_D )                    
    \\ &\qquad {} 
    - ( x_{A \cap B} (n_C - x_C) (n_D - x_D) + (n_{A \cap B} - x_{A \cap B}) x_C x_D )                                    
    \\ &\qquad {} 
    - ( x_A x_B (n_{C \cap D} - x_{C \cap D})  + (n_A - x_A) (n_B - x_B) x_{C \cap D} ) 
    \\ &\qquad {} 
    + ( x_{A \cap B} (n_{C \cap D} - x_{C \cap D}) + (n_{A \cap B} - x_{A \cap B}) x_{C \cap D} )                         
    \\ &\qquad {} 
    + ( x_{A \cap D} (n_B - x_B) (n_C - x_C) + (n_{A \cap D} - x_{A \cap D}) x_B x_C )
    \\ &\qquad {} 
    + ( x_A x_D (n_{B \cap C} - x_{B \cap C})  + (n_A - x_A) (n_D - x_D) x_{B \cap C} )          
    \\ &\qquad {} 
    - ( x_{A \cap D} (n_{B \cap C} - x_{B \cap C}) + (n_{A \cap D} - x_{A \cap D}) x_{B \cap C} )
\end{aligned}$$
The denominator has $2^6$ terms.


## In general with overlapping leaf sets

The general procedure for correcting a statistic for sample size
goes something like this.
Suppose that our statistic counts any branch that separates
a collection of samples, one each from $A_1$, \ldots, $A_k$,
from another collection of samples, one each from $B_1$, \ldots, $B_m$.
For a given subset $I$ of $\{1, \ldots, k\}$,
let $p_I = x_{\cap_I A_i} / n_{\cap_I A_i}$ be the proportion of samples below the branch 
among all those in the intersection of $A_i$ across $i \in I$;
and likewise let $q_J = x_{\cap_J B_j} / n_{\cap_J B_j}$.
The statistic without correction will be calculated with
$\prod_i p_i \prod_j (1-q_j) + \prod_i (1-p_i) \prod_j q_j$;
to correct it we can do inclusion-exclusion to get both a numerator and a denominator.
The numerator sums over partitions of $\{1,\ldots,k\}$ into $I_1, \ldots, I_{\ell}$
and of $\{1,\ldots,m\}$ into $J_1, \ldots, J_{m}$
a constant multiplied by 
$\prod_i p_{I_i} \prod_j (1-q_{J_j}) + \prod_i (1-p_{I_i}) \prod_j q_{J_j}$.


# Multiple mutations

**Observation:** branches are equivalent to splits (i.e., bipartitions) are equivalent to biallelic sites.

Currently, the statistic is defined by

-  a list of sets of samples $A_1,...,A_n$, and
-  a function $f()$ that takes a list of tuples and returns a number;
-  then each mutation (or unit of branchlength) counts towards the statistic weighted by $f()$ of the vector of numbers of individuals in each set of samples inheriting from that mutation (or branch), $x_1,...,x_n$.

This definition applies perfectly fine to possibly recurrent mutations - the only issue is that we have to do a bit more work to determine what the $x$'s are for a particular mutation, since if it has occurred more than once then these aren't the same as the $x$'s for the branch it lies on.

What about sites with more than one allele? There is not currently a consensus about how to compute single-site statistics like divergence using e.g. triallelic sites, and there's no obvious (simple) best answer, so this definition is as good as any. Concretely, let's take the definition above after replacing "each mutation" by "each allele". So, at sites with $k$ alleles, we'd compute the vector of $x$'s for each of the $k$ alleles, apply $f()$ to each of these, and sum them.

This is slightly different than before, because with biallelic sites if we know the $x$'s for a particular mutation, we also know them for the alternate allele, and we exploited this in defining particular $f()$s. We can't just sum over derived mutations, either: define divergence to be "density of differing sites"; let $A_1$ and $A_2$ both contain exactly one sample, and let $f(x,y) = xor(x,y)$. This counts derived alleles inherited by $A_1$ or $A_2$ but not both. Biallelic sites are counted correctly, but triallelic sites where $A_1$ and $A_2$ each inherit different derived alleles would get counted twice, incorrectly. The solution is to let $f(x,y) = xor(x,y)/2$, and sum over all alleles, so that

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
