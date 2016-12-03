# Algorithm to output a tree sequence from forwards-time computation:


## Example

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

The algorithm does:

1. Begin with an individual `a` (and another anonymous one) at `t=0`.

```
active:   []
inactive: []
individuals: { a:0.0 }
```

2. `(a,?,1.0)->b` and `(a,?,1.0)->c` at `t=1`

```
#            left right node children  time
active:   [ ( 0.0, 1.0,  a,    (b,c),  0.0 ) ]
inactive: []
individuals: { a:0.0, b:1.0, c:1.0 }
```

3. `(b,a,0.9)->d` and `(a,c,0.1)->e` and then `a` dies at `t=2`

```
# offspring:
#            left right node children  time
active:   [ ( 0.0, 0.1,  a,  (b,c,e),  0.0 ),
            ( 0.1, 0.9,  a,    (b,c),  0.0 ),
            ( 0.9, 1.0,  a,  (b,c,d),  0.0 ),
            ( 0.9, 1.0,  b,     (d,),  1.0 ),
            ( 0.9, 1.0,  c,     (e,),  1.0 ) ]
inactive: []
individuals: { a:0.0, b:1.0, c:1.0, d:2.0, e:2.0 }

# death:
active:   [ ( 0.9, 1.0,  b,     (d,),  1.0 ),
            ( 0.9, 1.0,  c,     (e,),  1.0 ) ]
inactive: [ ( 0.0, 0.1,  a,  (b,c,e),  0.0 ),
            ( 0.1, 0.9,  a,    (b,c),  0.0 ),
            ( 0.9, 1.0,  a,  (b,c,d),  0.0 ) ]
individuals: { a:0.0, b:1.0, c:1.0, d:2.0, e:2.0 }
```

4. `(d,e,0.7)->f` at `t=3`

```
#            left right node children  time
active:   [ ( 0.9, 1.0,  b,     (d,),  1.0 ),
            ( 0.9, 1.0,  c,     (e,),  1.0 ),
            ( 0.0, 0.7,  d,     (f,),  2.0 ),
            ( 0.7, 1.0,  e,     (f,),  2.0 ) ]
inactive: [ ( 0.0, 0.1,  a,  (b,c,e),  0.0 ),
            ( 0.1, 0.9,  a,    (b,c),  0.0 ),
            ( 0.9, 1.0,  a,  (b,c,d),  0.0 ) ]
individuals: { a:0.0, b:1.0, c:1.0, d:2.0, e:2.0, f:3.0 }
```

5. `(f,d,0.8)->g` and `(e,f,0.2)->h` at `t=4`.

```
#            left right node children  time
active:   [ ( 0.9, 1.0,  b,     (d,),  1.0 ),
            ( 0.9, 1.0,  c,     (e,),  1.0 ),
            ( 0.0, 0.7,  d,     (f,),  2.0 ),
            ( 0.7, 1.0,  e,     (f,),  2.0 ),
            ( 0.0, 0.2,  f,     (g,),  3.0 ),
            ( 0.2, 0.8,  f,    (g,h),  3.0 ),
            ( 0.8, 1.0,  d,     (g,),  2.0 ),
            ( 0.0, 0.2,  e,     (h,),  2.0 ),
            ( 0.8, 1.0,  f,     (h,),  3.0 ) ]
inactive: [ ( 0.0, 0.1,  a,  (b,c,e),  0.0 ),
            ( 0.1, 0.9,  a,    (b,c),  0.0 ),
            ( 0.9, 1.0,  a,  (b,c,d),  0.0 ) ]
individuals: { a:0.0, b:1.0, c:1.0, d:2.0, e:2.0, f:3.0, g:4.0, h:4.0 }
```

6. `(b,g,0.6)->i` and `(g,h,0.5)->j` and `(h,c,0.4)->k` at `t=5`.
7. We sample `i` and `j`.


