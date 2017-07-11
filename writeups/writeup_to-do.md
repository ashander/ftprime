Notes from July 11 2017 with Jerome, Jaime, Peter

application paper:
- fwd time simulators currently have to drag along irrelevant neutral mutations
- (1) you can just throw these down on the ARG afterwards; (2) to the extent that stats we compute are estimating a quantity about an ARG then using ARGs or trees directly removes the noise of mutation
- definition/spec for nodes and edgesets
- steps to record the tree sequence
- way to attach trees together
- showing how much sequence you can simulate this way vs other ways
- how long to run? 10 N gens or so but you don't need to worry so much if you can combine a thin slice of forwards time simulation on top of a coalescent simulation

tree-sequence methods paper:
- intro: tree sequences are great (msprime) and we want to be able to do more things from them (methods paper). further there are now new methods for estimating tree sequences from real data (Kelleher in prep)
- defining a class of stats that can be computed on trees; alg for computing them
- simplify alg
- arg and tree-sequence correspondence (from the ARGweaver comparison)
