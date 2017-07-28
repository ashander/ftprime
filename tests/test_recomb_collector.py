import pytest
import simuPOP as sim
import math
from collections import Counter
import random

import ftprime
import msprime


# increases reproducibility by
# reset IDs and setting seed for each test
# ('autouse' makes every test in this file run this)
@pytest.fixture(scope='function', autouse=True)
def reset_id_tagger():
    print('Setup simupop seed')
    sim.setOptions(seed=111)


def check_tables(args):
    nodes = args.nodes
    assert(nodes.num_rows == args.num_nodes)
    edgesets = args.edgesets
    # check edgesets are in order and all parents are recorded
    node_times = nodes.time
    last_time = 0.0
    for p in edgesets.parent:
        assert(node_times[p] >= last_time)
        last_time = node_times[p]
        assert(p < args.num_nodes)
    for ch in edgesets.children:
        assert(ch < args.num_nodes)


@pytest.fixture(scope="function", params=[
    lambda recombinator, popsize, id_tagger: sim.RandomMating(
        ops=[
            id_tagger,
            recombinator
        ]),
    # Overlapping generations mating system -- popsize grows by 2x
    lambda recombinator, popsize, id_tagger: sim.HeteroMating([sim.RandomMating(
        ops=[
            id_tagger,
            recombinator
        ]),
        sim.CloneMating()],
        subPopSize=popsize * 2),
    # Overlapping generations mating system -- popsize grows by 5x
    lambda recombinator, popsize, id_tagger: sim.HeteroMating([sim.RandomMating(
        ops=[
            id_tagger,
            recombinator
        ]),
        sim.CloneMating()],
        subPopSize=popsize * 5),
    ])
def make_pop(request):
    # request.param stores a lambda function to make mating scheme
    # each test that uses this fixture will be run for both entries in 'params'
    mating_scheme_factory = request.param

    def _make_pop(popsize, nloci, locus_position, id_tagger, init_geno,
                  recomb_rate, generations, length):
        random.seed(123)
        pop = sim.Population(
                size=[popsize],
                loci=[nloci],
                lociPos=locus_position,
                infoFields=['ind_id'])
        # tag the first generation so we can pass it to rc
        id_tagger.apply(pop)

        first_gen = pop.indInfo("ind_id")
        init_ts = msprime.simulate(2*len(first_gen), 
                                   length=max(locus_position))
        haploid_labels = [(k,p) for k in first_gen 
                                for p in (0,1)]
        node_ids = {x:j for x, j in zip(haploid_labels, init_ts.samples())}
        rc = ftprime.RecombCollector(ts=init_ts, node_ids=node_ids,
                                     locus_position=locus_position)
        recombinator = sim.Recombinator(intensity=recomb_rate,
                                        output=rc.collect_recombs,
                                        infoFields="ind_id")
        pop.evolve(
            initOps=[
                sim.InitSex()
            ]+init_geno,
            preOps=[
                sim.PyOperator(lambda pop: rc.increment_time() or True),
                # Must return true or false. True keeps whole population (?)
            ],
            matingScheme=mating_scheme_factory(recombinator, popsize, id_tagger),
            postOps=[
                sim.PyEval(r"'Gen: %2d\n' % (gen, )", step=1)
            ],
            gen=generations
        )
        return pop, rc
    return _make_pop


@pytest.mark.parametrize(('generations', 'popsize'), [
    (3, 5),
    (3, 10),
    (3, 20),
    (5, 5),
    (5, 10),
    (5, 20),
    (10, 5),
    (10, 10),
    (10, 20),
])
def test_simupop(make_pop, generations, popsize):
    print("Popsize: ", popsize)
    # replications = 1
    nsamples = 2
    length = 10
    nloci = 5
    locus_position = list(set(list(range(0, length, int(length/nloci))) + [length]))
    locus_position.sort()
    nloci = len(locus_position)
    recomb_rate = 0.05

    init_geno = [sim.InitGenotype(freq=[0.9, 0.1], loci=sim.ALL_AVAIL)]

    id_tagger = sim.IdTagger(begin=0)
    id_tagger.reset(startID=1)  # must reset - creating a new one doesn't
    pop,rc = make_pop(popsize, nloci, locus_position, id_tagger, init_geno,
                   recomb_rate, generations, length)
    locations = [pop.subPopIndPair(x)[0] for x in range(pop.popSize())]
    rc.add_diploid_samples(nsamples, pop.indInfo("ind_id"))

    check_tables(rc.args)

    print(rc.args)

    ts = rc.args.tree_sequence()

    print("coalescence records:")
    for x in ts.records():
        print(x)

    ts.simplify()

    print("trees:")
    for x in ts.trees():
        print(x)

@pytest.mark.parametrize(('generations', 'popsize', 'locus_position'), [
    (3, 10, [0.0, 1.0]),
    (3, 10, [0.0, 0.5, 1.0]),
    (3, 10, [0.0, 0.5, 0.6, 1.0]),
    (3, 10, [0.0, 0.1, 1.0]),
    (3, 10, [0.0, 0.001, 0.01, 0.1, 0.5, 1.0]),
    (10, 20, [0.0, 0.5, 1.0]),
    (10, 20, [k/1000 for k in range(1001)])
])
def test_recombination(make_pop, generations, popsize, locus_position):
    print("Popsize: ", popsize)
    print("locus positions: ", locus_position)
    nsamples = 2
    length = 1
    nloci = len(locus_position)
    recomb_rate = 1.0

    init_geno = [sim.InitGenotype(freq=[0.9, 0.1], loci=sim.ALL_AVAIL)]

    id_tagger = sim.IdTagger(begin=0)
    id_tagger.reset(startID=1)  # must reset - creating a new one doesn't
    pop,rc = make_pop(popsize, nloci, locus_position, id_tagger, init_geno,
                   recomb_rate, generations, length)
    locations = [pop.subPopIndPair(x)[0] for x in range(pop.popSize())]
    rc.add_diploid_samples(nsamples, pop.indInfo("ind_id"))

    print(rc.args)

    # check for uniformity of recombination breakpoints
    edges = rc.args.edgesets

    breakpoints = [x for x in edges.left if (x > 0.0) and (x < 1.0)]
    breakpoints += [x for x in edges.right if (x > 0.0) and (x < 1.0)]
    breaks = set(breakpoints)

    def check_partition(x, p):
        # compare to binomial approximation (by normal approx)
        n = len(x)
        k = len([y for y in x if y < p])
        z = (k - n * p) / math.sqrt(n * p * (1-p))
        return z

    z = [check_partition(breaks, k/10) for k in range(1,10)]
    print(z)
    assert(all([abs(zz) < 3.5 for zz in z]))

