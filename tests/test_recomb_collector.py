import pytest
import simuPOP as sim

from ftprime import RecombCollector
# from http://simupop.sourceforge.net/manual_svn/build/userGuide_ch3_sec4.html

def check_record_order(args):
    for ind in args:
        if len(args[ind][1]) > 0:
            x = 0.0
            for cr in args[ind][1]:
                if x > cr.left or cr.left >= cr.right:
                    print("bad record order:", ind, args[ind])
                    raise ValueError
                x = cr.right


def check_tables(args):
    nodes = args.node_table()
    assert(nodes.num_rows == args.num_nodes)
    edgesets = args.edgeset_table()
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
    lambda recombinator, popsize: sim.RandomMating(
        ops=[
            sim.IdTagger(),
            recombinator
        ]),
    # Overlapping generations mating system
    lambda recombinator, popsize: sim.HeteroMating([sim.RandomMating(
        ops=[
            sim.IdTagger(),
            recombinator
        ]),
        sim.CloneMating()],
        subPopSize=popsize * 2)
    ])
def make_pop(request):
    # request.param stores a lambda function to make mating scheme
    # each test that uses this fixture will be run for both entries in 'params'
    mating_scheme_factory = request.param

    def _make_pop(popsize, nloci, locus_position, id_tagger, init_geno,
                  recomb_rate, rc, generations):
        sim.setOptions(seed=111)
        recombinator = sim.Recombinator(intensity=recomb_rate,
                                        output=rc.collect_recombs,
                                        infoFields="ind_id")
        pop = sim.Population(
                size=[popsize],
                loci=[nloci],
                lociPos=locus_position,
                infoFields=['ind_id'])
        pop.evolve(
            initOps=[
                sim.InitSex(),
                id_tagger
            ]+init_geno,
            preOps=[
                sim.PyOperator(lambda pop: rc.increment_time() or True),
                # Must return true or false. True keeps whole population (?)
            ],
            matingScheme=mating_scheme_factory(recombinator, popsize),
            postOps=[
                sim.PyEval(r"'Gen: %2d\n' % (gen, )", step=1)
            ],
            gen=generations
        )
        return pop
    return _make_pop

# # occasional failures marked below
# # ValueError: Parent 3's birth time has not been recorded with .add_individual()
# @pytest.mark.parametrize(('generations', 'popsize'), [
#     (3, 5),  # stochastic fail
#     (3, 10), # stochastic fail
#     (3, 20),
#     (5, 5),
#     (5, 10),
#     (5, 20),
#     (10, 5),
#     (10, 10),
#     (10, 20),
# ])

@profile
def test_simupop(make_pop, generations, popsize):
    print("Popsize: ", popsize)
    # replications = 1
    nsamples = 2
    length = 10
    nloci = 5
    locus_position = list(range(0, length, int(length/nloci)))
    recomb_rate = 0.05

    rc = RecombCollector(
            nsamples=nsamples, generations=generations, N=popsize,
            ancestor_age=10, length=length, locus_position=locus_position)

    init_geno = [sim.InitGenotype(freq=[0.9, 0.1], loci=sim.ALL_AVAIL)]

    id_tagger = sim.IdTagger(begin=0)
    id_tagger.reset(startID=1)  # must reset - creating a new one doesn't
    pop = make_pop(popsize, nloci, locus_position, id_tagger, init_geno,
                   recomb_rate, rc, generations)
    locations = [pop.subPopIndPair(x)[0] for x in range(pop.popSize())]
    rc.add_diploid_samples(pop.indInfo("ind_id"), locations)

    # check_record_order(rc.args)
    # check_tables(rc.args)

    # for x in rc.args:
    #     print(rc.args[x])
    # print(rc.args.node_table())
    # print(rc.args.edgeset_table())

    # ts = rc.args.tree_sequence()

    # print("coalescence records:")
    # for x in ts.records():
    #     print(x)

    # ts.simplify(samples=list(range(nsamples)))

    # print("trees:")
    # for x in ts.trees():
    #     print(x)


if __name__ == '__main__':
    class MockRequest(object):
        def __init__(self):
            pass
    mr = MockRequest()
    mr.param = lambda recombinator, popsize: sim.HeteroMating([sim.RandomMating(
        ops=[
            sim.IdTagger(),
            recombinator
        ]),
        sim.CloneMating()],
        subPopSize=popsize * 2)

    test_simupop(make_pop(mr), generations=200, popsize=1000)
