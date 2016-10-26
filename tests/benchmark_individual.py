from ftprime import Individual, Chromosome

def reference_gametes():
    p1_id = 1
    p2_id = 3
    c1 = Chromosome(p1_id, .5)
    c2 = Chromosome(p2_id, .5)
    i = Individual((c1, c2))
    i_gams = [c for c in i.gametes(5)]
    return i_gams


def test_benchmark_gametes(benchmark):
    benchmark(reference_gametes)
