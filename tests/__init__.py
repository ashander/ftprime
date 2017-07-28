import unittest

class FtprimeTestCase(unittest.TestCase):
    """
    Superclass of test cases containing common utilities.
    """
    random_seed = 123456

    def check_trees(self, tsa, tsb, npos=20):
        # check trees at a bunch of positions agree
        self.assertEqual(tsa.sequence_length, tsb.sequence_length)
        for x in tsa.dump_tables():
            print(x)
        for x in tsb.dump_tables():
            print(x)
        self.assertEqual(tsa.num_nodes, tsb.num_nodes)
        self.assertEqual(len(tsa.samples()), len(tsb.samples()))
        for a, b in zip(tsa.samples(), tsb.samples()):
            self.assertEqual(a, b)
        check_positions = [tsa.sequence_length*k/npos for k in range(npos)]
        tsa_t = tsa.trees()
        tsb_t = tsb.trees()
        ta = next(tsa_t)
        tb = next(tsb_t)
        for pos in check_positions:
            while pos > ta.interval[1]:
                ta = next(tsa_t)
            while pos > tb.interval[1]:
                tb = next(tsb_t)
            self.check_tree(ta, tb, tsa.samples())
        return True

    def check_tree(self, ta, tb, samples):
        for u in samples:
            for v in samples:
                self.assertEqual(ta.mrca(u,v), tb.mrca(u,v))
        return True

    def assertArrayEqual(self, x, y):
        self.assertEqual(len(x), len(y))
        for k in range(len(x)):
            self.assertEqual(x[k], y[k])
