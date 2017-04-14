import msprime


def test_time_attribute():
    node = msprime.NodeTable()
    node.time
    print(msprime.__version__)


if __name__ == '__main__':
    test_time_attribute()
