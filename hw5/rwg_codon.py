# File Location: /soe/thjmatthsoe/BME205_15/hw5
# Local py: #!/usr/bin/env python2.7
# Last Refactored: 10/29/15
# Contact: thjmatth@ucsc.edu
# Author: Thomas J Matthew
# Ref:
# eli.thegreenplace.net/2010/01/22/weighted-random-generation-in-python/
# http://stackoverflow.com/questions/3679694/a-weighted-version-of-random-choice

from random import randint as ri
from bisect import bisect_right as br


class RandomCodonGenerator(object):
    """
    RandomCodonGenerator is a weighted random generator.

    instantiates a dictionary object with codons as keys and counts as
    values. Every call to the class returns a random codon weighted by the
    counts at instantiation.
    """

    def __init__(self, counts):

        self.totals = []  # sum of the first k+1 values in 'counts'
        self.counts_list = counts.items()  # lists dictionary for ordering
        running_total = 0

        for item in self.counts_list:
            running_total += item[1]
            self.totals.append(running_total)

        self.total = self.totals[-1]

    def gen_rand_codon(self):

        rnd = ri(0, self.total - 1)
        index = br(self.totals, rnd)
        return self.counts_list[index][0]

    def __call__(self):
        return self.gen_rand_codon()