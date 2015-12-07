#!/usr/bin/env python2.7
# File Location: /soe/thjmatthsoe/BME205_15/hw4
# Local py: #!/usr/bin/env python2.7
# Last Refactored: 11/12/15 (count_kmers modified for hw7)
# Contact: thjmatth@ucsc.edu
# Author: Thomas J Matthew
"""
markov.py makes kmer:count dicts from sequence data for downstream analysis

Functions:
    count_kmer (arguments in paren):
    Uses read_fasta function imported from fastIO to get sequence data
    from file (input) that is filtered on alphabet to ignore all other
    characters but the specified set (alphabet), breaking seqs into kmers of
    length order + 1 to capture the kmer context (order), while specifying
    start characters (start_char) and stop character (stop_char)
"""
import collections
from fastIO import read_fasta
from copy import deepcopy as dc
from itertools import product as prod


def count_kmer(fafile, alphabet='ABCDEFGHIJKLMNOPQRSTUVWXYZ', order=0,
               start_char=None, stop_char=None, kmin=None, kmax=None,
               pseudo=0):
    """
    counts kmers from a parsed sequence file

    :param inputlist: file, fasta sequences
    :param alphabet: str, all allowed bases
    :param order: int, markov order defaulting to 0 (iid model)
    :param start_char: str, single start char defaulting to '^'
    :param stop_char: str, single stop char defaulting to '$'
    :param kmin: int, minimum kmer length
    :param kmax: int, maximum kmer length
    :return: collections.counter object, kmer counts
    """
    seq = None
    seq_len = 0

    kmer_counts = collections.Counter(seq)
    for seqID, com, seq in read_fasta(fafile, alphabet, case='keep'):
        if (start_char and stop_char) and not (kmin and kmax):
            if order == 0:
                seq = seq + stop_char
            if order >= 1:
                seq = order * start_char + seq + stop_char * order
            for i in range(len(seq) - order):
                kmer_counts[seq[i:i + order + 1]] += 1

        # Below used with script score_palindromes
        else:
            seq_len += len(seq)
            for kmer_len in range(kmin, kmax + 1):
                for start in range(len(seq) - kmer_len + 1):
                    kmer_counts[seq[start:start + kmer_len]] += (1 + pseudo)

    return kmer_counts, seq_len

def count_kmer_table(table):
    """
    counts kmers from a table of kmer'\t'count, typically from stdin

    :param table: str, table of kmer counts as kmer'\t'count
    :return: dict, kmer counts as {'kmer':int}
    """
    # dictionary of {kmer:value}
    kmer_dict = {}
    for line in table:
        kmer_count = line.split()
        if len(kmer_count) == 2:
            kmer_count[1] = int(kmer_count[1])
            kmer_dict[kmer_count[0]] = kmer_count[1]
        else:
            raise Exception(
                "kmer count table from stdin must be kmer-tab-count")

    return kmer_dict


def gen_kmers(kmin, kmax, alphabet="ACGT"):
    """
    generates possible k-mers of range(kmin, kmax + 1)

    :param kmin: int, minimum kmer length
    :param kmax: int, maximum kmer length
    :param alphabet: str, accepted sequence alphabet for DNA, RNA, Amino Acids

    :return list of str, possible kmers
    """

    for n in xrange(kmin, kmax + 1):
        kmer_lis = [''.join(mer) for mer in prod(alphabet, repeat=n)]

    return kmer_lis

def counts_pseudo(kmer_count, order, alphabet, start_char, stop_char,
    pseudo):
    """
    counts all observed kmers (counted in training data) and possible kmers
    (incremented by psedocount), adds intersection as in a dict of dicts

    :param kmer_count: dict, kmer counts as {'kmer':int}
    :param order: int, context length for Markov chain
    :param alphabet: str, accepted sequence alphabet for DNA, RNA, Amino Acids
    :param start_char: str, single start char defaulting to '^'
    :param stop_char: str, single stop char defaulting to '$'
    :param pseudo: float or int, non-negative pseudocount value
    :return: dict of dict, context counts as{'context':{'finalchar':int}}
    """

    # tracks every final letter of each kmer
    char_count = {}
    for letter in alphabet + stop_char:  # creates pseudocounts from alphabet
        char_count[letter] = pseudo

    # tracks contexts of each final character
    contexts = {}

    # cartesian product of all possible kmers
    if order >= 0:
        for start_pos in range(order + 1):
            for mer in prod(alphabet, repeat=order - start_pos):
                contexts[start_char * start_pos + ''.join(mer)] = dc(
                    char_count)

    for kmer, freq in kmer_count.items():
        # prevents stop_char in context or start_char last (invalid seqs)
        if((order and kmer[0:-1][-1] != stop_char and kmer[-1] !=
            start_char) or order == 0):
            if kmer[0:-1] not in contexts.keys():
                raise Exception(
                    "K-mer context {} contains letters not in alphabet!"
                    .format(kmer[0:-1]))
            if kmer[-1] not in contexts[kmer[0:-1]].keys():
                raise Exception("K-mer contains {} but is not in alphabet!"
                    .format(kmer[-1]))
            contexts[kmer[0:-1]][kmer[-1]] += freq

    return contexts