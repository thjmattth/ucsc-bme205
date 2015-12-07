#!/usr/bin/env python2.7
# File Location: /soe/thjmatthsoe/BME205_15/hw8
# Local py: #!/usr/bin/env python2.7
# Contact: thjmatth@ucsc.edu
# Author: Thomas J Matthew
"""

align.py defines global & local sequence alignment classes with affine gap cost

Classes:
    Alignment:      initializes affine gap costs
    GlobalAlign:    Needleman-Wunsch Algorithm
                    global alignment of two sequences with affine gap cost
    LocalAlign:     Smith-Waterman Algorithm
                    local alignment of two sequences with affine gap cost

Recurrance Relation:
    M_ij = best alignment ending with Rj alined to Cj
         = S_R[i],C[j]
           + max (  0
                    M_i-1,j-1
                    Ir_i-1,j-1
                    Ic_i-1,j-1 )
    Ir_ij = best alignment of R_0....R_i
                        to    C_0....C_j with R_i aligned to gap (no match C)
          = max  (  M_i-1,j - g
                    Ir_i-1,j - e
                    Ic_i-1,j - d
    Ic_ij = best alignment of R_0....R_i
                        to    C_0....C_j with C_i aligned to gap (no match R)
          = max  (  M_i,j-1 - g
                    Ir_i,j-1 - d
                    Ic_i,j-1 - e

Boundary Conditions:
        M_-1,j  =   - infinity
        M_i,-1  =   - infinity
        Ir_-1,j =   - infinity
        Ir_i,-1 =   - ( g * e(i) )
                    R[0]--------R[i]
                    ---------------
                    alignment is there R-0 through R_i with nothing
                    starting a gap at the beginning and extending gap
        Ic_i,-1 =   - infinity
        Ic_-1,j =   - ( g * e(i) )
                    ---------------
                    C[0]--------C[i]
Notation:
    M  :  Match
    Ir :  Insert in Row
    Ic :  Insert in Column
    S  :  Substitution score
    g  :  Start Gap
    e  :  Extend Gap
    d  :  Double Gap (gap continues on other strand)

Resources:
    Biological Sequence Analysis:
    Probabilistic Models of Proteins and Nucleic Acids
    from Cambridge University Press
    by R. Durbin, S. Eddy, A. Krogh, and G. Mitchison
    Page 17-28

"""

from __future__ import division, print_function
import numpy as np
from operator import itemgetter


class Alignment(object):

    """
    initializes alignment attributes and methods


    """
    def __init__(self, subst, start=12, extend=1, double=3):
        """
        initializes affine gap costs

        :param subst: dict, substitutions as tuple with score {('A', 'A'): 4...}
        :param start: gap start penalty, defaults to 12
        :param extend: gap extend penalty, defaults to 1
        :param double: int, gap extension penalty across strands defaults to 3
        """
        self.subst = subst
        self.start = start
        self.extend = extend
        self.double = double

        self.neginf = -10 ** 10


class GlobalAlign(Alignment):
    """ 
    global alignment of two sequences with affine gap cost (Needleman-Wunsch)

    Functions:
        align:
                fills M, IR, and IC alignment matrices as np arrays
                returns highest-scoring alignment

        traceback_col_seq:
                performs traceback on global alignment
                master-slave alignment sequence in A2M
    """

    def align(self, row_seq, col_seq):
        """
        align as global

        :param row_seq: str, sequence 1 (usually the master)
        :param col_seq: str, sequence 2 (usually the slave)
        :return: np.array, int scores in matrix 0 + row_seq by 0 + col_seq
        """

        # set dimensions of zeroed numpy array for base case and filling
        matrix_dim = [len(row_seq) + 1, len(col_seq) + 1]

        # http://www.biorecipes.com/DynProgBasic/code.html
        # http://stackoverflow.com/q/568962/4154548
        self.M = np.zeros(shape=matrix_dim, dtype=int)
        self.IR = np.zeros(shape=matrix_dim, dtype=int)
        self.IC = np.zeros(shape=matrix_dim, dtype=int)

        for i, char in enumerate(row_seq):
            self.M[i + 1][0] = self.neginf
            self.IR[i + 1][0] = -int(self.start + (i * self.extend))
            self.IC[i + 1][0] = self.neginf
            
        for j, char in enumerate(col_seq):
            self.M[0][j + 1] = self.neginf
            self.IR[0][j + 1] = self.neginf
            self.IC[0][j + 1] = -int(self.start + (j * self.extend))

        for i, row_char in enumerate(row_seq):
            for j, col_char in enumerate(col_seq):
                self.M[i + 1][j + 1] = self.subst[(row_char, col_char)] \
                                   + max(
                                        self.M[i][j],
                                        self.IR[i][j],
                                        self.IC[i][j])
                self.IR[i + 1][j + 1] = max(
                                        self.M[i][j + 1] - self.start,
                                        self.IR[i][j + 1] - self.extend,
                                        self.IC[i][j + 1] - self.double)
                self.IC[i + 1][j + 1] = max(
                                        self.M[i + 1][j] - self.start,
                                        self.IR[i + 1][j] - self.double,
                                        self.IC[i + 1][j] - self.extend)

        self.row_seq, self.col_seq = row_seq, col_seq

        global_max = max(self.M[i + 1][j + 1],
                        self.IR[i + 1][j + 1],
                        self.IC[i + 1][j + 1])

        return global_max

    def traceback_col_seq(self):
        """
        performs traceback on global alignment

        :return str, master-slave alignment sequence in A2M
        """
        # records traceback as a list aligned bases in reverse
        trace_list = []

        # traceback position
        row = len(self.row_seq)
        col = len(self.col_seq)

        #  matrix of trace with highest score (M:0, IR:1, IC:1)
        matricies = [self.M[row][col], self.IR[row][col], self.IC[row][col]]
        enumat = [pair for pair in enumerate(matricies)]
        hi_score_matrix, hi_score = max(enumat, key=itemgetter(1))

        # http://www.biorecipes.com/DynProgBasic/code.html
        while row != 0 or col != 0:

            # hi_score_matrix is M
            if hi_score_matrix == 0:
                trace_list += self.col_seq[col - 1].upper()
                subst = self.subst[(self.row_seq[row - 1],
                                    self.col_seq[col - 1])]
                row -= 1
                col -= 1

                if hi_score == subst + self.M[row][col]:
                    hi_score_matrix = 0
                    hi_score = self.M[row][col]
                elif hi_score == subst + self.IR[row][col]:
                    hi_score_matrix = 1
                    hi_score = self.IR[row][col]
                elif hi_score == subst + self.IC[row][col]:
                    hi_score_matrix = 2
                    hi_score = self.IC[row][col]

            # hi_score_matrix is IR
            elif hi_score_matrix == 1:
                trace_list += '-'
                row -= 1
                if hi_score + self.start == self.M[row][col]:
                    hi_score_matrix = 0
                    hi_score = self.M[row][col]
                elif hi_score + self.extend == self.IR[row][col]:
                    hi_score_matrix = 1
                    hi_score = self.IR[row][col]
                elif hi_score + self.double == self.IC[row][col]:
                    hi_score_matrix = 2
                    hi_score = self.IC[row][col]

            # hi_score_matrix is IC
            elif hi_score_matrix == 2:
                trace_list += self.col_seq[col - 1].lower()
                col -= 1
                if hi_score + self.start == self.M[row][col]:
                    hi_score_matrix = 0
                    hi_score = self.M[row][col]
                elif hi_score + self.double == self.IR[row][col]:
                    hi_score_matrix = 1
                    hi_score = self.IR[row][col]
                elif hi_score + self.extend == self.IC[row][col]:
                    hi_score_matrix = 2
                    hi_score = self.IC[row][col]

        forward_trace = ''.join(trace_list[::-1])

        return forward_trace


class LocalAlign(Alignment):

    """
    local alignment of two sequences with affine gap cost (Smith-Waterman)

    Functions:
        align:
                fills M, IR, and IC alignment matrices as np arrays
                returns highest-scoring alignment

        traceback_col_seq:
                performs traceback on global alignment
                master-slave alignment sequence in A2M

    """

    def align(self, row_seq, col_seq):
        """
        align as local

        :param row_seq: str, sequence 1 (usually the master)
        :param col_seq: str, sequence 2 (usually the slave)
        :return: np.array, int scores in matrix 0 + row_seq by 0 + col_seq
        """

        # set dimensions of zeroed numpy array for base case and filling
        matrix_dim = [len(row_seq) + 1, len(col_seq) + 1]

        # http://www.biorecipes.com/DynProgBasic/code.html
        # http://stackoverflow.com/q/568962/4154548
        self.M = np.zeros(shape=matrix_dim, dtype=int)
        self.IR = np.zeros(shape=matrix_dim, dtype=int)
        self.IC = np.zeros(shape=matrix_dim, dtype=int)

        # stores first maximum alignment pair
        self.best_row, self.best_col = 0, 0

        for i, char in enumerate(row_seq):
            self.M[i + 1][0] = self.neginf
            self.IR[i + 1][0] = self.neginf
            self.IC[i + 1][0] = self.neginf
            
        for j, char in enumerate(col_seq):
            self.M[0][j + 1] = self.neginf
            self.IR[0][j + 1] = self.neginf
            self.IC[0][j + 1] = self.neginf

        for i, row_char in enumerate(row_seq):
            for j, col_char in enumerate(col_seq):
                self.M[i + 1][j + 1] = self.subst[(row_char,col_char)] \
                                            + max(0,
                                            self.M[i][j],
                                            self.IR[i][j],
                                            self.IC[i][j])
                self.IR[i + 1][j + 1] = max(
                                            self.M[i][j + 1] - self.start,
                                            self.IR[i][j + 1] - self.extend,
                                            self.IC[i][j + 1] - self.double)
                self.IC[i + 1][j + 1] = max(
                                            self.M[i + 1][j] - self.start,
                                            self.IR[i + 1][j] - self.double,
                                            self.IC[i + 1][j] - self.extend)

                if self.M[i + 1][j + 1] > self.M[self.best_row][self.best_col]:
                    self.best_row = i + 1
                    self.best_col = j + 1

        self.row_seq = row_seq
        self.col_seq = col_seq

        return self.M[self.best_row][self.best_col]

    def traceback_col_seq(self):
        """
        performs traceback on local alignment

        :return str, master-slave alignment sequence in A2M
        """
        # records traceback as a list aligned bases in reverse
        trace_list = []

        # traceback position
        row = len(self.row_seq)
        col = len(self.col_seq)

        # appends insertions/deletions to best alignment in the local case
        while row > self.best_row:
            trace_list += '-'
            row -= 1
        while col > self.best_col:
            trace_list += self.col_seq[col - 1].lower()
            col -= 1

        # matrix path of highest score (M:0, IR:1, IC:1)
        hi_score_matrix = 0
        hi_score = self.M[self.best_row][self.best_col]

        # trace indicator
        trace_done = False

        # http://www.biorecipes.com/DynProgBasic/code.html
        while (row != 0 or col != 0) and not trace_done:

            # hi_score_matrix is M
            if hi_score_matrix == 0:
                trace_list += self.col_seq[col - 1].upper()
                subst = self.subst[(self.row_seq[row - 1],
                                    self.col_seq[col - 1])]
                row -= 1
                col -= 1

                if hi_score == subst:
                    trace_done = True
                elif hi_score == subst + self.M[row][col]:
                    hi_score_matrix = 0
                    hi_score = self.M[row][col]
                elif hi_score == subst + self.IR[row][col]:
                    hi_score_matrix = 1
                    hi_score = self.IR[row][col]
                elif hi_score == subst + self.IC[row][col]:
                    hi_score_matrix = 2
                    hi_score = self.IC[row][col]

            # hi_score_matrix is IR
            elif hi_score_matrix == 1:
                trace_list += '-'
                row -= 1

                if hi_score + self.start == self.M[row][col]:
                    hi_score_matrix = 0
                    hi_score = self.M[row][col]
                elif hi_score + self.extend == self.IR[row][col]:
                    hi_score_matrix = 1
                    hi_score = self.IR[row][col]
                elif hi_score + self.double == self.IC[row][col]:
                    hi_score_matrix = 2
                    hi_score = self.IC[row][col]

            # hi_score_matrix is IC
            elif hi_score_matrix == 2:
                trace_list += self.col_seq[col-1].lower()
                col -= 1
                if hi_score + self.start == self.M[row][col]:
                    hi_score_matrix = 0
                    hi_score = self.M[row][col]
                elif hi_score + self.double == self.IR[row][col]:
                    hi_score_matrix = 1
                    hi_score = self.IR[row][col]
                elif hi_score + self.extend == self.IC[row][col]:
                    hi_score_matrix = 2
                    hi_score = self.IC[row][col]

        while col > 0:
            trace_list += self.col_seq[col - 1].lower()
            col -= 1
        while row > 0:
            trace_list += '-'
            row -= 1

        forward_trace = ''.join(trace_list[::-1])

        return forward_trace
