#!/usr/bin/env python2.7
# File Location: /soe/thjmatthsoe/BME205_15/hw2
# Local py: #!/usr/bin/env python2.7
# Refactored: 11/21/15
# Contact: thjmatth@ucsc.edu
# Author: Thomas J Matthew
"""
fastIO helps with parsing, trimming, & writing common sequence files

Functions:
read_fasta:
    generates (seqID, com, seq) tuple from fasta or quality files
read_fastq:
    generates (seqID, com, seq) tuple from fastq files
read_fasta_with_quality:
    yields 'zipped' together parsed fasta and quality information
write_fastfile:
    writes a fasta/Q/qual file according to mostly widely accepted form
read_codon:
    parses codon preference table
trim_gen:
    trims sequence and quality data on first occurrence of low quality base
read_subst:
    parses base substitution matrix into dictionary for score lookups
read_codon_url:
    parses codon preference table in GCG-format from url
ncbi_genetic_code_table:
    parses genetic code versions from the NCBI
trim_gen:
    trims sequence and quality data on first occurrence of low quality base
"""

import sys
from itertools import izip
import string
import textwrap
import collections


def read_fasta(fainput, alphabet='ABCDEFGHIJKLMNOPQRSTUVWXYZ',
               case='keep', quality=False):
    """
    parses fasta / quality files and yields sequence information by generator

    :param fainput: file, opened fasta file object
    :param alphabet: str, accepted sequence alphabet for DNA, RNA, Amino Acids
    :param quality: bool, False is fasta and True is quality file
    :return: generator, tuple of str(seqID), str(com), str(seq), int(qual)
    """
    seqID, com = None,''
    if not quality: seq = '' 
    if quality: seq = []
    for line in fainput:
        if len(line) > 0:
            line = line.rstrip()
            if line.startswith('>'):
                if seqID is not None:
                    if not quality: seq = filter(lambda x: x.upper() in
                                                                alphabet,seq)
                    yield (seqID, com, seq)
                    if not quality: seq = ''
                    if quality: seq = [] 
                seqID = line[1:]
                sep = None  # sep is separator value
                if ' ' in seqID: sep = ' '
                elif ',' in seqID: sep = ','
                if sep is not None:
                    seqID, com =  seqID.split(sep,1)
            else:
                if quality == False:
                    seq += line
                    if case == 'upper': seq = seq.upper()
                    if case == 'lower': seq = seq.lower()
                if quality == True:
                    seq = seq + [int(value) for value in line.split()]

    if seqID is not None:
        if not quality: seq= filter(lambda x: x.upper() in alphabet,seq)
        yield (seqID, com, seq)

def read_fastq(fqinput, offset):
    """
    parses fastq files and yields sequence information by generator

    :param fqinput: file, opened fastq file object
    :param offset: int, quality encoding offset (usually 33 or 64)
    :return: generator, tuple(str(seqID)), str(com), str(seq), list(int(qual))
    """
    seqID, com, seq, qual = None,'','',''

    step = 1 # step is a variable that organizes the order fastq parsing
    # step = 1 parser at header line to extract seqID and/or comm
    # step = 2 parser at sequence line
    # step = 3 parser at quality line
    for line in fqinput:
        line = line.strip()
        if step == 1 and line.startswith('@'): # Step system from Nedda Saremi
            if seqID is not None:
                # Converts phred encoding from str to list of offset integer
                qual = [ord(char)-offset for char in qual]
                sep = None
                if string.whitespace in seqID: sep = string.whitespace
                if sep is not None:
                   seqID, com = seqID.split(sep,1)
                yield seqID, com, seq, qual
                seqID, com, seq, qual = None, '', '', ''  # Resets sequence
            seqID = line[1:]
            step = 2
            continue
        if step==2 and not line.startswith('@') and not line.startswith('+'):
            seq = seq + line
            continue
        if step == 2 and line.startswith('+'):
            step = 3
            continue
        if step == 3:
            qual = qual + line
            if len(seq) == len(qual):
                step = 1
                continue
            else:
                continue
    if seqID is not None:
        # Parses final sequence of file
        if len(qual) > 0:
            qual = [ord(char) - offset for char in qual]
        sep = None
        if string.whitespace in seqID:
            sep = string.whitespace
        if sep is not None:
            seqID, com = seqID.split(sep, 1)
        yield seqID, com, seq, qual


def read_fasta_with_quality(in_fasta, in_qual):
    """
    yields 'zipped' together parsed fasta and quality information

    :param in_fasta: generator of seqID, com, seq,
    :param in_qual: generator of seqID, com, qual
    :return: generator, tuple, of seqID, comment, seq, qual
    """
    for fa,fq in izip(read_fasta(in_fasta),read_fasta(in_qual, quality=True)):
        seqID, comment, seq, qual = fa[0],fa[1],fa[2],fq[2]
        yield(seqID, comment, seq, qual)


def write_fastfile(seqID, com, seq, qual, start_char, outlis):
    """
    writes a fasta/Q/qual file according to mostly widely accepted form

    :param seqID: str, identification of sequence
    :param com: str, comment on sequence
    :param seq: list of chars, sequence
    :param qual: str, quality of sequence
    :param start_char: str, '@' for fastQ ID, '>' for fasta ID
    :param outlis: list of strings, name of write file(s)
    :return: file object, used with sys.stdout or sys.stderr
    """
    if start_char == '@':
        fQname = outlis[0]
        fQname.write("{a}{b} {c}\n{d}\n+\n{e}\n".format(a=start_char,
                                                         b=seqID, c=com,
                                                         d=seq, e=qual))
    elif start_char == '>' and len(outlis) == 2:
        fAname, qname = outlis[0], outlis[1]
        seqwrap = textwrap.fill(seq, width=80)
        qualwrap = textwrap.fill(str(qual), width=80)
        fAname.write("{a}{b} {c}\n{d}\n".format(a=start_char, b=seqID,
                                                  c=com, d=seqwrap))
        qname.write("{a}{b} {c}\n{d}\n".format(a=start_char, b=seqID,
                                                  c=com, d=qualwrap))
    else:
        sys.stderr.write('Specify --out_fasta with --out_qual argument')


def read_codon(table):
    """
    parses codon preference table in GCG-format from file

    :param table: file, tab-delim with header, codons(col2) counts(col3)
    :return: dict, {str(codon) : float(count)}
    """
    codon_dict = {}
    for line in table:
        codon_count = line.split()
        if len(codon_count) == 5:
            codon_count[2] = float(codon_count[2])
            codon_dict[codon_count[1]] = codon_count[2]
        # else:
        #     raise Exception(
        #         "kmer count table from stdin must be kmer-tab-count")

    return codon_dict


def trim_seq(seqID, com, seq, qual, threshold):
    """
    trims sequence & quality data at first quality score

    :param gen: generator, of str(seqID), str(com), str(seq), int(qual)
    :param threshold: int, base quality threshold
    :return: generator, of seqID, comment, *trimmed* seq, *trimmed* qual
    """
    i = 0
    for q in qual:
        #i indexes satisfactory bases up to, or equaling, score
        if q < threshold:
            break
        else:
            i += 1
    yield seqID, com, seq[:i], qual[:i]


def read_subst(substfile):
    """
    parses base substitution matrix into dictionary for score lookups

    :param substfile: file, matrix with row & col labels and whitespace delim
    :return: dict, substitutions as tuple with score {('A', 'A'): 4...}
    :return: str, alphabet of joined col headers
    """

    sub_matrix = []
    for line in substfile:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        sub_matrix.append(line.split())

    col_label = sub_matrix.pop(0)
    alphabet = ''.join(col_label)

    row_label = [x[0] for x in sub_matrix]
    sub_matrix = [x[1:] for x in sub_matrix]

    subst_dict = {}
    for i in range(len(col_label)):
        for j in range(len(sub_matrix)):
            subst_dict[(col_label[i], row_label[j])] = int(sub_matrix[i][j])

    return subst_dict, alphabet


def read_codon_url(url):
    """
    parses codon preference table in GCG-format from url

    :param url: file, tab-delim with header, codons(col2) counts(col3)
    :return: counter, {str(codon) : float(count)}
    """

    # codon table bounded with <PRE> and </PRE> in HTML
    read_duration = None

    counts = collections.Counter()
    for line in url:
        if line.startswith('<PRE>'):
            read_duration = True
        elif line.startswith('</PRE>'):
            read_duration = False
        elif read_duration:
            split_line = line.split()

            if len(split_line) != 0 and split_line[0] != 'AmAcid':
                freq = float(split_line[2])
                codon = split_line[1]
                counts[codon] = freq

    return counts


def ncbi_genetic_code_table(prt, gene_code='Standard'):
    """
    parses genetic code versions from the NCBI

    :param prt: file, genetic codes from ftp://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt.
    :param gene_code: str, which genetic code to use, matching 1st entry name
    :return: dict, {str(codon) : float(count)}
    """

    # determines bounds of genetic code
    read_duration = False
    genetic_code_choices = set()
    codon_table = {}
    for line in prt:
        if not line.startswith('--'):
            line = line.rstrip()
            if line.startswith('  name "'):
                line = line.lstrip('  name "').rstrip('" ,')
                genetic_code_choices.add(line)

                if gene_code not in genetic_code_choices:
                    gene_code = 'Standard'
                    sys.stderr.write('WARNING: genetic_code not valid '
                                      'continuing with default Standard code')

                if line == gene_code:
                    read_duration = True
                    line = prt.next()
                else:
                    read_duration = False

            if read_duration:
                if line.startswith('  ncbieaa'):
                    line = line.lstrip('  ncbieaa  "').rstrip('",')
                    ncbieaa = line
                elif line.startswith('  -- Base1'):
                    line = line.lstrip('  -- Base1')
                    Base1 = line
                elif line.startswith('  -- Base2'):
                    line = line.lstrip('  -- Base2')
                    Base2 = line
                elif line.startswith('  -- Base3'):
                    line = line.lstrip('  -- Base3')
                    Base3 = line

    for index, aminobase in enumerate(ncbieaa):
        codon_table[Base1[index] + Base2[index] + Base3[index]] = aminobase

    return codon_table


###
# SURPLUS FUNCTIONS

def trim_gen(gen, threshold):
    """
    trims sequence & quality data at first quality score from generator

    call it like this from main:
    trimmed = trim_gen(needtrim, options.min_qual)

    :param gen: generator, of str(seqID), str(com), str(seq), int(qual)
    :param threshold: int, base quality threshold
    :return: generator, of seqID, comment, *trimmed* seq, *trimmed* qual
    """
    for seqID, com, seq, qual in gen:
        i = 0
        for q in qual:
            #i indexes satisfactory bases up to, or equaling, score
            if q < threshold:
                break
            else:
                i += 1
        yield seqID, com, seq[:i], qual[:i]