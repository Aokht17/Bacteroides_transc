import os
import re
import pandas as pd
from Bio.Seq import Seq
from Bio import pairwise2


def find_coord(fasta, prefix, plus_motif=r"GGA.....TCC", minus_motif=r"CCT.....AGG"):
    """
    finds coordinates of given motifs in both strands
    writes output in 2 separate tsv files
    :param fasta: full path to genome fasta
    :param prefix: prefix for output files
    :param plus_motif: palindromic consensus motif for + strand
    :param minus_motif: palindromic consensus motif for - strand
    :return: None
    """
    with open(fasta, 'r', encoding='utf-8') as inp:
        data = str()
        for line in inp:
            if not line.startswith('>'):
                data += line.strip()

    length = len(data)
    plus_hits = re.findall(plus_motif, data)
    minus_hits = re.findall(minus_motif, data)
    plus_hits = list(set(plus_hits))
    minus_hits = list(set(minus_hits))
    step = len(plus_motif)-1
    for motif in plus_hits:
        if motif in data:
            n = data.find(motif)
            with open(f'{prefix}_motif_plus_coord.tsv', 'a') as out:
                out.write(f'{n}\t{n+step}\n')
            while n < length and n != -1:
                n = data.find(motif, n+1)
                if n != -1:
                    with open(f'{prefix}_motif_plus_coord.tsv', 'a') as out:
                        out.write(f'{n}\t{n+step}\n')

    for motif in minus_hits:
        if motif in data:
            n = data.find(motif)
            with open(f'{prefix}_motif_minus_coord.tsv', 'a') as out:
                out.write(f'{n}\t{n+step}\n')
            while n < length and n != -1:
                n = data.find(motif, n+1)
                if n != -1:
                    with open(f'{prefix}_motif_minus_coord.tsv', 'a') as out:
                        out.write(f'{n}\t{n+step}\n')


def kilobase_take(fasta, coord, orient, prefix):
    """
    extracts 1 kb regions downstream of each pair of coordinates given
    writes output in 2 separate fasta files
    :param fasta: full path to genome fasta
    :param coord: tsv file with coordinate pairs
    :param orient: + or - according to coord file
    :param prefix: prefix for output files
    :return: None
    """
    with open(fasta, 'r', encoding='utf-8') as inp:
        data = str()
        for line in inp:
            if not line.startswith('>'):
                data += line.strip()
    coordinate = pd.read_csv(coord, sep='\t', names=['start', 'end'])
    if orient == '+':
        for end in coordinate['end'].to_list():
            gen_slice = data[end+1: end+1001]
            with open(f'{prefix}_reg_net_plus_1k.fasta', 'a') as out:
                out.write(f'>{end}\n{gen_slice}\n')
    if orient == '-':
        for start in coordinate['start'].to_list():
            gen_slice = data[start-1001: start-1]
            with open(f'{prefix}_reg_net_minus_1k.fasta', 'a') as out:
                out.write(f'>{start}\n{gen_slice}\n')


def fragment_extract(fasta, coord, prefix):
    """
    extracts genome fragment for each pair of coordinates given
    writes output in fasta file
    :param fasta: full path to genome fasta
    :param coord: tsv file with coordinate pairs
    :param prefix: prefix for output files
    :return: None
    """
    with open(fasta, 'r', encoding='utf-8') as inp:
        data = str()
        for line in inp:
            if not line.startswith('>'):
                data += line.strip()
    coordinate = pd.read_csv(coord, sep='\t', names=['start', 'end'])
    for _, row in coordinate.iterrows():
        gen_slice = data[row.start-1:row.end]
        with open(f'{prefix}_fragments.fasta', 'a') as out:
            out.write(f'>{row.start}-{row.end}\n{gen_slice}\n')


def upstream_extractor(gene, cds_file, fasta, len_up):
    """
    extracts upstream region of given gene
    writes output to {gene}_seq.txt
    :param gene: gene of interest
    :param cds_file: path to cdf file
    :param fasta: path to genome fasta file
    :param len_up: the length of upstream region
    :return: None
    """
    with open(cds_file, 'r') as cds:
        for lines in cds:
            line = cds.readline()
            if line.startswith('>') and gene in line:
                strand = '-' if 'location=complement' in line else '+'
                line = cds.readline()
                while not line.startswith('>'):
                    with open(f'{gene}_seq.txt', 'a') as out:
                        out.write(f'{line}')
                    line = cds.readline()
                break
    try:
        with open(f'{gene}_seq.txt', 'r') as file:
            seq = str()
            for i in file.readlines():
                seq += i.strip()
        with open(fasta, 'r') as fasta:
            data = str()
            for line in fasta:
                if not line.startswith('>'):
                    data += line.strip()
        # take reverse complement if location=complement
        start = data.find(seq) if strand == '+' else data.find(str(Seq(seq).reverse_complement()))
        end = start + len(seq)
        os.remove(f'{gene}_seq.txt')
        if strand == '-':
            return data[end+1:end+len_up+1]
        else:
            return data[start-len_up-1:start-1]
    except:
        return None


def id_below_maxid_perc(el1, el2, max_percent_id):
    """ Aligns two sequence elements and determines whether they are
        more than %ID identical (false) or not (true)
        Scoring: Match:+2, Mismatch:-1, GapO: -2, GapE: -0.2
    """

    al = pairwise2.align.globalms(el1, el2, 2, 0, -2, -.5, one_alignment_only=True, penalize_end_gaps=False)

    matches = 0
    gapless = 0
    # for each position in the alignment
    for ch_pair in zip(al[0][0], al[0][1]):
        # if this is a non-gapped position
        if '-' not in ch_pair:
            # if it's a match, count it
            if ch_pair[0] == ch_pair[1]:
                matches = matches + 1
            gapless = gapless + 1

    perID = float(matches) / float(gapless)

    # return true or false depending on percent identity
    if perID * 100 <= float(max_percent_id):
        return True
    else:
        return False


def remove_redundants(promoters_list, max_percent_identity):
    """ Takes a list of promoter sequences and returns a list of promoters with
        <= max_percent_identity.
        The identitiy values are computed using id_below_maxid_perc function defined
        above
    """

    removed_ids = []
    nr_promoters = []

    for i in range(0, len(promoters_list) - 1):
        if promoters_list[i] not in removed_ids:
            nr_promoters.append(promoters_list[i])
            for j in range(i + 1, len(promoters_list)):
                if not id_below_maxid_perc(promoters_list[i], promoters_list[j], max_percent_identity):
                    removed_ids.append(promoters_list[j])

    return nr_promoters


def unique_filter(file, max_percent_identity=75):
    with open(file, 'r') as inp:
        seq_list = []
        for lines in inp:
            seq = ''
            line = inp.readline()
            if not line.startswith('>'):
                seq += line.strip()
            seq_list.append(seq)
        filt_list = remove_redundants(seq_list, max_percent_identity)
        return filt_list
