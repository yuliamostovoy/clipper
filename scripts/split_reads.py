import sys
import re
import pysam


def main():
    if len(sys.argv) != 2:
        sys.stderr.write("Usage: program.py bamfile\n")
        sys.exit(1)

    bamfile = pysam.AlignmentFile(sys.argv[1], 'rb')
    # min_clipping = 500

    # print("ref\tpos\tpos+1\tmapq\talignment-length\tsoft-clipping-length\treadname")
    for read in bamfile:
        if read.is_secondary:  # should leave only primary and supplementary alignments
            continue
        if not read.has_tag('SA'):
            continue
        alignments = parse_bam_record(read)  # list of dicts, entries: r, rs, re, qs, qe, mq, read_length, aligned_length, strand

        if read.is_reverse:
            alignments = aln_reverse(alignments)
        index_qs = alignments[-1]['qs']
        index_qe = alignments[-1]['qe']
        alns_sorted = sorted(alignments, key=lambda k: min(k['qs'], k['qe']))

        for i in range(0,len(alns_sorted)):
            aln = alns_sorted[i]
            if aln['qs'] == index_qs and aln['qe'] == index_qe and aln['r'] == read.reference_name:  # check if this is the index alignment
                if i > 0:  # if a previous alignment exists
                    aln_SA = alns_sorted[i - 1]
                    if aln_SA['strand'] == '+':
                        SA = "%s:%d%s" % (aln_SA['r'], aln_SA['re'], '-')
                    else:
                        SA = "%s:%d%s" % (aln_SA['r'], aln_SA['rs'], '+')
                    print("%s\t%d\t%d\t.\t.\t%s\t%d\t%d\t%s\t%s" % (aln['r'], aln['rs'], aln['rs'] + 1, '-', aln['mq'], aln['aligned_length'], read.query_name, SA))
                if i < len(alns_sorted)-1:  # if a subsequent alignment exists
                    aln_SA = alns_sorted[i + 1]
                    if aln_SA['strand'] == '+':
                        SA = "%s:%d%s" % (aln_SA['r'], aln_SA['rs'], aln_SA['strand'])
                    else:
                        SA = "%s:%d%s" % (aln_SA['r'], aln_SA['re'], aln_SA['strand'])
                    print("%s\t%d\t%d\t.\t.\t%s\t%d\t%d\t%s\t%s" % (aln['r'], aln['re'], aln['re'] + 1, '+', aln['mq'], aln['aligned_length'], read.query_name, SA))


def aln_reverse(alignments):
    # flip the orientation of all alignments for a read
    ori_flip = {'+': '-', '-': '+'}
    for aln in alignments:
        aln['qs'] = aln['read_length'] - aln['qs']
        aln['qe'] = aln['read_length'] - aln['qe']
        aln['strand'] = ori_flip[aln['strand']]
    return alignments


def parse_bam_record(record):
    chrom = record.reference_name
    rstart = record.reference_start
    mq = record.mapping_quality
    cigar = record.cigartuples

    if record.is_reverse:
        strand = "-"
    else:
        strand = "+"

    alignments = []  # list of dicts

    alignments += parse_SA_field(record.get_tag("SA"))
    alignments.append(read_cigar(cigar, chrom, rstart, strand, mq))

    read_length = alignments[-1]['read_length']

    for aln in alignments:
        if aln['read_length'] != read_length:
            sys.stderr.write("Warning: read length of primary and supplementary alignments do not match for this read (calculated using cigar strings)\n" % record.query_name)
            sys.exit(1)

    return alignments


def parse_SA_field(sa):
    alignments = []
    aligns = sa.rstrip(';').split(";")
    for align in aligns:
        fields = align.split(",")
        if len(fields) >= 6:
            chrom = fields[0]
            rstart = int(fields[1])
            raw_cigar = fields[3]
            strand = fields[2]
            mq = int(fields[4])
            alignments.append(read_cigar(parse_cigar(raw_cigar), chrom, rstart, strand, mq))
        else:
            sys.stderr.write("ignoring alternate alignment because it doesn't have all 6 columns: %s\n" % align)

    return alignments  # returns an alignment list of dicts


def parse_cigar(raw_cigar):
    # returns cigar as list of tuples just like in pysam
    translate = {'M': 0, 'I': 1, 'D': 2, 'N': 3, 'S': 4, 'H': 5, 'P': 6, '=': 7, 'X': 8, 'B': 9}
    cigar = []
    for num, element in re.findall('(\d+)(\D+)', raw_cigar):
        if element not in translate:
            sys.stderr.write("CIGAR element not found: %s\n" % element)
            sys.exit(1)
        cigar.append((translate[element], int(num)))
    return cigar


def read_cigar(cigar, chrom, rstart, strand, mq):
    # returns an alignment dict
    coordinates = cigar_coords(cigar)

    alignment = {'r': chrom, 'rs': rstart, 're': rstart + coordinates['ref_alignment_length']}

    if strand == "+":
        alignment['qs'] = coordinates['front_padding_length']
        alignment['qe'] = coordinates['front_padding_length'] + coordinates['read_alignment_length']
    else:
        alignment['qe'] = coordinates['end_padding_length']
        alignment['qs'] = coordinates['end_padding_length'] + coordinates['read_alignment_length']

    alignment['read_length'] = coordinates['front_padding_length'] + coordinates['read_alignment_length'] + coordinates['end_padding_length']
    alignment['mq'] = mq
    alignment['aligned_length'] = coordinates['read_alignment_length']
    alignment['strand'] = strand

    # skipping the part where we add a 'path' field to alignment because we don't need it here
    return alignment


def cigar_coords(cigar):
    coords = {'read_alignment_length': 0, 'ref_alignment_length': 0, 'front_padding_length': 0, 'end_padding_length': 0}

    no_matches_yet = True
    for element in cigar:
        num = element[1]
        if element[0] == 5:  # H
            if no_matches_yet:
                coords['front_padding_length'] += num
            else:
                coords['end_padding_length'] += num
        elif element[0] == 4:  # S
            if no_matches_yet:
                coords['front_padding_length'] += num
            else:
                coords['end_padding_length'] += num
        elif element[0] in [0, 7, 8]:  # M, =, X
            no_matches_yet = False
            coords['read_alignment_length'] += num
            coords['ref_alignment_length'] += num
        elif element[0] == 1:  # I
            no_matches_yet = False
            coords['read_alignment_length'] += num
        elif element[0] in [2, 3, 6]:  # D, N, P
            no_matches_yet = False
            coords['ref_alignment_length'] += num
        else:
            sys.stderr.write("Don't recognize cigar character: %s, assuming it advances both query and reference, like a match or mismatch\n" % element[0])
            coords['read_alignment_length'] += num
            coords['ref_alignment_length'] += num
    return coords


if __name__ == "__main__": main()
