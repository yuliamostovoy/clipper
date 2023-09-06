import sys
import os
import numpy
import scipy.cluster.hierarchy as hcluster
import argparse


def main():
    parser = argparse.ArgumentParser(description="Filter intermediate file derived from bigclipper_processbam.py\n")
    parser.add_argument('intermediate_file', help='Intermediate file produced by bigclipper_processbam.py')
    parser.add_argument('-d', '--min_dist', dest='min_dist', help="Minimum reference distance between split alignment positions of a single read (increase to find longer-range split-reads) (bp) [default=1]", type=int, default=1)
    parser.add_argument('-c', '--min_cluster_count', dest='min_cluster_count', help="Minimum number of reads in a cluster [default=5]", type=int, default=5)
    parser.add_argument('-s', '--cluster_dist', dest='cluster_dist', help="Maximum distance between supplementary alignment positions to cluster them together as a single \"supplementary alignment region\" [default=50]", type=int, default=50)
    parser.add_argument('-u', '--max_unique_breakends', dest='max_unique_breakends', help="Maximum number of supplementary alignment regions within a cluster (reduces false positives by filtering high-copy-number repeats) [default=10]", type=int, default=10)

    args = parser.parse_args()

    # preprocess intermediate file
    prefix = args.intermediate_file.rstrip('bed').rstrip('.')
    os.system("bedtools cluster -d %d -s -i %s | cut -f1,2,3,6,7,8,9,10,11 > %s_clustered" % (args.cluster_dist, args.intermediate_file, prefix))

    infile = open("%s_clustered" % prefix,'r')
    outfile = open("%s_d%d_c%d_s%d_u%d.vcf" % (prefix, args.min_dist, args.min_cluster_count, args.cluster_dist, args.max_unique_breakends), 'w')

    current_cluster = []
    current_cluster_ID = None

    # print header
    outfile.write("##fileformat=VCFv4.1\n")
    outfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    for line in infile:
        line = line.strip().split('\t')
        if line[-1] != current_cluster_ID:
            if len(current_cluster) >= args.min_cluster_count:
                process_cluster(current_cluster, args.min_dist, args.max_unique_breakends, args.cluster_dist, outfile)
            current_cluster = [line]
            current_cluster_ID = line[-1]
        else:
            current_cluster.append(line)

    if len(current_cluster) >= args.min_cluster_count:  # process last one
        process_cluster(current_cluster, args.min_dist, args.max_unique_breakends, args.cluster_dist, outfile)

def process_cluster(cluster, min_dist, max_unique_breakends, thresh, outfile):
    ref = cluster[0][0]
    pos = int(cluster[0][1])
    ID = int(cluster[0][8])

    max_mapq = 0

    other_breakends = []
    for line in cluster:
        ori1 = line[3]
        other_breakends.append(line[7])
        max_mapq = max(max_mapq, int(line[4]))
    unique_breakends = list(set(other_breakends))

    breakend_info = []  # list of dicts with chr, pos, ori, count
    for breakend in unique_breakends:
        breakend_ori = breakend[-1]
        count = other_breakends.count(breakend)
        breakend = breakend[:-1].split(':')  # trim off the orientation character, convert to [chr, pos]
        breakend_ref = breakend[0]
        breakend_pos = int(breakend[1])
        if breakend_ref != ref or abs(pos-breakend_pos)>=min_dist:
            breakend_info.append({'ref':breakend_ref, 'pos':breakend_pos, 'ori':breakend_ori, 'count':count})
    if not breakend_info:  # none of the other breakends passed the min_dist filter
        return

    # cluster the other breakends
    cluster_info = []  # like breakend_info but clustered by position
    for breakend_ref in set([x['ref'] for x in breakend_info]):
        positions = {'+': [], '-': []}
        for breakend in filter(lambda x: x['ref'] == breakend_ref, breakend_info):
            for i in range(breakend['count']):
                positions[breakend['ori']].append(breakend['pos'])
        for ori in positions:
            if len(positions[ori]) > 1:
                ndata = [[x, x] for x in positions[ori]]
                data = numpy.array(ndata)
                clusters = hcluster.fclusterdata(data, thresh, criterion="distance")
                for i in range(len(set(clusters))):
                    index = list(filter(lambda j: clusters[j] == i + 1, range(len(clusters))))
                    pos_in_cluster = []
                    for j in index:
                        pos_in_cluster.append(positions[ori][j])
                    cluster_info.append({'ref': breakend_ref, 'pos1': min(pos_in_cluster), 'pos2': max(pos_in_cluster), 'ori': ori, 'count': len(pos_in_cluster)})
            elif len(positions[ori]) == 1:
                cluster_info.append({'ref': breakend_ref, 'pos1': positions[ori][0], 'pos2': positions[ori][0], 'ori': ori, 'count': 1})

    if len(cluster_info) > max_unique_breakends:
        return

    if not cluster_info:
        return
    cluster_info = sorted(cluster_info, key=lambda x: x['count'], reverse=True)  # sort by frequency

    print_cluster(ref, pos, ID, max_mapq, cluster_info, outfile)


def print_cluster(ref, pos, ID, qual, cluster_info, outfile):
    # VCF format is 1-based!
    outfile.write("%s\t%d\t%d\t.\t.\t%d\tPASS\tALNS=" % (ref, pos+1, ID, qual))
    for hit in cluster_info:
        outfile.write("%s:%d-%d(%d,%s);" % (hit['ref'], hit['pos1']+1, hit['pos2']+1, hit['count'], hit['ori']))
    outfile.write("\n")


if __name__ == "__main__": main()
