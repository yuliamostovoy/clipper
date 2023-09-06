# BigClipper

A tool to identify pile-ups of clipped reads whose clipped portions align elsewhere in the genome, useful for identifying long-range rearrangements that may be missed by other SV callers. Written for long-read genome sequencing data, e.g. from PacBio and Oxford Nanopore.

## Requirements
Bedtools

Python >= 3

**Python libraries:**
  * pysam
  * numpy
  * scipy
  * argparse

## Installation
To create a conda environment for running BigClipper locally:

```
git clone https://github.com/yuliamostovoy/bigclipper.git
cd bigclipper
conda env create -f environment.yml
conda activate bigclipper
```

Alternatively, a docker image and WDL workflow are available for running BigClipper on the cloud.

## Usage
BigClipper runs in two phases. The first -- which is slower and only needs to be run once -- reads a BAM file, detects soft-clipped reads with supplementary alignments, and outputs their information to an intermediate file. The second phase parses the intermediate file with a set of parameters and outputs pile-ups of clipped reads that pass the filtering parameters. The second phase can be run multiple times with different parameters.

Phase 1: process the BAM file
```
python bigclipper/scripts/bigclipper_processbam.py [BAM] [output_prefix]
```
Phase 1 output: `[output_prefix]_intermediate.bed`

Phase 2: process the intermediate file and output clusters
```
python bigclipper/scripts/bigclipper_processbam.py [-h] [-d MIN_DIST] [-c MIN_CLUSTER_COUNT] [-u MAX_UNIQUE_BREAKENDS] [-s CLUSTER_DIST] intermediate_file

positional arguments:
  intermediate_file           Intermediate file produced by bigclipper_processbam.py

optional arguments:
  -h, --help            show help message and exit
  -d MIN_DIST, --min_dist MIN_DIST
                        Minimum reference distance between split alignment positions of a single read (increase to find longer-range split-reads) (bp) [default=1]
  -c MIN_CLUSTER_COUNT, --min_cluster_count MIN_CLUSTER_COUNT
                        Minimum number of reads in a cluster [default=5]
  -s CLUSTER_DIST, --cluster_dist CLUSTER_DIST
                        Maximum distance between supplementary alignment positions to cluster them together as a single "supplementary alignment region" [default=50]
  -u MAX_UNIQUE_BREAKENDS, --max_unique_breakends MAX_UNIQUE_BREAKENDS
                        Maximum number of supplementary alignment regions within a cluster (reduces false positives by filtering high-copy-number repeats) [default=10]

```

## Output
VCF file containing filtered clusters of split reads, including the locations of supplementary alignments in the cluster
