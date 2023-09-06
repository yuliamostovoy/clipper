version 1.0

######################################################################################
## A pipeline for running the clipper script
######################################################################################

workflow RunClipper {
    input {
        File aligned_bam
        String prefix
        Int? min_dist
        Int? min_cluster_count
        Int? max_unique_breakends
        Int? max_cluster_dist
    }

    parameter_meta {
        aligned_bam:            "aligned bam"
        prefix:                 "e.g. sample name"
        min_dist:               "Minimum distance between where read splits align on the reference (bp) (default=1)"
        min_cluster_count:      "Minimum number of reads in a cluster (default=5)"
        max_unique_breakends:   "Maximum number of unique breakends for a cluster (default=10)"
        max_cluster_dist:       "Max distance between supplementary breakends to cluster them together (default=50)"
    }

    call ProcessBam { input: aligned_bam = aligned_bam, prefix=prefix }

    call GetClusters {
      input:
        intermediate_file = ProcessBam.intermediate_file,
        prefix = prefix,
        min_dist = select_first([min_dist, 1]),
        min_cluster_count = select_first([min_cluster_count, 5]),
        max_unique_breakends = select_first([max_unique_breakends, 10]),
        max_cluster_dist = select_first([max_cluster_dist, 50])
    }

    output {
        File clipped_vcf = GetClusters.clustervcf
    }
}

task ProcessBam {
    input {
        File aligned_bam
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(aligned_bam, "GB"))+20

    command <<<
        set -euxo pipefail

        python /bigclipper_processbam.py ~{aligned_bam} ~{prefix}
    >>>

    output {
        File intermediate_file = "~{prefix}.bed"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        0,
        docker:             "quay.io/ymostovoy/clipper:1.1"
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task GetClusters {
    input {
        File intermediate_file
        String prefix
        Int min_dist
        Int min_cluster_count
        Int max_unique_breakends
        Int max_cluster_dist
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = ceil(size(clusterfile, "GB"))+50

    command <<<
        python /bigclipper_getclusters.py ~{intermediate_file} -d ~{min_dist} -c ~{min_cluster_count} -s ~{max_cluster_dist} -u ~{max_unique_breakends}
    >>>

    output { File clustervcf = "~{prefix}_clipped_reads_d~{min_dist}_c~{min_cluster_count}_s~{max_cluster_dist}_u~{max_unique_breakends}.vcf" }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        0,
        docker:             "quay.io/ymostovoy/clipper:1.1"
    }

     RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
     runtime {
         cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
         memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
         disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
         bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
         preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
         maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
         docker:                 select_first([runtime_attr.docker,            default_attr.docker])
     }
}

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
    String? docker
}
