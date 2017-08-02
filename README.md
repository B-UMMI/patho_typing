patho_typing
============
*In silico pathogenic typing pathogenic directly from raw Illumina reads*  

<https://github.com/B-UMMI/patho_typing>



Requirements
------------

 - Fastq file



Dependencies
------------

 - *ReMatCh* >= v3.2 (<https://github.com/B-UMMI/ReMatCh>)



Installation
------------
    git clone https://github.com/B-UMMI/patho_typing.git



Usage
-----
    usage: patho_typing.py [-h] [--version] -f /path/to/input/file.fq.gz
                           [/path/to/input/file.fq.gz ...] -s Yersinia
                           enterocolitica [-o /path/to/output/directory/] [-j N]
                           [--trueCoverage] [--noCheckPoint] [--minGeneCoverage N]
                           [--minGeneIdentity N] [--minGeneDepth N]
                           [--doNotRemoveConsensus] [--debug]

    In silico pathogenic typing directly from raw Illumina reads

    optional arguments:
      -h, --help            show this help message and exit
      --version             Version information

    Required options:
      -f /path/to/input/file.fq.gz [/path/to/input/file.fq.gz ...], --fastq /path/to/input/file.fq.gz [/path/to/input/file.fq.gz ...]
                            Path to single OR paired-end fastq files. If two files
                            are passed, they will be assumed as being the paired
                            fastq files (default: None)
      -s Yersinia enterocolitica, --species Yersinia enterocolitica
                            Species name (default: None)

    General facultative options:
      -o /path/to/output/directory/, --outdir /path/to/output/directory/
                            Path to the directory where the information will be
                            stored (default: .)
      -j N, --threads N     Number of threads to use (default: 1)
      --trueCoverage        Assess true coverage before continue typing (default:
                            False)
      --noCheckPoint        Ignore the true coverage checking point (default:
                            False)
      --minGeneCoverage N   Minimum typing percentage of target reference gene
                            sequence covered to consider a gene to be present
                            (value between [0, 100]) (default: None)
      --minGeneIdentity N   Minimum typing percentage of identity of reference
                            gene sequence covered to consider a gene to be present
                            (value between [0, 100]). One INDEL will be considered
                            as one difference (default: None)
      --minGeneDepth N      Minimum typing gene average coverage depth of present
                            positions to consider a gene to be present (default
                            15, or 1/3 of average sample coverage assessed by true
                            coverage analysis) (default: None)
      --doNotRemoveConsensus
                            Do not remove ReMatCh consensus sequences (default:
                            False)
      --debug               DeBug Mode: do not remove temporary files (default:
                            False)



Outputs
-------
**run.*.log**  
*ReMatCh* running log file.  

**patho_typing.report.txt**  
List with possible pathotypes found.

***ReMatCh* folder**  
For each *ReMatCh* run (typing and trueCoverage, the latter when required), the *rematchModule_report.txt* is kept.  
 - *rematchModule_report.txt* - Report file containing gene information: 1) gene name, 2) percentage of target gene sequence covered, 3) Mean target gene coverage depth of present positions, 4) percentage of target gene sequence with lower coverage depth, 5) number of positions in target gene sequence containing multiple alleles, 6) percentage identity of target gene sequence covered. The general sample information will also be stored: number of absent genes, number of genes with multiple alleles among the genes present and the mean sample coverage depth (only considering the genes present).



Contact
-------
Miguel Machado  
<mpmachado@medicina.ulisboa.pt>
