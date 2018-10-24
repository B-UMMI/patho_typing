# patho_typing

*In silico pathogenic typing directly from raw Illumina reads*  

---

* [Rational](#rational)
* [Input requirements](#input-requirements)
* [Dependencies](#dependencies)
  * [Install dependencies](#install-dependencies)
* [Install patho_typing](#install-patho_typing)
* [Usage](#usage)
* [Outputs](#outputs)
* [Citation](#citation)
* [Contact](#contact)

## Rational
<html>
 <div align="right">
  <a href="#seq_typing">Back to top</a><br>
 </div>
</html>

**patho_typing** is a tool for _in silico_ pathogenic typing sample's reads through a read mapping approach using a set of reference sequences and defined rules for sequences presence/absence.  
Sample's reads are mapped to the given reference sequences using [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), parsed with [Samtools](http://www.htslib.org/) and analysed via [ReMatCh](https://github.com/B-UMMI/ReMatCh). Based on the length of the sequence covered, it's depth of coverage and sequence nucleotide identity, **patho_typing** scores those for presence or absence, following defined thresholds. According to the combination of sequences present, a pathotype is returned following a set of rules for sequences presence/absence. Some of the sequences can be either present or absent.
Reference sequences definition and presence/absence rules delineation are required pathotyping classification.

## Input requirements
<html>
 <div align="right">
  <a href="#seq_typing">Back to top</a><br>
 </div>
</html>

 - Fastq file

## Dependencies
<html>
 <div align="right">
  <a href="#seq_typing">Back to top</a><br>
 </div>
</html>

* Python 3
* [ReMatCh](https://github.com/B-UMMI/ReMatCh)

### Install dependencies
<html>
 <div align="right">
  <a href="#seq_typing">Back to top</a><br>
 </div>
</html>


ReMatCh:
```bash
git clone https://github.com/B-UMMI/ReMatCh.git
cd ReMatCh
python3 setup.py install
```
*__NOTE__*:  
If you don't have permission for global system installation, try the following _install_ command instead:  
`python3 setup.py install --user`

## Install patho_typing
<html>
 <div align="right">
  <a href="#seq_typing">Back to top</a><br>
 </div>
</html>

```bash
git clone https://github.com/B-UMMI/patho_typing.git
cd patho_typing
python3 setup.py install
```
*__NOTE__*:  
If you don't have permission for global system installation, try the following _install_ command instead:  
`python3 setup.py install --user`

## Usage
<html>
 <div align="right">
  <a href="#seq_typing">Back to top</a><br>
 </div>
</html>

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


## Outputs
<html>
 <div align="right">
  <a href="#seq_typing">Back to top</a><br>
 </div>
</html>

**run.*.log**  
*ReMatCh* running log file.  

**patho_typing.report.txt**  
List with possible pathotypes found.

***ReMatCh* folder**  
For each *ReMatCh* run (typing and trueCoverage, the latter when required), the *rematchModule_report.txt* is kept.  
 - *rematchModule_report.txt* - Report file containing gene information: 1) gene name, 2) percentage of target gene sequence covered, 3) Mean target gene coverage depth of present positions, 4) percentage of target gene sequence with lower coverage depth, 5) number of positions in target gene sequence containing multiple alleles, 6) percentage identity of target gene sequence covered. The general sample information will also be stored: number of absent genes, number of genes with multiple alleles among the genes present and the mean sample coverage depth (only considering the genes present).

## Citation
<html>
 <div align="right">
  <a href="#seq_typing">Back to top</a><br>
 </div>
</html>

MP Machado, J Halkilahti, M Pinto, JP Gomes, M Ramirez, M Rossi, JA Carrico. _patho_typing_ **GitHub** https://github.com/B-UMMI/patho_typing

## Contact
<html>
 <div align="right">
  <a href="#seq_typing">Back to top</a><br>
 </div>
</html>

Miguel Machado  
<mpmachado@medicina.ulisboa.pt>
