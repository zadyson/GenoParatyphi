# genoparatyphi
Detect AMR mutations in Salmonella Paratyphi A genomes based on VCF or BAM files (mapped to Paratyphi A AKU_12601 reference genome)


This script detects mutations associated with AMR in *Salmonella* Paratyphi A genomes including the QRDR (quinolone-resistance determining region) of genes *gyrA* and *parC* reported previouly ["Laboratory and Molecular Surveillance of Paediatric Typhoidal Salmonella in Nepal: Antimicrobial Resistance and Implications for Vaccine Policy", Britto et al 2018, PLoS NTDs], and mutations in *acrB* reported by Hooda et al 2019 ["Molecular mechanism of azithromycin resistance among typhoidal Salmonella strains in Bangladesh identified through passive pediatric surveillance", Hooda et al 2019].

The name genoparatyphi refers to the *Salmonella* Typhi genotyping framework described in this paper : ["An extended genotyping framework for Salmonella enterica serovar Typhi, the cause of human typhoid", Wong et al, 2016, Nature Communications, available at: https://github.com/katholt/genotyphi/.

Genoparatyphi does not currently type Paratyphi A genomes by lineage, nor does it detect acquired genes associated with AMR, for the later we recommend SRST2, available at: https://github.com/katholt/srst2/.



Inputs are BAM or VCF files (mapped to Paratyphi A AKU_12601 reference genome, accession number FM200053).

For short read data, we recommend using the raw alignments (BAM files) as input (via --bam), since VCFs generated by different pipelines may have filtered out important data at the AMR loci. Note you also need to supply the reference sequence that was used. If you don't have BAMs or you are really confident in the quality of your SNP calls in existing VCF files, you can provide your own VCFs (generated by read mapping) via --vcf.

For assemblies, we recommend using [ParSNP](http://harvest.readthedocs.org/) (version 1.0.1) to align genomes to the AKU_12601 reference. The resulting multi-sample VCF file(s) can be passed directly to this script via --vcf_parsnp.

Dependencies: Python 2.7.5+ ([SAMtools](http://samtools.sourceforge.net/) and [BCFtools](https://samtools.github.io/bcftools/) are also required if you are working from BAM files.  Genoparatyphi has been tested with versions 1.1, 1.2, and 1.3 of both SAMtools and bcftools, subsequently we advise using the same version of both of these dependencies together i.e. SAMtools v1.2 and bcftools v1.2).



[Basic Usage](https://github.com/zadyson/genoparatyphi/#basic-usage---own-bam-recommended-if-you-have-reads)

[Options](https://github.com/zadyson/genoparatyphi/#options)

[Outputs](https://github.com/zadyson/genoparatyphi/#outputs)

[Generating input BAMS from reads (with example)](https://github.com/zadyson/genoparatyphi/#generating-input-bams-from-reads)

[Generating input VCFs from assemblies (with example)](https://github.com/zadyson/genoparatyphi/#generating-input-vcfs-from-assemblies)

## Basic Usage - own BAM (recommended if you have reads)

Note the BAM files must be sorted (e.g. using samtools sort)

```
python genoparatyphi.py --mode bam --bam *.bam --ref FM200053.fasta --ref_id FM200053.1 --output genotypes.txt
```

## Basic Usage - own VCF

```
python genoparatyphi.py --mode vcf --vcf *.vcf --ref_id FM200053 --output genotypes.txt
```

## Basic Usage - assemblies aligned with ParSNP (recommended if you only have assembly data available and no reads)

```
python genoparatyphi.py --mode vcf_parsnp --vcf parsnp.vcf --output genotypes.txt
```

## Options

### Required options

```
-- mode            Run mode, either bam, vcf or vcf_parsnp
```

### Mode specific options

#### --mode bam

Requires [SAMtools](http://samtools.sourceforge.net/) and [BCFtools](https://samtools.github.io/bcftools/)

```
--bam              Path to one or more BAM files, generated by mapping reads to the Paratyphi A AKU_12601 (FM200053)
                   Note the SNP coordinates used here for genotyping are relative to Paratyphi A AKU_12601 (FM200053) 
                   so the input MUST be a BAM obtained via mapping to this reference sequence.

--ref              Reference sequence file used for mapping (fasta).

--ref_id           Name of the Paratyphi A AKU_12601 (FM200053) chromosome reference in your VCF file.
                   This is the entry in the first column (#CHROM) of the data part of the file.
                   This is necessary in case you have mapped to multiple replicons (e.g. chromosome and
                   plasmid) and all the results appear in the same VCF file.

--samtools_location     Specify the location of the folder containing the SAMtools installation if not standard/in path [optional]

--bcftools_location     Specify the location of the folder containing the bcftools installation if not standard/in path [optional]
```

#### --mode vcf

```
--vcf              Path to one or more VCF files, generated by mapping reads to the Paratyphi A AKU_12601 (FM200053)
                   Note the SNP coordinates used here for genotyping are relative to Paratyphi A AKU_12601 (FM200053)
                   so the input MUST be a VCF obtained via mapping to this reference sequence.
```

#### --mode vcf_parsnp

```
--vcf              Path to one or more VCF files generated by mapping assemblies to the Paratyphi A AKU_12601 (FM200053) reference genome
                   using ParSNP (--ref_id is optional, default value is '1').
```

### Other options

```
--phred                 Minimum phred quality to count a variant call vs Paratyphi A AKU_12601 (FM200053) as a true SNP (default 20)

--min_prop              Minimum proportion of reads required to call a SNP (default 0.1)

--output                Specify the location of the output text file (default genotypes_[timestamp].txt)

```

## Outputs

Output is to standard out, in tab delimited format.

TBC.


## Generating input BAMS from reads

We recommend using [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) to align reads to the Paratyphi A AKU_12601 (FM200053) reference, and [SAMtools](http://http://samtools.sourceforge.net/) to convert the *.sam file to the *.bam format.  The resulting bam file(s) can be passed directly to this script via --bam.  Please note the differences in the commands listed below when using SAMtools v1.1/1.2 vs. SAMtools v1.3.1 to generate bam files.

For example:

```
# Download Paratyphi A AKU_12601 (FM200053) and unzip the reference genome

wget -O ...

gunzip ...

# Separate the chromosome sequence from the plasmids with the emboss toolkit
seqretsplit ...

mv ...

# Replace the header line of the Paratyphi A AKU_12601 (FM200053) file with the reference id i.e. chage

...

# to
...

# Download reads for ...

wget ...
wget ...

# Use bowtie to map reads to the Paratyphi A AKU_12601 (FM200053) reference genome

For example, to align paired end reads to the Paratyphi A AKU_12601 (FM200053) reference genome sequence:

bowtie2-build ...

bowtie2 -p 2 -x ... -1 ... -2 ... -S ...
 
samtools view -bS ... > ...

samtools sort ... output

(or, 'samtools sort unsorted_output.bam > output.bam' for SAMtools v1.3.1 instead of SAMtools v1.2/1.1)

# Call ... genotypes from the resulting BAM(s)

python genoparatyphi.py --mode bam --bam output.bam --ref ... --ref_id ... --output genotypes_test.txt

```

#### Output

tBC

```
TBC

```

## Generating input VCFs from assemblies

We recommend using [ParSNP](http://harvest.readthedocs.org/) to align genomes to the Paratyphi A AKU_12601 (FM200053) reference. The resulting multi-sample VCF file(s) can be passed directly to this script via --vcf_parsnp.

For example:

```
# Download Paratyphi A AKU_12601 (FM200053) reference genome

wget ...

gunzip ...
mv ... ...

# Download two example Paratyphi A genomes for genotyping

wget ...
wget ...

gunzip ...
gunzip ...

mkdir genomes/
mv ... genomes/...
mv ... genomes/...

# Use ParSNP to generate variant calls (VCF) for these genomes against the Paratyphi A AKU_12601 (FM200053) reference sequence

parsnp -r ... -d genomes/ -o output

# Call ... from the resulting VCF

python genoparatyphi.py --mode vcf_parsnp --vcf output/parsnp.vcf --output genotypes_parsnptest.txt

```
#### Output

TBC

```
TBC

```

