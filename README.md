<p align="center">
  <img src="images/GENEX_logo.png" alt="Left Logo" height="320"/>
  &nbsp;&nbsp;&nbsp;
  <img src="images/aeDNA_logo.png" alt="Right Logo" height="320"/>
</p>

# GENome EXogenous (GENEX) sequence detection

This is a computational workflow for detecting coordinates of microbial-like or human-like sequences in eukaryotic and procaryotic reference genomes. The workflow accepts a reference genome in FASTA-format and outputs coordinates of microbial-like (human-like) regions in BED-format. The workflow builds a Bowtie2 index of the reference genome and aligns pre-computed microbial (GTDB v.214 or NCBI RefSeq release 213) or human (hg38) pseudo-reads to the reference, then custom scripts are used for detection of the positions of covered regions and quantification of most abundant microbial species, the latter is only when screening for microbial-like sequences in eukaryotic referenes.

The workflow was developed by Nikolay Oskolkov, Lund University, Sweden, within the NBIS SciLifeLab long-term support project, and PhD student Chenyu.Jin, PI Tom van der Valk, Centre for Palaeogenetics, Stockholm, Sweden.

If you use the workflow for your research, please cite our manuscript:

    Nikolay Oskolkov, Chenyu Jin, Samantha López Clinton, Benjamin Guinet, Flore Wijnands, 
    Ernst Johnson, Verena E Kutschera, Cormac M Kinsella, Peter D Heintzman, Tom van der Valk, 
    Improving taxonomic inference from ancient environmental metagenomes by masking 
    microbial-like regions in reference genomes, 
    GigaScience, Volume 14, 2025, giaf108, https://doi.org/10.1093/gigascience/giaf108

Questions regarding the dataset should be sent to nikolay.oskolkov@scilifelab.se
Question regarding the nextflow workflow refers to Chenyu.Jin (amend.jin@gmail.com)

## Quick start

Required programs\
`nextflow`\
`conda`\
`singularity` OR `apptainer`\
`seqkit`

1) Please clone this repository and install the workflow tools as follows:

    git clone https://github.com/NikolayOskolkov/MCWorkflow
    cd MCWorkflow

2) Now we need to creat a directory of which the path is given to the parameter `pseudo_reads_file_dir` in `nextflow.config`. Please note that in this gitub reporsitory, we provide a small subset of microbial pseudo-reads for demonstration purposes, the full dataset is available at the SciLifeLab Figshare `https://doi.org/10.17044/scilifelab.28491956`.

After downloading the needed fna.gz (e.g., GTDB_sliced_seqs_sliding_window.fna.gz) in the `pseudo_reads_file_dir` directory, you can then run which might take > 6 hours to obtain the subsetted database:\
`seqkit split -s 10000000 GTDB_sliced_seqs_sliding_window.fna.gz`

3) All inputs are specified in `nextflow.config`. To `nextflow run`, you first need to modify:\
   `input_dir`: the path to the directory with all fasta files (gzipped or not)\
   `type_of_pseudo_reads`: "GTDB" # or "RefSeq", "human" depends on which database you want to use to mask\
    `pseudo_reads_file_dir`: where it contains all the subsets of sliced GTDB or other databases\
    `n_allowed_multimappers`: the number of allowed multimapper. Based on the test done in the paper, 10 is recommended.\
    `output_dir`: the path to the directory where you want all the outputs \
    `work_dir`: where your MCWorkflow directory is\
    `fna2name`: the contig to species name correspondance file. It's GTDB_fna2name.txt for GTDB.\

4) Then you can run the workflow as:

    nextflow run main.nf -profile apptainer,conda -c nextflow.config,dardel.config -resume -with-trace
You can use  `dardel.config` if you want to submit jobs on SLURM. Or use the config file of your cluser.\
`-with-trace` is used if you are interested to know the resources (memory and time) used by each processes within the workflow. In the workflow the cpu and ram are pre-specified and will retry with higher number of cpu and alloted time if there is out of memory issues.

## Interepreting results

```
-rw-r--r-- 1 A B     55 Feb 12 14:07 mito.subset_GTDB_contigs_abund_sorted_GTDB_mito.subset.txt
-rw-r--r-- 1 A B     82 Feb 12 14:07 mito.subset_GTDB_contigs_boc_sorted_GTDB_mito.subset.txt
-rw-r--r-- 1 A B    270 Feb 12 14:07 mito.subset_GTDB_coords_micr_like_regions_mito.subset.bed
-rw-r--r-- 1 A B   1152 Feb 12 14:07 mito.subset_GTDB_coords_micr_like_regions_mito.subset.txt
-rw-r--r-- 1 A B    828 Feb 12 14:07 mito.subset_GTDB_microbes_abundant_GTDB_mito.subset.txt
-rw-r--r-- 1 A B 234066 Feb 12 14:22 mito.subset.masked.fna
```

### The major file is `*.bed` which contains the regions to mask, and the `*.masked.fna` file as the masked reference genome.


Other outputs to understand which microbes/organisms/contigs matches the masked region on reference genome:

`*coords_micr_contam_*.txt`:
(1) name (id) of the eukaryotic reference genome profiles, 
(2) contig / scaffold / chromosome id within the eukaryotic reference genome containing microvbial-like region,
(3) start coordinate of the detected microbial-like region,
(4) end coordinate of the detected microbial-like region,
(5) genomic length of the microbial-like region,
(6) total number of reads aligned to the detected microbial-like region,
(7) average number of reads supporting each position within the detected microbial-like region,
(8) the next five columns represent top abundant microbial species for each detected microbial-like region (the number of reads is reported for each of the top abundant microbes); if fewer than five unique microbes are dicovered within the microbial-like region, the rest of the columns contain recods "NA_reads_NA".

```
*coords_micr_contam*.txt:
ORGANISM	CONTIG	START	END	LENGTH	N_READS	AVERAGE_DEPTH	MICR1	MICR2	MICR3
GCF_002220235.fna.gz	NC_023997.1	231657	231716	59	1	1	1_reads_UBA796_sp002707085	NA_reads_NA	NA_reads_NA
GCF_002220235.fna.gz	NC_023997.1	231807	232156	349	26	4	26_reads_UBA796_sp002707085	NA_reads_NA	NA_reads_NA
GCF_002220235.fna.gz	NC_023997.1	232752	232891	139	9	4	9_reads_UBA796_sp002707085	NA_reads_NA	NA_reads_NA
GCF_002220235.fna.gz	NC_023997.1	232905	232964	59	1	1	1_reads_UBA796_sp002707085	NA_reads_NA	NA_reads_NA
GCF_002220235.fna.gz	NC_023997.1	233015	233104	89	4	3	4_reads_UBA796_sp002707085	NA_reads_NA	NA_reads_NA
GCF_002220235.fna.gz	NC_023997.1	233145	233369	224	13	3	13_reads_UBA796_sp002707085	NA_reads_NA	NA_reads_NA
GCF_002220235.fna.gz	NC_023997.1	233373	233377	4	1	1	1_reads_UBA796_sp002707085	NA_reads_NA	NA_reads_NA
GCF_002220235.fna.gz	NC_023997.1	233395	233494	99	3	2	3_reads_UBA796_sp002707085	NA_reads_NA	NA_reads_NA
GCF_002220235.fna.gz	NC_023997.1	233515	233694	179	11	4	11_reads_UBA796_sp002707085	NA_reads_NA	NA_reads_NA
Further, congigs within the eukaryotic reference genome sorted by the total number of aligned microbial pseuso-reads and contig-wide breadth of coverage are reported in contigs_abund_sorted_GTDB_GCF_002220235.fna.gz.tx and contigs_boc_sorted_GTDB_GCF_002220235.fna.gz.txt files, respectively.
```

```
head *contigs_abund_sorted*.txt
ORGANISM	CONTIG	N_READS
GCF_002220235.fna.gz	NC_023997.1	347
GCF_002220235.fna.gz	NC_024004.1	274
GCF_002220235.fna.gz	NC_024008.1	221
GCF_002220235.fna.gz	NC_023992.1	193
GCF_002220235.fna.gz	NC_024006.1	186
GCF_002220235.fna.gz	NC_024001.1	175
GCF_002220235.fna.gz	NC_024000.1	123
GCF_002220235.fna.gz	NC_024003.1	59
```

```
head *contigs_boc_sorted_GTDB_*.txt
ORGANISM	CONTIG	N_READS	LENGTH	BOC
GCF_002220235.fna.gz	NC_023997.1	347	712459	0.00871348
GCF_002220235.fna.gz	NC_023992.1	193	465570	0.00777327
GCF_002220235.fna.gz	NC_024004.1	274	1019276	0.00525471
GCF_002220235.fna.gz	NC_024008.1	221	1352724	0.00400451
GCF_002220235.fna.gz	NC_024001.1	175	937610	0.00396433
GCF_002220235.fna.gz	NC_024006.1	186	1091008	0.00389915
GCF_002220235.fna.gz	NC_024000.1	123	895536	0.00294349
GCF_002220235.fna.gz	NC_024003.1	59	989707	0.00161462
```

The N_READS column in both files reports the number of microbial pseudo-reads aligned to each contig / scaffold / chromosome within the eukaryotic reference genome. In addition, the LENGTH and BOC columns in the "boc_sorted" file report the length of each contig / scaffold / chromosome in the eukaryotic reference and the breadth of coverage (BOC - fraction of reference nucleotides covered at least once), respectively. Finally, microbes_abundant_GTDB_GCF_002220235.fna.gz.txt file contains GTDB microbial species (their reference genome ids) sorted by their total abundance within the whole eukaryotic reference genome. Note that this file will be absent when screening a reference genome for human contamination since all aligned pseudo-reads belong to human in this case.

```
head microbes_abundant_GTDB_GCF_002220235.fna.gz.txt
  1578 5oFfr3yp0G.fna
```

In this particular example, we had only one GTDB microbial species UBA796 sp002707085 (with the reference 5oFfr3yp0G.fna) belonging to the phylum Myxococcota in the GTDB databse. All the 1578 micobial pseudo-reads aligned to the green algae reference genome belonged to that single microbe. The bam-file containng the alignments is PseudoReads_aligned_to_GCF_002220235.fna.gz.bam, let us quickly check that there are 1578 aligned reads in total.
