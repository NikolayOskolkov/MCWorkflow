<p align="center">
  <img src="images/GENEX_logo.png" alt="Left Logo" height="320"/>
  &nbsp;&nbsp;&nbsp;
  <img src="images/aeDNA_logo.png" alt="Right Logo" height="320"/>
</p>

# GENome EXogenous (GENEX) sequence detection

This is a computational workflow for detecting coordinates of microbial-like or human-like sequences in eukaryotic and procaryotic reference genomes. The workflow accepts a reference genome in FASTA-format and outputs coordinates of microbial-like (human-like) regions in BED-format. The workflow builds a Bowtie2 index of the reference genome and aligns pre-computed microbial (GTDB v.214 or NCBI RefSeq release 213) or human (hg38) pseudo-reads to the reference, then custom scripts are used for detection of the positions of covered regions and quantification of most abundant microbial species, the latter is only when screening for microbial-like sequences in eukaryotic referenes.

The workflow was developed by Nikolay Oskolkov, Lund University, Sweden, within the NBIS SciLifeLab long-term support project, PI Tom van der Valk, Centre for Palaeogenetics, Stockholm, Sweden.

If you use the workflow for your research, please cite our manuscript:

    Nikolay Oskolkov, Chenyu Jin, Samantha López Clinton, Benjamin Guinet, Flore Wijnands, 
    Ernst Johnson, Verena E Kutschera, Cormac M Kinsella, Peter D Heintzman, Tom van der Valk, 
    Improving taxonomic inference from ancient environmental metagenomes by masking 
    microbial-like regions in reference genomes, 
    GigaScience, Volume 14, 2025, giaf108, https://doi.org/10.1093/gigascience/giaf108

Questions regarding the dataset should be sent to nikolay.oskolkov@scilifelab.se
Question regarding the nextflow workflow refers to Chenyu.Jin(amend.jin@gmail.com).

## Quick start

Required programs\
`nextflow`\
`conda`\
`singularity` OR `apptainer`\
`seqkit`

(1) Please clone this repository and install the workflow tools as follows:

    git clone https://github.com/NikolayOskolkov/MCWorkflow
    cd MCWorkflow

(2) Now we need to creat a directory of which the path is given to the parameter `pseudo_reads_file_dir` in `nextflow.config`. Please note that in this gitub reporsitory, we provide a small subset of microbial pseudo-reads for demonstration purposes, the full dataset is available at the SciLifeLab Figshare `https://doi.org/10.17044/scilifelab.28491956`.

After downloading the needed fna.gz (e.g., GTDB_sliced_seqs_sliding_window.fna.gz) in the `pseudo_reads_file_dir` directory, you can then run which might take > 6 hours to obtain the subsetted database:\
`seqkit split -s 10000000 GTDB_sliced_seqs_sliding_window.fna.gz`

(3) All inputs are specified in `nextflow.config`. To `nextflow run`, you first need to modify:\
   `input_dir`: the path to the directory with all fasta files (gzipped or not)\
   `type_of_pseudo_reads`: "GTDB" # or "RefSeq", "human" depends on which database you want to use to mask\
    `pseudo_reads_file_dir`: where it contains all the subsets of sliced GTDB or other databases\
    `n_allowed_multimappers`: the number of allowed multimapper. Based on the test done in the paper, 10 is recommended.\
    `output_dir`: the path to the directory where you want all the outputs \
    `work_dir`: where your MCWorkflow directory is\
    `fna2name`: the contig to species name correspondance file. It's GTDB_fna2name.txt for GTDB.\

(4) Then you can run the workflow as:

    nextflow run main.nf -profile apptainer,conda -c nextflow.config,dardel.config -resume -with-trace
You can use  `dardel.config` if you want to submit jobs on SLURM. Or use the config file of your cluser.\
`-with-trace` is used if you are interested to know the resources (memory and time) used by each processes within the workflow. In the workflow the cpu and ram are pre-specified and will retry with higher number of cpu and alloted time if there is out of memory issues.
