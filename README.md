<p align="center">
  <img src="images/GENEX_logo.png" alt="Left Logo" height="400"/>
  &nbsp;&nbsp;&nbsp;
  <img src="images/aeDNA_logo.png" alt="Right Logo" height="400"/>
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

Please note that in this gitub reporsitory, we provide a small subset of microbial pseudo-reads for demonstration purposes, the full dataset is available at the SciLifeLab Figshare https://doi.org/10.17044/scilifelab.28491956.

Questions regarding the dataset should be sent to nikolay.oskolkov@scilifelab.se

## Quick start

Please clone this repository and install the workflow tools as follows:

    git clone https://github.com/NikolayOskolkov/MCWorkflow
    cd MCWorkflow
    conda env create -f environment.yaml
    conda activate MCWorkflow

Then you can run the workflow as:

    ./micr_cont_detect.sh GCF_002220235.fna.gz data GTDB 4 \
    GTDB_sliced_seqs_sliding_window.fna.gz 10

Here, `GCF_002220235.fna.gz` is the eukaryotic reference to be screened for microbial-like sequeneces, `data` is the directory containing the eukaryotic reference, `GTDB` is the type of pseudo-reads to be used for detecting exogenous regions in the eukaryotic reference (can be `GTGB`, `RefSeq` or `human`), `4` is the number of available threads in your computational environment, `GTDB_sliced_seqs_sliding_window.fna.gz` is the pre-computed pseudo-reads (small subset is provided in this github repository, the full datasets can be downloaded from the SciLifeLab Figshare https://doi.org/10.17044/scilifelab.28491956), and `10` is the number of allowed Bowtie2 multi-mappers.


Please also read the very detailed `vignette.html` and follow the preparation steps described there. The vignette `vignette.html` walks you through the explanations of the workflow parameters and interpretation of the output files.



## Nextflow implementation

Alternatively, you can specify the workflow input files and parameters in the `nextflow.config` and run it using Nextflow:

    nextflow run main.nf

The Nextflow implementation is preferred for scalability and reproducibility purposes. Please place your reference genomes (fasta-files) to be screened for exogenous regions in the `data` folder. An example of the config-file, `nextflow.config`, can look like this:

    params {
        input_dir = "data"                                             // folder with multiple reference genomes (fasta-files)
        type_of_pseudo_reads = "GTDB"                                  // type of pseudo-reads to be used for screening the input reference genome, can be "GTDB", "RefSeq" or "human"
        threads = 4                                                    // number of available threads
        input_pseudo_reads = "GTDB_sliced_seqs_sliding_window.fna.gz"  // name of pre-computed file with pseudo-reads, can be "GTDB_sliced_seqs_sliding_window.fna.gz", "RefSeq_sliced_seqs_sliding_window.fna.gz" or "human_sliced_seqs_sliding_window.fna.gz"
        n_allowed_multimappers = 10                                    // number of multi-mapping pseudo-reads allowed by Bowtie2, do not change this default number unless you know what you are doing
    }

Please modify it to adjust for the number of available threads in your computational environment and the type of analysis, i.e. detecting microbial-like or human-like sequeneces in the reference genome, you would like to perform.
