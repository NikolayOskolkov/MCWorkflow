# Microbial Contamination Workflow

This is a vignette demonstrating the computational workflow for detecting coordinates of microbial-like sequences in eukaryotic reference genomes. The workflow accepts a reference genome in FASTA-format and outputs coordinates of microbial-like regions in BED-format. The workflow builds a Bowtie2 index of the eukaryotic reference genome and aligns pre-computed microbial GTDB v.214 (https://gtdb.ecogenomic.org/) pseudo-reads to the reference, then custom scripts are used for detection of the positions of covered regions and quantification of most abundant microbial contaminants.

The workflow was developed by Nikolay Oskolkov, Lund University, Sweden, within the NBIS SciLifeLab long-term support project, PI Tom van der Valk, Centre for Palaeogenetics, Stockholm, Sweden.

If you use the workflow for your research, please cite our manuscript:

    Nikolay Oskolkov, Chenyu Jin, Samantha López Clinton, Flore Wijnands, Ernst Johnson, 
    Benjamin Guinet, Verena Kutschera, Cormac Kinsella, Peter D. Heintzman and Tom van der Valk, 
    Disinfecting eukaryotic reference genomes to improve taxonomic inference from environmental 
    ancient metagenomic data, https://www.biorxiv.org/content/10.1101/2025.03.19.644176v1, 
    https://doi.org/10.1101/2025.03.19.644176

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

Alternatively, you can specify the workflow input files and parameters in the `nextflow.config` and run it using Nextflow:

    nextflow run main.nf

Please also read the very detailed `vignette.html` and follow the preparation steps described there. The vignette `vignette.html` walks you through the explanations of the workflow parameters and interpretation of the output files.
