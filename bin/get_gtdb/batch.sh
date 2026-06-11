#!/bin/bash -l
#SBATCH -A project_code
#SBATCH -t 5:00:00
#SBATCH -N 1
#SBATCH -c 3
#SBATCH -p pelle
#SBATCH -J download
#SBATCH --mail-type=BEGIN,END

cp download_figshare.sh api.py $TMPDIR

cd $TMPDIR

bash download_figshare.sh

mv GTDB_sliced_seqs_sliding_window.fna.gz /location/for/files
