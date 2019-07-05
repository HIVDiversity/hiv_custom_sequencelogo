# hiv_custom_sequencelogo
This script is designed to take a full gene (or partial gene) amino acid sequence alignment of HIV sequences (must include HXB2) and create a custom sequence logo using a list of sites provided in a csv file


# Requirements
python v3.6>

conda install -c anaconda pandas 

conda install -c bioconda weblogo

* A sequence alignment containing HXB2

* A csv file listing the sites to include in a column

    * This file must contain the column heading "sites"

