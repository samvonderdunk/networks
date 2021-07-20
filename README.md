# networks
Gene Regulatory Network visualization
[Published 20 July 2021]

##################################################################

>>> Read my blog on how to visualize GRNs using igraph in the script networks.R: https://dunkandcover.blogspot.com/2021/07/custom-visualization-of-gene-regulatory.html

##################################################################

The script is called networks.R
Example files are provided as MA_B4_genes.net, MA_B4_genes.annot, MA_B4_genome.net, MA_B4_genome.annot.

Make the different network representation as follows (i.e. choose the right option, and provide the right data files).

1. Circular genome:
networks.R ~/location-to-current-dir-where-we-find-the-data-files/ MA_B4_genome.net MA_B4_genome.annot genome

2. Linear genome:
networks.R ~/location-to-current-dir-where-we-find-the-data-files/ MA_B4_genome.net MA_B4_genome.annot linear

3. Gene network:
networks.R ~/location-to-current-dir-where-we-find-the-data-files/ MA_B4_genes.net MA_B4_genes.annot genes

##################################################################
