Overview
--------
`scurgen` is a command line tool for plotting genomic data using scpace-filling
curves such as the Hilbert Curve.

.. WARNING::
   This project is in its infancy.  Documentation is sparse, functionality is 
   evolving, and installation/dependencies are not documented.  

   
The plot module
---------------

This module allows you to make a basic plot of genomic data using a Hilbert 
curve.::

	python scurgen.py plot  --chrom chr10 \
	                        --color Greens \
	                        --format png \
	                        --dim 128 \
	                        data/cpg-islands.hg19.chr10.bed

yields a PNG image named cpg-islands.hg19.chr10.bed.png:

	.. image:: https://raw.github.com/arq5x/scurgen/master/cpg-islands.hg19.chr10.bed.png

And let's do the same with 5' UTRs.::

	python scurgen.py plot  --chrom chr10 \
	                        --color Greens \
	                        --format png \
	                        --dim 128 \
	                        data/refseq.chr10.exons.bed

Voila:

	.. image:: https://raw.github.com/arq5x/scurgen/master/refseq.chr10.exons.bed.png
	

	

By default, we assume the chromosomes are humand (build 37 or hg19).  However, 
other genomes are allowed::

	python scurgen.py plot  --chrom chr10 \
	                        --color Greens \
	                        --format png \
	                        --genome mm9 \
	                        data/cpg-islands.hg19.chr10.bed
							

By default, we count each interval the same.  However, we can use weights such
as the score column (4th col.) in BEDGRAPH files::

	python scurgen.py plot  --chrom chr10 \
	                        --color Greens \
	                        --format png \
	                        --inc_column 4 \
	                        my.bedgraph

Full help::

	python scurgen.py plot --help
	usage: scurgeon plot [-h] [--genome GENOME] [--chrom CHROM]
	                     [--inc_col INC_COL] [--dim MATRIXDIM] [--color COLORMAP]
	                     [--format FORMAT] [--dpi DPI]
	                     file

	positional arguments:
	  file               The name of the file to be plotted.

	optional arguments:
	  -h, --help         show this help message and exit
	  --genome GENOME    The genome the dataset comes from (e.g., hg19)?
	  --chrom CHROM      The chrom that should be plotted (e.g., chr1)
	  --inc_col INC_COL  Use a specific column for incrementing file.
	  --dim MATRIXDIM    The dimensions of the curve. A power of 2.
	  --color COLORMAP   The name of the color map that should be used.
	  --format FORMAT    The type of output figure to create.
	  --dpi DPI          The resolution (in DPI) of the output.