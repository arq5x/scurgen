Overview
--------
`scurgen` is a command line tool for plotting genomic data using scpace-filling
curves such as the Hilbert Curve.

.. WARNING::
   This project is in its infancy.  Documenation is sparse, functionality is 
   evolving, and installation/dependencies are not documented.  

   
The plot module
---------------

This module allows you to make a basic plot of genomic data using a Hilbert 
curve.::

	python scurgen.py plot  --chrom chr10 \
                            --color Greens \
                            --format png \
                            data/cpg-islands.hg19.chr10.bed

yields a PNG image named cpg-islands.hg19.chr10.bed.png:

	.. image:: https://raw.github.com/arq5x/scurgen/master/cpg-islands.hg19.chr10.bed.png
	

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
