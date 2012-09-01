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

	python scurgen.py plot --chrom chr10 \
                            --color Greens \
                            --format png \
                            cpg-islands.hg19.chr10.bed

yields:

	.. image:: https://raw.github.com/arq5x/scurgeon/master/cpg-islands.hg19.chr10.bed.png
	
