Overview
--------
`scurgen` (/skur'jən/) is a command line tool for plotting genomic data 
using scpace-filling curves such as the Hilbert Curve.

.. WARNING::
   This project is in its infancy.  Documentation is sparse, functionality is 
   evolving, and installation/dependencies are not documented.  

   
The plot module
---------------

This module allows you to make a basic plot of genomic data using a Hilbert 
curve.::

	python scurgen.py plot  --chrom chr10 \
	                        --cmap Greens \
	                        --format png \
	                        --dim 128 \
	                        data/cpg-islands.hg19.chr10.bed

yields a PNG image named cpg-islands.hg19.chr10.bed.png:

	.. image:: https://raw.github.com/arq5x/scurgen/master/cpg-islands.hg19.chr10.bed.png

And let's do the same with 5' UTRs.::

	python scurgen.py plot  --chrom chr10 \
	                        --cmap Blues \
	                        --format png \
	                        --dim 128 \
	                        data/refseq.chr10.exons.bed

Voila:

	.. image:: https://raw.github.com/arq5x/scurgen/master/refseq.chr10.exons.bed.png
	

Now, when one combines the two, one can see the co-localization (or lack of)
of the two sets in the genome.  The following image was merged manually with
image software, but I am working on dependable methods for doing this
programmatically:

	.. image:: https://raw.github.com/arq5x/scurgen/master/merge.png
	

By default, we assume the chromosomes are humand (build 37 or hg19).  However, 
other genomes are allowed::

	python scurgen.py plot  --chrom chr10 \
	                        --cmap Greens \
	                        --format png \
	                        --genome mm9 \
	                        data/cpg-islands.hg19.chr10.bed
							

Support for quantitative data
-----------------------------
In the case of quantitative experiments such a ChIP-seq, we often have
quantitative values associated with genomic intervals (e.g., "peaks").
The `inc-col` option allows one to use a specific column in a text file
such as BEDGRAPH to increment the scores associated with each cell in 
the curve matrix.  Here is an example of using the pValue column from
an ENCODE narrow peak (http://genome.ucsc.edu/FAQ/FAQformat#format12) file::

	python scurgen.py plot  --chrom chr10 \
	                        --cmap Reds \
	                        --format png \
	                        --dim 512 \
	                        --inc_col 8 \
	                        data/wgEncodeSydhTfbsGm12878Nrf1IggmusPk.narrowPeak
	
yielding:

	.. image:: https://raw.github.com/arq5x/scurgen/master/wgEncodeSydhTfbsGm12878Nrf1IggmusPk.narrowPeak.png 

Full help for plot module
-------------------------

	python scurgen.py plot --help
	usage: scurgen plot [-h] [--genome STRING] [--chrom STRING]
	                    [--inc_col INC_COL] [--dim INTEGER] [--min_mask INTEGER]
	                    [--cmap STRING] [--format STRING] [--dpi INTEGER]
	                    file

	positional arguments:
	  file                The name of the file to be plotted.

	optional arguments:
	  -h, --help          show this help message and exit
	  --genome STRING     The genome the dataset comes from (e.g., hg19)?
	  --chrom STRING      The chrom that should be plotted (e.g., chr1)
	  --inc_col INC_COL   Use a specific column for incrementing file.
	  --dim INTEGER       The dimensions of the curve. A power of 2.
	  --min_mask INTEGER  The minimum value allowed in a cell before it is set to
	                      the masking color (white).
	  --cmap STRING       The name of the matplotlib color map that should be
	                      used. See scipy.org/Cookbook/Matplotlib/Show_colormaps
	                      for options.
	  --format STRING     The type of output figure to create.
	  --dpi INTEGER       The resolution (in DPI) of the output.