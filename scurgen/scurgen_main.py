#!/usr/bin/env python

import sys
import os
import  matplotlib.pyplot as plt
from matplotlib import cm as cm
import matplotlib.colors as colors
from pylab import *
import Image
import argparse
import hilbert as hb


def plot(parser, args):
    hm = hb.HilbertMatrix(args.file, args.genome, args.chrom, args.dim,
                          args.inc_col)
    hm.mask_low_values(args.min_mask)
    cmap = cm.get_cmap(args.cmap, 12)
    
    fig = plt.figure()
    hilbert = plt.imshow(hm.matrix, interpolation='nearest', cmap=cmap)
    plt.colorbar()
    
    out_file = args.file + "." + args.format
    fig.savefig(out_file, dpi=args.dpi, transparent=True)
    
    plt.show()


def combine(parser, args):
    hm1 = hb.HilbertMatrix(args.file1, args.genome, args.chrom, args.dim, \
                           args.inc_col1)
    hm2 = hb.HilbertMatrix(args.file2, args.genome, args.chrom, args.dim, \
                           args.inc_col2)
                           
    hm1.mask_low_values()
    hm2.mask_low_values()
    
    # TO DO, set 0 == white (e.g., for centromere.)
    
    cmap1 = cm.get_cmap(args.cmap1, 5)
    cmap2 = cm.get_cmap(args.cmap2, 5)
    

    
    cmap1.set_bad('w')
    cmap2.set_bad('w')
    
    # build fig1
    fig = plt.figure()
    hilb1 = plt.imshow(hm1.matrix, interpolation='nearest', cmap=cmap1)
    fig.savefig("_a.png", dpi=args.dpi, transparent=True)
    
    # build fig2
    fig = plt.figure()

    hilb2 = plt.imshow(hm2.matrix, interpolation='nearest', cmap=cmap2)
    fig.savefig("_b.png", dpi=args.dpi, transparent=True)
    
    # merge them
    background = Image.open("_a.png")
    foreground = Image.open("_b.png")
    
    foreground.load()
    background.load()
    
    (r, g, b, a) = foreground.split()
    foreground = Image.merge("RGB", (r, g, b))
    mask = Image.merge("L", (a,))
    background.paste(foreground, (0, 0), mask)
    #background.paste(foreground, (0, 0), foreground)
    background.save("merge.png")
    
    # cleanup temp files
    os.remove("_a.png")
    os.remove("_b.png")


def main():
    #########################################
    # create the top-level parser
    #########################################
    parser = argparse.ArgumentParser(prog='scurgen')
    subparsers = parser.add_subparsers(title='[sub-commands]')
    
    #########################################
    # scurgeon 'plot'
    #########################################
    parser_plot = subparsers.add_parser('plot',  help='plot a single dataset')
    parser_plot.add_argument('file', metavar='file',  
                            help='The name of the file to be plotted.')

    parser_plot.add_argument('--genome', dest='genome', metavar='STRING',
                        help='The genome the dataset comes from (e.g., hg19)?',
                        default='hg19')

    parser_plot.add_argument('--chrom', dest='chrom', metavar='STRING', 
                        help='The chrom that should be plotted (e.g., chr1)')
    
    parser_plot.add_argument('--inc_col', dest='inc_col', \
                                metavar='INC_COL', 
                                help='Use a specific column for incrementing \
                                      file.',
                                type=int,
                                default=None)

    parser_plot.add_argument('--dim', dest='dim', metavar='INTEGER', 
                        help='The dimensions of the curve.  A power of 2.',
                        type=int,
                        default=256)
    
    parser_plot.add_argument('--min_mask', dest='min_mask', metavar='INTEGER', 
                        help='The minimum value allowed in a cell before \
                        it is set to the masking color (white).',
                        type=int,
                        default=0)
    
    parser_plot.add_argument('--cmap', dest='cmap', metavar='STRING', 
                        help='The name of the matplotlib color map that \
                              should be used. See \
                              scipy.org/Cookbook/Matplotlib/Show_colormaps \
                              for options.')
                        
    parser_plot.add_argument('--format', dest='format', metavar='STRING', 
                        help='The type of output figure to create.',
                        default='png')
                        
    parser_plot.add_argument('--dpi', dest='dpi', metavar='INTEGER', 
                        help='The resolution (in DPI) of the output.',
                        type=int,
                        default=150)
    
    parser_plot.set_defaults(func=plot)
    
    #########################################
    # scurgeon 'combine'
    #########################################
    parser_combine = \
        subparsers.add_parser('combine',  help='combine two datasets')

    parser_combine.add_argument('-1', dest='file1', metavar='FILE1',
                        help='The first file to be plotted.')

    parser_combine.add_argument('-2', dest='file2', metavar='FILE2',
                        help='The second file to be plotted.')

    parser_combine.add_argument('--cmap1', dest='cmap1', metavar='COLORMAP', 
                        help='The name of the color map for file 1.')

    parser_combine.add_argument('--cmap2', dest='cmap2', metavar='COLORMAP', 
                        help='The name of the color map for file 2.')
                        
    parser_combine.add_argument('--inc_col1', dest='inc_col1', \
                                metavar='INC_COL1', 
                                help='Use a specific column for incrementing \
                                      file 1.',
                                type=int,
                                default=-1)

    parser_combine.add_argument('--inc_col2', dest='inc_col2', \
                                metavar='INC_COL2', 
                                help='Use a specific column for incrementing \
                                      file 2.',
                                type=int,
                                default=-1)
                        
    parser_combine.add_argument('--genome', dest='genome', metavar='GENOME',
                        help='The genome the dataset comes from (e.g., hg19)?',
                        default='hg19')
    parser_combine.add_argument('--chrom', dest='chrom', metavar='CHROM', 
                        help='The chrom that should be plotted (e.g., chr1)')
    
    parser_combine.add_argument('--dim', dest='dim', metavar='MATRIXDIM', 
                        help='The dimensions of the curve.  A power of 2.',
                        type=int,
                        default=256)
                          
    parser_combine.add_argument('--format', dest='format', metavar='FORMAT', 
                        help='The type of output figure to create.',
                        default='png')
                        
    parser_combine.add_argument('--dpi', dest='dpi', metavar='DPI', 
                        help='The resolution (in DPI) of the output.',
                        default=150)

    parser_combine.set_defaults(func=combine)
    
    
    #########################################
    # scurgeon
    #########################################    
    args = parser.parse_args()
    args.func(parser, args)
    
    
if __name__ == "__main__":
    main()