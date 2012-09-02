#!/usr/bin/env python

import sys
import  matplotlib.pyplot as plt
from matplotlib import cm as cm
import matplotlib.colors as colors
from pylab import *
import Image
import argparse
import hilbert as hb


def plot(parser, args):
    hm = hb.HilbertMatrix(args.file, args.genome, args.chrom, args.dim)
    hm.mask_zeros()
    cmap = cm.get_cmap(args.cmap, 5)
    
    fig = plt.figure()
    hilbert = plt.imshow(hm.matrix, interpolation='nearest', cmap=cmap)
    
    out_file = args.file + "." + args.format
    fig.savefig(out_file, dpi=args.dpi, transparent=True)
    plt.show()


def combine(parser, args):
    hm1 = hb.HilbertMatrix(args.file1, args.genome, args.chrom, args.dim, \
                           args.inc_col1 - 1)
    hm2 = hb.HilbertMatrix(args.file2, args.genome, args.chrom, args.dim, \
                           args.inc_col2 - 1)
                           
    hm1.mask_zeros()
    hm2.mask_zeros()
    
    # TO DO, set 0 == white (e.g., for centromere.)
    
    cmap1 = cm.get_cmap(args.cmap1, 20)
    cmap2 = cm.get_cmap(args.cmap2, 20)
    
    cmap1._init()
    cmap2._init() # create the _lut array, with rgba values
    
    # create your alpha array and fill the colormap with them.
    # here it is progressive, but you can create whathever you want
    alphas = np.linspace(0, 0.8, cmap2.N+3)
    cmap2._lut[:,-1] = alphas
    alphas = np.linspace(0, 0.8, cmap1.N+3)
    cmap1._lut[:,-1] = alphas
    
    cmap1.set_bad('w')
    cmap2.set_bad('w')
    
    # build fig1
    fig = plt.figure()

    hilb1 = plt.imshow(hm1.matrix, interpolation='nearest', \
                       cmap=cmap1, origin='lower')
    fig.savefig("_a.png", dpi=args.dpi, transparent=True)
    
    # build fig2
    fig = plt.figure()
    hilb2 = plt.imshow(hm2.matrix, interpolation='nearest', \
                       cmap=cmap2, origin='lower')
    fig.savefig("_b.png", dpi=args.dpi, transparent=True)
    
    # merge them
    background = Image.open("_a.png")
    foreground = Image.open("_b.png")
    
    #composite = Image.blend(background, foreground, 0.5)
    
    background.paste(foreground, (0, 0), foreground)
    background.show()
    #composite.show()
    plt.show()



def main():
    #########################################
    # create the top-level parser
    #########################################
    parser = argparse.ArgumentParser(prog='scurgeon')
    subparsers = parser.add_subparsers(title='[sub-commands]')
    
    #########################################
    # scurgeon 'plot'
    #########################################
    parser_plot = subparsers.add_parser('plot',  help='plot a single dataset')
    parser_plot.add_argument('file', metavar='file',  
                            help='The name of the file to be plotted.')

    parser_plot.add_argument('--genome', dest='genome', metavar='GENOME',
                        help='The genome the dataset comes from (e.g., hg19)?',
                        default='hg19')

    parser_plot.add_argument('--chrom', dest='chrom', metavar='CHROM', 
                        help='The chrom that should be plotted (e.g., chr1)')
    
    parser_plot.add_argument('--inc_col', dest='inc_col', \
                                metavar='INC_COL', 
                                help='Use a specific column for incrementing \
                                      file.',
                                default=-1)

    parser_plot.add_argument('--dim', dest='dim', metavar='MATRIXDIM', 
                        help='The dimensions of the curve.  A power of 2.',
                        type=int,
                        default=256)
    
    parser_plot.add_argument('--color', dest='cmap', metavar='COLORMAP', 
                        help='The name of the color map that should be used.')
                        
    parser_plot.add_argument('--format', dest='format', metavar='FORMAT', 
                        help='The type of output figure to create.',
                        default='png')
                        
    parser_plot.add_argument('--dpi', dest='dpi', metavar='DPI', 
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