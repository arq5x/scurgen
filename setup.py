import ez_setup
ez_setup.use_setuptools()

import glob
import os
import sys
from setuptools import setup
from distutils.extension import Extension

if 'setuptools.extension' in sys.modules:
    m = sys.modules['setuptools.extension']
    m.Extension.__dict__ = m._Extension.__dict__

version_py = os.path.join(os.path.dirname(__file__), 'scurgen', 'version.py')
version = open(version_py).read().strip().split('=')[-1].replace('"','')
long_description = """
``scurgen`` is a tool for exploring genomic data with space filling curves'
"""

setup(
        name="scurgen",
        version=version,
        install_requires=['numpy', 'pybedtools>=0.6.1', 
                          'matplotlib>=1.0.0', 'PIL>=1.1.7', 
                          'bx-python>=0.7.1'],
        requires = ['python (>=2.5, <3.0)'],
        packages=['scurgen',
                  'scurgen.scripts',
                  'scurgen.data'],
        author="Aaron Quinlan",
        description='A tool for exploring genomic data with space \
                    filling curves',
        long_description=long_description,
        url="none",
        package_dir = {'scurgen': "scurgen"},
        zip_safe = False,
        scripts = ['scurgen/scripts/scurgen'],
        author_email="arq5x@virginia.edu",
        classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: GNU General Public License (GPL)',
            'Topic :: Scientific/Engineering :: Bio-Informatics']
    )
