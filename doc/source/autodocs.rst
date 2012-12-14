
API documentation
=================
Classes
=======
Base classes
------------
.. autosummary::
    :toctree: autodocs
    :template: auto_template

    scurgen.hilbert.HilbertBase
    scurgen.hilbert.HilbertNormalized
    scurgen.hilbert.HilbertMatrix
    scurgen.hilbert.HilbertMatrixBigWig

Plotting classes
----------------
.. autosummary::
    :toctree: autodocs
    :template: auto_template

    scurgen.plotting.HilbertPlot
    scurgen.plotting.HilbertGUI




Functions
=========

Manipulating Hilbert curve coordinates
--------------------------------------

.. autosummary::
    :toctree: autodocs

    scurgen.hilbert.d2xy
    scurgen.hilbert.d2rc
    scurgen.hilbert.rc2d
    scurgen.hilbert.xy2d
    scurgen.hilbert.get_interval_from_string

Helper functions
----------------

.. autosummary::
    :toctree: autodocs

    scurgen.plotting.data_dir


Debugging functions
-------------------
.. autosummary::
    :toctree: autodocs

    scurgen.plotting.debug_plot
    scurgen.plotting.plot_hilbert
    scurgen.plotting.gui_main
    scurgen.plotting._debug_HilbertPlot
    scurgen.plotting._debug_HilbertGUI
