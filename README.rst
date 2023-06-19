.. These are examples of badges you might want to add to your README:
   please update the URLs accordingly

    .. image:: https://api.cirrus-ci.com/github/<USER>/coffee.svg?branch=main
        :alt: Built Status
        :target: https://cirrus-ci.com/github/<USER>/coffee
    .. image:: https://readthedocs.org/projects/coffee/badge/?version=latest
        :alt: ReadTheDocs
        :target: https://coffee.readthedocs.io/en/stable/
    .. image:: https://img.shields.io/coveralls/github/<USER>/coffee/main.svg
        :alt: Coveralls
        :target: https://coveralls.io/r/<USER>/coffee
    .. image:: https://img.shields.io/pypi/v/coffee.svg
        :alt: PyPI-Server
        :target: https://pypi.org/project/coffee/
    .. image:: https://img.shields.io/conda/vn/conda-forge/coffee.svg
        :alt: Conda-Forge
        :target: https://anaconda.org/conda-forge/coffee
    .. image:: https://pepy.tech/badge/coffee/month
        :alt: Monthly Downloads
        :target: https://pepy.tech/project/coffee
    .. image:: https://img.shields.io/twitter/url/http/shields.io.svg?style=social&label=Twitter
        :alt: Twitter
        :target: https://twitter.com/coffee

.. image:: https://img.shields.io/badge/-PyScaffold-005CA0?logo=pyscaffold
    :alt: Project generated with PyScaffold
    :target: https://pyscaffold.org/
.. image:: https://readthedocs.org/projects/coffee-gr/badge/?version=latest
    :target: https://coffee-gr.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status


======
coffee
======


    Add a short description here!


Installation
------------

To use this software, you will need the following dependencies:

- fitsio C library
- chealpix C library
- mpich or Open MPI
- pyFFTW
- gnuplot + appropriate terminal if needed
- HDF library

Please make sure these dependencies are installed before using this software.

.. rubric:: Getting started

One needs to point the ``PYTHONPATH`` environment variable to the packages directory. Then the system code can be run with::

    python3.10 -O TwoDAdvection_setup.py

Note that adding a ``-d`` flag to this allows you to see a visualization during the computation. You should see a wave propagating from right to left with varying amplitude.


.. _pyscaffold-notes:

Note
====

This project has been set up using PyScaffold 4.4. For details and usage
information on PyScaffold see https://pyscaffold.org/.
