

Installation
============
GQT depends on htslib and hdf5.

Install htslib
::

    git clone https://github.com/samtools/htslib.git
    cd htslib
    make


Install hdf5
::

We recommend downloading one of the statically-linked binary distributions from
here: http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8.7/obtain5187.html.

Modify GQT Makiefile
::

Next, modify the `HTS_ROOT` and `HDF_ROOT` variables in src/c/Makfile to
reflect ther locations.

Compile GQT
::

    cd gqt/
    make

Test GQT
::

    cd gqt/src/test/unit
    make
    cd ../func
    bash functional_tests.sh

