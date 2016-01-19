==============
Notes on CTAUX
==============

CTAUX Installation instructions
=========================

#. Prerequisite: ``FFTW3`` ::
     export FFTW3_DIR=/where/FFTW3/installed
#. Prerequisite: ``NFFT`` ::
     export NFFT3_DIR=/where/NFFT3/installed
#. Prerequisite: ``ALPSCore`` ::
     export ALPSCore_DIR=/where/ALPSCore/installed
#. Prerequisite: ``DMFT`` ::
     export DMFT_DIR=/where/DMFT/source
     dmft_build=/where/DMFT/built
#. Download::
     ctaux_src=$PWD/CTAUX
     git clone git@github.com:egull/CTAUX.git $ctaux_src
#. Generate::
     ctaux_build=$ctaux_src/000build
     mkdir -p $ctaux_build
     cd $ctaux_build
     cmake -DLIBCLUSTER_DIR=$dmft_build/libcluster $ctaux_src 
     


DMFT installation instructions
==============================

#. Prerequisite: ``ALPSCore`` ::
     export ALPSCore_DIR=/where/ALPSCore/installed
#. Prerequisite: ``GFTools`` (will later become a part of ALPSCore)::
     gftools_src=$PWD/GFTools
     git@github.com:ALPSCore/GFTools.git $gftools_src
     gftools_build=$gftools_src/000build
     mkdir -p $gftools_build
     cd $gftools_build
     cmake -DCMAKE_INSTALL_PREFIX=$PWD $gftools_src
     make
     make test
     make install # ! required
     export GFTools_DIR=$gftools_build

#. Download::
     dmft_src=$PWD/DMFT
     git clone ssh://github.com/CQMP/DMFT.git $dmft_src
     git checkout GFTools
#. Generate::
     dmft_build=$dmft_src/000build
     mkdir -p $dmft_build
     cd $dmft_build
     cmake -DCMAKE_INSTALL_PREFIX=$PWD $dmft_src
#. Build and install ::
     make install
     export DMFT_DIR=$dmft_src



General considerations about CMake modules
==========================================

Ideally all we need is to provide an external target, transitively
containing all dependencies. For now, all we need is to be able to
locate the library and the header file. For an uninstalled package,
``${CMAKE_BINARY_DIR}`` is likely to contain the library and
``${CMAKE_SOURCE_DIR}`` is likely to contain the headers. Moreover, it
should be always possible to install the package to its build
directory.


