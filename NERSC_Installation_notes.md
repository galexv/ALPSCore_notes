The purpose of this document is to provide you with a step-by-step tutorial for installing ALPSCore and additional libraries you might need for your solver, in this case on [NERSC](https://www.nersc.gov) Cray machines, for example, on [Cori](https://www.nersc.gov/systems/cori/).

Here we assume that:

1. You have an access to the machine.
2. You use `bash` shell.
3. You have access to the GitHub repositories with the codes.
4. You are fine to use the default Cray development environment: compiler, linear algebra library, and MPI library.

Let's begin!

# 1. Prepare the directory structure and load the necessary software modules

1. First of all, load the necessary modules for CMake, Boost, HDF5, and FFTW3.
```bash
$ module add cmake boost cray-hdf5 fftw
```
and declare C and C++ compilers to be used:
```bash
$ export CC=$(which cc)
$ export CXX=$(which CC)
```

2. It is convenient to decide where you will download, build and install the software. For example, it can be somewhere in your home directory, on in a "scratch" directory (pointed by environment variable `SCRATCH` on Cori), or in your project directory. 

You may be able to find the names of your project directories by the following command:
```bash
$ (for f in $(id -Gn); do d=/global/project/projectdirs/$f; [ -r $d ] && [ -w $d ] && [ -x $d ] && ls -ld $d; done)
```

3. So, let's assume your project directory name is `/global/project/projectdirs/m1234567` and it is in its subdirectory we want to install the code. Let's put this directory name in a shell variable `my_base_dir`, for convenience, and make sure that the directory is created.

```bash
$ my_base_dir=/global/project/projectdirs/m1234567/$USER/alpscore_stuff
$ mkdir -pv $my_base_dir
```

4. Assign the directories for source, build and installation.
    * **ALPSCore directories:** Note that `ALPSCore_DIR` variable *must* be exported and *must* point to the directory where ALPSCore is to be installed.
    ```bash
    $ alpscore_src=$my_base_dir/ALPSCore
    $ alpscore_build=$alpscore_src/000build
    $ export ALPSCore_DIR=$my_base_dir/install
    ```
    * **NFFT3 directories:**
    ```bash
    $ nfft3_src=$my_base_dir/nfft3
    $ export NFFT3_DIR=$my_base_dir/install/nfft3
    ```

# 2. Download and install ALPSCore
(See [ALPSCore installation Wiki page](https://github.com/ALPSCore/ALPSCore/wiki/Installation) for details.)

Download:
```bash
$ git clone https://github.com/ALPSCore/ALPSCore.git $alpscore_src
```

Build and install. Note that on Cray we have to compile ALPSCore statically:
```bash
$ mkdir -pv $alpscore_build
$ cd $alpscore_build
$ cmake $alpscore_src -DALPS_BUILD_TYPE=static -DCMAKE_INSTALL_PREFIX=$ALPSCore_DIR
$ make -j4
$ make install
```

Verify that installation is successful:
```bash
$ (cd $ALPSCore_DIR/lib && ls -1 libalps*)
```
should produce the output like this:
```
libalps-accumulators.a
libalps-gf.a
libalps-hdf5.a
libalps-mc.a
libalps-params.a
libalps-utilities.a
```

# 3. Download and install NFFT3
Download [NFFT3 library](https://www-user.tu-chemnitz.de/~potts/nfft/) from https://www-user.tu-chemnitz.de/~potts/nfft/ and unpack it.
```bash
$ mkdir -pv $nfft3_src
$ cd $nfft3_src
$ wget https://www-user.tu-chemnitz.de/~potts/nfft/download/nfft-3.3.2.tar.gz
$ tar -xzf nfft-3.3.2.tar.gz
```

Prepare the build. Note that we have to inform NFFT3 of the location of the FFTW3 library, and also we disable the generation of the shared libraries.
```bash
$ cd nfft-3.3.2/
$ ./configure --prefix=$NFFT3_DIR CC=$CC --with-fftw3=$FFTW_DIR/.. --enable-shared=no --enable-static=yes
```

Unfortunately, there seems to be a bug either in GNU Autotools or in NFFT3: the `configure` script generates spurious link arguments that will cause the build to fail. We have to fix it and regenerate the build:
```bash
$ sed -e 's/-lne required -lm//' <config.status >config.status.fixed
$ make distclean
$ cp config.status.fixed config.status
$ chmod +x config.status
$ ./config.status
```

Then build and install:

```bash
$ make -j4
$ make install
```

Verify that installation is successful:
```bash
$ (cd $NFFT3_DIR && ls -1dF lib/* include/*)
```
should produce the output like this:
```
include/nfft3.h
include/nfft3mp.h
lib/libnfft3.a
lib/libnfft3.la*
lib/pkgconfig/
```


