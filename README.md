# WMSA

A method of Multiple Sequence Alignment, which the writer is Wym6912.

![version](https://anaconda.org/wym6912/wmsa/badges/version.svg) ![last updated](https://anaconda.org/wym6912/wmsa/badges/latest_release_date.svg) ![supported platforms](https://anaconda.org/wym6912/wmsa/badges/platforms.svg) ![conda](https://anaconda.org/wym6912/wmsa/badges/installer/conda.svg)


## Recommended environment

Liunx-based systems, like `Ubuntu`, `CentOS` .

If you are a `Windows 10` user, you can use it by  [`WSL`]([Manually download Windows Subsystem for Linux (WSL) Distros | Microsoft Docs](https://docs.microsoft.com/en-us/windows/wsl/install-manual)) .

## How to use this program

```bash
git clone github.com/wym6912/WMSA --recursive 
# add --recursive arugment to get cd-hit and mafft subprograms
make all THREADS=16 # not -j16 because the THREADS=16 can be used by the subprograms
make install # program is installed on /usr/bin by default
wmsa -H # help message on this program
```

You also can download source code on the release.

If you are an MACOS user, please use the following `make` command:

```bash
make all THREADS=16 NOOPENMP=yes
```

in order to avoid the bug on `openmp` library.

You can also use our program in `conda` :

```bash
conda install wmsa -c wym6912 
```

## How to upgrade this program

```bash
git pull
git submodule foreach 'git pull'
```

...or you can download source code on the release.

## How to remove this program

```bash
make uninstall
```

