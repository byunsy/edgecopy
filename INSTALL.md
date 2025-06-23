# Edgecopy Installation Instructions

Edgecopy is written in Python (≥ v3.6) and is available on Linux and macOS. 

## Prerequisites:

To get started with Edgecopy, you'll need to install a few dependencies and the Edgecopy package itself.

1. **Configure Conda Channels**:

The following commands help Conda to look in `bioconda` and `conda-forge` for packages when you install them.
```
conda config --add channels bioconda
conda config --add channels conda-forge
```

2. **Install Parascopy**:

Edgecopy utilizes several modules from [Parascopy](https://github.com/tprodanov/parascopy). To install it, run the following commands:
```
conda install -c bioconda parascopy
```
It typically takes 1-4 minutes to install Parascopy using `conda`.

3. **Install ExomeDepth**:

Edgecopy also relies on [ExomeDepth](https://github.com/vplagnol/ExomeDepth) for the `edgecopy depth` module. Install it with:
```
conda install -c bioconda r-exomedepth
```
In most cases, installing ExomeDepth can take around 2-5 minutes.

## Additional Dependencies:
You'll also need these additional packages installed:
```
conda install pandas
conda install pyreadr
conda install networkx
conda install r-optparse
conda install bedtools
```
Each package will take approximately one minute to install (five minutes in total for all five pacakges) using `conda`.

## Install Edgecopy:
Once the prerequisites are installed, you can set up Edgecopy (in under a minute) by cloning the GitHub repository and installing it with `pip`:

```
git clone https://github.com/byunsy/edgecopy.git
cd edgecopy
pip install -e .
```
