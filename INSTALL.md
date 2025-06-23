# Edgecopy Installation Instructions

To get started with Edgecopy, you'll need to install a few dependencies and the Edgecopy package itself.

## Prerequisites:

To get started with Edgecopy, you'll need to install a few dependencies and the Edgecopy package itself.

1. **Install Parascopy**:

Edgecopy depends on [Parascopy](https://github.com/tprodanov/parascopy). To install it, run the following commands:
```
conda config --add channels bioconda
conda config --add channels conda-forge
conda install -c bioconda parascopy
```

2. **Install ExomeDepth**:

Edgecopy also relies on [ExomeDepth](https://github.com/vplagnol/ExomeDepth). Install it with:
```
conda install -c bioconda r-exomedepth
```

## Additional Dependencies:
You'll also need these additional packages installed:
```
conda install pandas
conda install pyreadr
conda install networkx
conda install r-optparse
conda install bedtools
```

## Install Edgecopy:
After installing the prerequisites, you can install Edgecopy itself by cloning the GitHub repository and using `pip`:

```
git clone https://github.com/byunsy/edgecopy.git
cd edgecopy
pip install -e .
```