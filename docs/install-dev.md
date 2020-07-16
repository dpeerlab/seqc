# Setup for Development

Last verified: Jun 4, 2020

## Create Conda Environment

```bash
conda create -n seqc-debug python=3.8.2 pip
conda activate seqc-dev
```

## Install Dependencies

```bash
pip install Cython
pip install numpy
pip install bhtsne
```

For Mac (Mojave 10.14.6), install the following additional components. You must have `brew` to install.

```
brew install cairo
brew install pango
```

## Install SEQC (editable mode)

```bash
pip install --editable .
```

## Install STAR

curl -OL https://github.com/alexdobin/STAR/archive/2.5.3a.tar.gz
tar -xf 2.5.3a.tar.gz
cp STAR-2.5.3a/bin/MacOSX_x86_64/STAR /usr/local/bin/

## Install samtools

conda install -c bioconda samtools=1.3.1

## Install Packages for Testing

```bash
pip install nose
```

## Install Packages for Linting and Formating

```bash
pip install pylint
pip install autopep8
pip install black
```
