# Installation for SUSE

This was tested with AWS SUSE Linux Enterprise Server 15 SP1 (HVM).

## Install gcc & c++

```bash
sudo zypper in gcc-c++
```

## Install Miniconda

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

For more information:
- https://docs.conda.io/en/latest/miniconda.html
- https://conda.io/projects/conda/en/latest/user-guide/install/linux.html#install-linux-silent

Log out log back in.

## Create a Virtual Environment

```bash
conda create -n seqc python=3.8.2 pip
conda activate seqc
```

## Install dependencies

```
pip install Cython
pip install numpy
pip install bhtsne

conda install -c anaconda hdf5
conda install -c bioconda samtools
conda install -c bioconda star
```

## Install SEQC

```
wget https://github.com/dpeerlab/seqc/archive/v0.2.5.tar.gz
tar xvzf v0.2.5.tar.gz
cd seqc-0.2.5/
pip install .
```
