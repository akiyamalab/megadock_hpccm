# MEGADOCK-HPCCM

## Description

MEGADOCK-HPCCM is a HPC container making workflow for [MEGADOCK](https://github.com/akiyamalab/MEGADOCK) application on HPC environments by using [HPCCM (HPC Container Maker)](https://github.com/NVIDIA/hpc-container-maker/) framework.
It generates container specification (recipe) files both in [Dockerfile](https://docs.docker.com/engine/reference/builder/) and [Singularity definition](https://sylabs.io/guides/3.3/user-guide/definition_files.html) format from one simple python recipe.
Dependent libraries (GPU, OpenMPI, FFTW, InfiniBand, Intel Omni-Path) will be properly configured by setting parameters to use MEGADOCK application on HPC environments.
It gives users an easier way to use MEGADOCK application containers when considering the specification differences between the hosts and containers in multiple HPC environments.

## Requirements

- [NVIDIA GPU devices, drivers](https://www.nvidia.com/)
- [HPC Container Maker](https://github.com/NVIDIA/hpc-container-maker/)
- [Docker](https://www.docker.com/) (if you use)
- [Singularity](https://sylabs.io/) (if you use)

## Repository Overview
```
.
├── data                                #  
│   └── ...                             #  sample of input file (.pdb, .table)
├── sample                              #  
│   ├── Dockerfile                      #  for general Docker environments
│   ├── singularity_ompi-2-1-2_opa.def  #  for TSUBAME3.0 (ompi=2.1.2, opa=True)
│   └── singularity_ompi-2-1-6_ofed.def #  for ABCI (ompi=2.1.6, ofed=True)
├── script                              #  
|   └── makeTable.sh                    #  script for generating input table (.table)
├── megadock-scfa20                     #  source code of MEGADOCK 5.0 (alpha)
├── megadock_hpccm.py                   #  MEGADOCK-HPCCM for HPCCM framework
├── Makefile                            #  
└── README.md                           #  this document

# The following directories will be generated after running scripts
.
├── table                           # directory for storing metadata files
└── out                             # directory for storing output files
```

----

## Quick Links

- [MEGADOCK-HPCCM](#megadock-hpccm)
  - [Description](#description)
  - [Requirements](#requirements)
  - [Repository Overview](#repository-overview)
  - [Quick Links](#quick-links)
  - [Docker Environment](#docker-environment)
    - [Requirements](#requirements-1)
    - [1. Setting up (HPCCM)](#1-setting-up-hpccm)
    - [2. Generate Dockerfile](#2-generate-dockerfile)
    - [3. Build Docker image](#3-build-docker-image)
    - [4. Test with sample protein-protein pairs](#4-test-with-sample-protein-protein-pairs)
    - [5. Test with ZLAB benchmark dataset](#5-test-with-zlab-benchmark-dataset)
  - [Singularity Environment](#singularity-environment)
    - [Requirements](#requirements-2)
    - [1. Setting up (HPCCM)](#1-setting-up-hpccm-1)
    - [2. Generate Singularity Definition](#2-generate-singularity-definition)
    - [3. Build Singularity image](#3-build-singularity-image)
    - [4. Test with sample protein-protein pairs](#4-test-with-sample-protein-protein-pairs-1)
    - [5. Test with ZLAB benchmark dataset](#5-test-with-zlab-benchmark-dataset-1)

----

## Docker Environment

### Requirements

- pip, python (for HPCCM)
- docker ( > 19.03 )
  - or `nvidia-docker` for gpu support

### 1. Setting up (HPCCM)

```sh
# install hpccm
sudo pip install hpccm

# clone MEGADOCK-HPCCM repository
git clone https://github.com/akiyamalab/megadock_hpccm.git
cd megadock_hpccm
```

### 2. Generate Dockerfile

``` sh
# generate 'Dockerfile' from hpccm recipe
hpccm --recipe megadock_hpccm.py --format docker > Dockerfile

## or adding 'userarg' for specifying library versions
hpccm --recipe megadock_hpccm.py --format docker --userarg ompi=2.1.2 fftw=3.3.8 > Dockerfile

## Available userargs:
##  ompi=${ompi_version} : version of OpenMPI library
##  fftw=${fftw_version} : version of FFTW library
```

### 3. Build Docker image

```sh
# build a container image from Dockerfile
docker build . -f Dockerfile -t megadock:hpccm
```

### 4. Test with sample protein-protein pairs

```sh
# run with host gpus
docker run --rm -it --gpus all \
  -v `pwd`/data:/data  megadock:hpccm \
  mpirun --allow-run-as-root -n 2 /workspace/megadock-gpu-dp -tb /data/SAMPLE.table
```

### 5. Test with ZLAB benchmark dataset

```sh
# clone MEGADOCK-HPCCM repository
git clone https://github.com/akiyamalab/megadock_hpccm.git
cd megadock_hpccm

# download benchmark dataset (ZDOCK Benchmark 5.0)
mkdir -p data
wget https://zlab.umassmed.edu/benchmark/benchmark5.tgz
tar xvzf benchmark5.tgz -C data
rm -f benchmark5.tgz

# create docking table using script (only 100 pairs)
INTERACTIVE=1 TABLE_ITEM_MAX=100 RUNTIME_RELATIVE_ROOT=/ script/makeTable.sh . data/benchmark5/structures/ \*_r_b.pdb \*_l_b.pdb test100pairs

# Note: 
# - unset ${TABLE_ITEM_MAX} variable to unlimit the number of docking calculations (all-to-all)
# - if you need to change the repository root path when runtime, use ${RUNTIME_RELATIVE_ROOT} to modify path in generating the table.

# run
docker run --rm -it --gpus all \
  -v `pwd`/data:/data -v `pwd`/table:/table -v `pwd`/out:/out \
  megadock:hpccm \
    mpirun --allow-run-as-root -n 2 -x OMP_NUM_THREADS=20 \
      /workspace/megadock-gpu-dp -tb /table/test100pairs/test100pairs.table
```

----

## Singularity Environment

### Requirements

- pip, python (for HPCCM)
- singularity
  - require `singularity exec` command on HPC system
  - require privilege for `sudo singularity build` or `singularity build --fakeroot` (>= 3.3)

Note: Following commands should be executed on your local environment where you have system privilege.

### 1. Setting up (HPCCM)

```sh
# install hpccm
sudo pip install hpccm

# clone MEGADOCK-HPCCM repository
git clone https://github.com/akiyamalab/megadock_hpccm.git
cd megadock_hpccm
```

### 2. Generate Singularity Definition

``` sh
# generate 'singularity.def' from hpccm recipe
hpccm --recipe megadock_hpccm.py --format singularity > singularity.def

## or adding 'userarg' for specifying library versions
hpccm --recipe megadock_hpccm.py --format singularity --userarg ompi=2.1.6 fftw=3.3.8 ofed=True > singularity.def

## Available userargs:
##  ompi=${ompi_version} : version of OpenMPI library
##  fftw=${fftw_version} : version of FFTW library
##  ofed=${True|False} : flag for install 'Mellanox OpenFabrics Enterprise Distribution for Linux'
##  opa=${True|False} : flag for install Intel Ompni-Path dependencies
```

### 3. Build Singularity image

```sh
# build a container image from Dockerfile
sudo singularity build megadock-hpccm.sif singularity.def

## or '.simg' format (singularity < 3.2)
sudo singularity build megadock-hpccm.simg singularity.def
```

### 4. Test with sample protein-protein pairs

- **Notes:**
  - Following commands should be running on HPC environment (compute-node with gpus).
  - Please replace `${SINGULARITY_IMAGE}` to **path to the container image file** on your environment.
  - **Please read the 'Singularity' section of system manual** which provided by your HPC system. We must add specific options for singularity runtime when using system resources.
    - e.g.) Volume option (`-B XXX`) for mounting system storage, applications, libraries, etc.

```sh
# clone MEGADOCK-HPCCM repository
git clone https://github.com/akiyamalab/megadock_hpccm.git
cd megadock_hpccm

# singularity exec 
singularity exec --nv ${SINGULARITY_IMAGE} \
  mpirun -n 2 /workspace/megadock-gpu-dp -tb data/SAMPLE.table
```

### 5. Test with ZLAB benchmark dataset

```sh
# clone MEGADOCK-HPCCM repository
git clone https://github.com/akiyamalab/megadock_hpccm.git
cd megadock_hpccm

# download benchmark dataset (ZDOCK Benchmark 5.0)
mkdir -p data
wget https://zlab.umassmed.edu/benchmark/benchmark5.tgz
tar xvzf benchmark5.tgz -C data
rm -f benchmark5.tgz

# create docking table using script (only 100 pairs)
INTERACTIVE=1 TABLE_ITEM_MAX=100 script/makeTable.sh . data/benchmark5/structures/ \*_r_b.pdb \*_l_b.pdb test100pairs

# Note: 
# - unset ${TABLE_ITEM_MAX} variable to unlimit the number of docking calculations (all-to-all)
# - if you need to change file path in compute-node, use ${RUNTIME_RELATIVE_ROOT} to modify path in generating the table.

# ${SINGULARITY_IMAGE}: path to the singularity image file

# singularity exec 
singularity exec --nv ${SINGULARITY_IMAGE} \
  mpirun -n 2 -x OMP_NUM_THREADS=20 \
  /workspace/megadock-gpu-dp -tb table/test100pairs/test100pairs.table

# singularity exec (with host MPI library)
mpirun -n 2 -x OMP_NUM_THREADS=20 \
  singularity exec --nv ${SINGULARITY_IMAGE} \
  /workspace/megadock-gpu-dp -tb table/test100pairs/test100pairs.table
```
