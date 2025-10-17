# ARTEMIS
An information theory based numerical framework for analysis of intramolecular communication using molecular dynamics (MD) data processed with the [PARENT](https://github.com/markusfleck/PARENT) or [PARENT_GPU](https://github.com/markusfleck/PARENT_GPU) software packages. This tool gives the possibility to analyze intra- and intermolecular communication networks using all-to-all mutual information (MI) matrix for amino acids in proteins, nucleotides in nucleic acids, or any other molecular subgroups defined in the system. Alternatively, for proteins and nucleic acids MI matrices can be generated from multiple sequence alignments (MSA) and further processed in the same way using ARTEMIS. 

A potentially wide application of ARTEMIS, as well as a description of the methodology, has been published as a pre-print on ChemRxiv.

![Framework scheme](framework_scheme.png) 

## 0) INSTALLATION
To download the console-app version, use:

```bash
    git clone https://github.com/nalsur-veallam/ARTEMIS.git
```

After that, enter the directory with the library

```bash
    cd ARTEMIS
```

And run:

```bash
    conda env create -f environment.yml
    conda activate artemis
```

After this, the ARTEMIS program will be installed in your conda environment. Run

```bash
    artemis -h
```

to view available modules.

## 1) MANUAL AND TUTORIAL

Uploaded via [link](https://nalsur-veallam.github.io/ARTEMIS/).
