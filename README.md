# ARTEMIS
A framework with scripts for the analysis and clustering of molecular systems according to data obtained using the [PARENT](https://github.com/markusfleck/PARENT) or [PARENT_GPU](https://github.com/markusfleck/PARENT_GPU) package.

You will need g++ and python libraries to work: numpy, pandas, seaborn, json, matplotlib, sklearn and scipy.

![Framework scheme](framework_scheme.png) 

## 0) INSTALLATION
To download the console-app version, use:

> git clone https://github.com/nalsur-veallam/ARTEMIS.git

After that, enter the directory with the library

> cd ARTEMIS

And compile the C++ part of the library using make:

> make

Next, you can install all the dependencies for Python. Use:

> pip3 install -r requirements.txt

To install the library completely, run:

> pip3 install -e .

After this, the ARTEMIS program will be installed on your computer. Run

> artemis -h

to view available modules.

## 1) MANUAL AND TUTORIAL

Uploaded via [link](https://nalsur-veallam.github.io/TestPages/).
