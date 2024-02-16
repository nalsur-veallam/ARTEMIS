# ARTEMIS
A framework with scripts for the analysis and clustering of molecular systems according to data obtained using the [PARENT](https://github.com/markusfleck/PARENT) or [PARENT_GPU](https://github.com/markusfleck/PARENT_GPU) package.

You will need g++ and python libraries to work: numpy, pandas, seaborn, json, matplotlib, sklearn and scipy.

![Framework scheme](framework_scheme.png) 

## 0) INSTALLATION
To download the console-app version, use:

> git clone https://github.com/nalsur-veallam/ARTEMIS.git -b console-app

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

## 1) EXAMPLE

To test use:

> cd example

> mkdir output

> artemis map v536e_1ps.par v536e_2ps.par -dt1 1 -dt2 2 --denoise -o output/map.json

> artemis map output/map.json --draw -norm -o output/map.pdf

> artemis allostery output/map.json v536e_groups.json --search -noseq 2 -o output/v536_intensity.pdf

> artemis allostery output/v536_intensity.json --draw -strc v536e.pdb -o output/v536_intensity.pse -noseq 2

> artemis allostery output/v536_intensity.json --analysis -o output/v536_intensity.pdf -noseq 2 -zscore
