# CLUSTERING
A framework with scripts for the analysis and clustering of molecular systems according to data obtained using the [PARENT](https://github.com/markusfleck/PARENT) or [PARENT_GPU](https://github.com/markusfleck/PARENT_GPU) package.

You will need g++ and python libraries to work: numpy, pandas, seaborn, json, pylab, sklearn and scipy

## QUICK START
Set the necessary parameters in the ff.sh and hh.sh scripts and run in turn:
> bash gen_map.sh
> bash clustering.sh
> bash analysis.sh

## SHORT DESCRIPTION
The last report with a description is in **"report.pdf"** file. Here is a brief description of the framework blocks:

### CONVERTER

The converter from a binary .par file produces a matrix of mutual information between amino acid residues in json format. The converter is completely written in C ++ in the form of two versions (for GPU and CPU PARENT, but in the end the only difference is that the CPU version takes into account all degrees of freedom, and the GPU uses only dihedral angles, and they read the binary file in the same way). It is advised to use the GPU version. No libraries other than STL are required to run this code. To run, specify your C++ compiler in the Makefile (g++ by default). Next use:
> make clean
> make 
> bin/get_map_gpu -f input.par -n project_name

where after the -f flag is the path to the binary file, and after the -n flag is the desired project name

### DRAWING THE MATRIX

Script src/python/draw_map.py will draw the resulting matrix. To run use:
> python src/python/draw_map.py -n project_name

draws a matrix to the ./output/map/ directory

### CONSTRUCTION OF METRIC GRAPHS FOR EVALUATION OF THE OPTIMAL NUMBER OF CLUSTERS

Script src/python/opt_num_of_clust.py will draw this graphs. To run use:
> python src/python/opt_num_of_clust.py -n project_name -min min_num_of_clust -max max_num_of_clust

draws a graphs to the ./output/map/ directory. The -min and -max flags are optional and reflect the minimum and maximum number of clusters in the system

### CLUSTERING 

The src/python/clustering.py script performs clustering according to the existing matrix and builds a hierarchical clustering dendrogram in the ./output/clusterization/ directory. To run use:
> python src/python/clustering.py -n project_name -nclust num_of_clust

where after the -nclust flag is the desired number of clusters.

