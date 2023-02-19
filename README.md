# UPD: Recently updated the repository, but did not change the README (CAREFULLY! YOU NEED TO PROVIDE RAM AT LEAST THE SIZE OF THE .PAR FILE!)


# CLUSTERING
A framework with scripts for the analysis and clustering of molecular systems according to data obtained using the [PARENT](https://github.com/markusfleck/PARENT) or [PARENT_GPU](https://github.com/markusfleck/PARENT_GPU) package.

You will need g++ and python libraries to work: numpy, pandas, seaborn, json, pylab, sklearn and scipy.

## QUICK START
Set the necessary parameters in the gen_map.sh, tools.sh, clustering.sh and analysis.sh scripts and run in turn:

> bash gen_map.sh

> bash tools.sh

> bash clustering.sh

> bash analysis.sh

## SHORT DESCRIPTION

![Framework scheme](framework_scheme.png) 

The last report with a description is in **"report.pdf"** file. Here is a brief description of the framework blocks:

## BUILDING A MATRIX OF MUTUAL INFORMATION (present in gen_map.sh):

### CONVERTER

The converter from a binary .par file produces a matrix of mutual information between amino acid residues in json format. The converter is completely written in C ++. No libraries other than STL are required to run this code. To run, specify your C++ compiler in the Makefile (g++ by default). Next use (present in gen_map.sh):

> make clean

> make 

> bin/get_map -f input.par -n project_name

where after the **-f** flag is the path to the binary file, and after the **-n** flag is the desired project name.

### DRAWING THE MATRIX

* Script src/python/draw_map.py will draw the resulting matrix. To run use:

> python3 src/python/draw_map.py -n project_name

draws a matrix to the ./output/${project_name}/map/ directory. You can also use the **-nodiag** flag to zero out all diagonal matrix elements for better matrix contrast.

* Script src/python/filtration.py can filter the matrix of mutual information by the exposure of amino acid residues. This script uses pymol to run: the user specifies cutoff as a fraction of the maximum exposure that an amino acid residue in the structure must take in order to be filtered. Pymol calculates the solvent contact surface for each residue and compares it to the maximum value for each residue. To run use:

> python3 src/python/filtration.py -n project_name -strc sctructure.pdb(.gro ...) -cutoff cutoff

draws a matrix to the ./output/${project_name}/map/ directory. The structure must match the molecular dynamics trajectory that was used to derive the PARENT value. Cutoff is specified as a decimal fraction of the total exposure required for a residue to be considered "on the surface" of the protein (example: -cutoff 0.3, to cut out residues with an exposure ratio greater than 30%)

* Script src/python/filtration.py can filter the matrix of mutual information by the exposure of amino acid residues. This script uses GROMACS gmx sasa data (You need the data returned with the -or flag - the average solvent-free surface per residue.): the user specifies cutoff as a fraction of the maximum exposure that an amino acid residue in the structure must take in order to be filtered. To run use:

> python3 src/python/sasa_filtration.py -n project_name -strc sctructure.pdb(.gro ...) -cutoff cutoff -sasa sasa_file.xvg

draws a matrix to the ./output/${project_name}/map/ directory. The structure must match the molecular dynamics trajectory that was used to derive the PARENT value. Cutoff is specified as a decimal fraction of the total exposure required for a residue to be considered "on the surface" of the protein (example: -cutoff 0.3, to cut out residues with an exposure ratio greater than 30%). The data from the GROMACS is input in the .xvg format file.

## CLUSTERING AND CLUSTER ANALYSIS (present in clustering.sh) :

### CONSTRUCTION OF METRIC GRAPHS FOR EVALUATION OF THE OPTIMAL NUMBER OF CLUSTERS

Script src/python/opt_num_of_clust.py will draw this graphs. To run use:

> python3 src/python/opt_num_of_clust.py -n project_name -min min_num_of_clust -max max_num_of_clust

draws a graphs to the ./output/${project_name}/clustering/ directory. The **-min** and **-max** flags are optional and reflect the minimum and maximum number of clusters in the system.

### CLUSTERING 

The src/python/clustering.py script performs clustering according to the existing matrix and builds a hierarchical clustering dendrogram in the ./output/${project_name}/clustering/ directory. To run use:

> python3 src/python/clustering.py -n project_name -nclust num_of_clust

where after the **-nclust** flag is the desired number of clusters.

### CREATING A PYMOL SESSION

The src/python/create_pse.py script on your pdb file, passed after the **-f** flag, creates a pymol session with the selection of clusters into groups and saves it to the ./output/${project_name}/clustering/ directory. To run use:

> python3 src/python/clustering.py -f input.pdb -n project_name -nclust num_of_clust

### CLUSTER ANALYSIS

The src/python/cluster_analysis.py script analyzes the clusters built by the clustering.py program and saves it to the ./output/${project_name}/clustering/ directory. The result of the program is a pdf file, which will contain: 
* a matrix of mutual information between clusters, 
* information about the number of residues in each of the clusters,
* information entropy of each cluster, 
* mutual information of each cluster with an active site (+ normalized to the number of clusters), 
* mutual information of each cluster with an allosteric site (+ normalized to the number of clusters), 
* mutual information of each cluster with the rest of the protein (+ normalized to the number of clusters), 
* the sum of all mutual information of each cluster within itself (+ normalized to the number of clusters), 
* the number of active site residues in each cluster, 
* the number of allosteric site residues in each cluster.
To run use:

> python3 src/python/cluster_analysis.py -f_act act_site.json -asn act_site_name -f_all all_site.json -allsn all_site_name -n project_name -nclust num_of_clust -noseq n(default 0)

The **-f_act** flag sends a file in the .json format with indication of the active site residues in the form of numbers corresponding to their mention in the structure file when calculating the corresponding MD trajectory. The **-f_all** flag is used to pass a .json file indicating the remnants of the allosteric site in the same order. To understand this, you can look at the data with the key real_numbers and the key names in the .json file with the matrix of mutual information. To build such lists by ligand from a structure, you can use the find_area.py script (see the TOOLS section (provided in the tools.py script). If you do not know the active and/or allosteric site, then pass files with empty lists using these flags (it will made more convenient in future versions).The **-asn** and **-allsn** flags pass keys corresponding to the lists of active and allosteric sites in the corresponding .json files.The **-nclust** flag passes the number of clusters that have already been clustered (this script uses data obtained from the clustering script .py). The **-noseq** flag is optional. After it, you can specify the number of residues that should not be taken into account during allosteric interaction along the protein sequence. That is, if this parameter is equal to one, then when calculating the intensity of allosteric communication for the i-th residue, it will not be taken into account its mutual information with i+1 and i-1 residues.

## ALLOSTERY ANALYSIS (present in analysis.sh) :

### CALCULATION OF THE INTENSITY OF RELATIONSHIPS WITH AN ACTIVE SITE

The src/python/allosteric_site_search.py script receives as input after the **-f** flag a json file with a list of residues in the active site, as well as the name of this list in the received file after the **-asn** flag. Output is made to the ./output/analysis/ directory. To run use:

> python src/python/allosteric_site_search.py -f input.json -n project_name -nclust num_of_clust -asn list_name

## OTHER USEFUL TOOLS (present in tools.sh) :

### COMPARISON WITH USER MATRIX

The calculation of the Frobenius norm of a given matrix and the user matrix and drawing the difference matrix of these normalized matrices is done using the src/python/matrix_comparison.py script, which receives a json file with the saved user matrix as input after the **-f** flag, as well as the name of the matrix in this file after the **-matname** flag, and after saves everything to the ./output/analysis/ directory. To run use:

> python src/python/matrix_comparison.py -f input.json -n project_name -nclust num_of_clust -matname matrix_name


