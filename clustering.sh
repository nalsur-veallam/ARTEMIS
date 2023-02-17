NAME="glu1" # Project name
PYTHON="python3" # Your python launch codeword version >=3
NUM_OF_CLUST=3 # Desired number of clusters in the system
SOURCE_PDB="test_system/1v4s_clean.pdb" # Path to pdb file to create a pymol session with clustering (you must have pymol installed!!!)
SOURCE_ACTIVE_SITE="test_system/glu_groups.json" # Path to a json file with a list of active site amino acid residues in the PARENT order 
                                                        # (can be understood from a json file obtained from a binary file using the "names" and "real _numbers" lists)
SOURCE_ALLOSTERIC_SITE="test_system/glu_groups.json" # Path to a json file with a list of allosteric site amino acid residues in the PARENT order
ALLOSTERIC_SITE_NAME="all_s" # List name in json file
ACTIVE_SITE_NAME="act_s" # List name in json file

mkdir output &> /dev/null
mkdir output/${NAME} &> /dev/null
mkdir output/${NAME}/clustering/ &> /dev/null

${PYTHON} src/python/opt_num_of_clust.py -n ${NAME} -max 20 # Graphs of the metric are drawn to obtain the optimal number of clusters in the directory ./output/opt_num_of_clust/
# If you want to change the range in which the metric graph is built, then add the -min and -max flags with the corresponding minimum and maximum number of clusters for the system
# Example: ${PYTHON} src/python/opt_num_of_clust.py -n ${NAME} -min 2 -max 10

${PYTHON} src/python/clustering.py -n ${NAME} -nclust ${NUM_OF_CLUST} # Calculates the clustering of the system and builds a dendrogram (in future versions I will add clustering by cutoff)

${PYTHON} src/python/create_pse.py -n ${NAME} -nclust ${NUM_OF_CLUST} -f ${SOURCE_PDB} # Just comment out the line if you don't want to get the session

${PYTHON} src/python/cluster_analysis.py -nclust ${NUM_OF_CLUST} -n ${NAME} -f_act ${SOURCE_ACTIVE_SITE} -asn ${ACTIVE_SITE_NAME} -allsn ${ALLOSTERIC_SITE_NAME} -f_all ${SOURCE_ALLOSTERIC_SITE} # Performs cluster analysis. If one of the sites is unknown, then just specify an empty list in the source (it will be fixed later)
# For cluster analysis, the -noseq parameter is also available, which specifies at what distance in the sequence to ignore the interaction (example: -noseq 1)

# All results are stored in the ./output/${NAME}/clustering/ directory
# If you do not need to execute any of the programs, then just comment out the corresponding line
