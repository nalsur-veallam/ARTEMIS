NAME="v536e" # Project name
PYTHON="python3" # Your python launch codeword version >=3
NUM_OF_CLUST=5 # Desired number of clusters in the system
SOURCE_PDB="test_system/v536e.pdb" # Path to pdb file to create a pymol session with clustering (you must have pymol installed!!!)

${PYTHON} src/python/clustering.py -n ${NAME} -nclust ${NUM_OF_CLUST} # Calculates the clustering of the system and builds a dendrogram (in future versions I will add clustering by cutoff)

${PYTHON} src/python/create_pse.py -n ${NAME} -nclust ${NUM_OF_CLUST} -f ${SOURCE_PDB} # Just comment out the line if you don't want to get the session

# All results are stored in the ./output/clusterization/ directory
