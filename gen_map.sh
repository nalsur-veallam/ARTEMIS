NAME="v536e" # Project name
SOURCE_PAR="test_system/v536e.par" # Path to binary file PARENT
PYTHON="python3" # Your python launch codeword version >=3

make clean
make

mkdir output &> /dev/null
mkdir output/map/ &> /dev/null
mkdir output/opt_num_of_clust/ &> /dev/null

bin/get_map -f ${SOURCE_PAR} -n ${NAME} # Obtaining a matrix of mutual information between residuals from a binary file
${PYTHON} src/python/draw_map.py -n ${NAME} # A matrix of mutual information on the remains is drawn in the directory ./output/map/
${PYTHON} src/python/opt_num_of_clust.py -n ${NAME} -max 10 # Graphs of the metric are drawn to obtain the optimal number of clusters in the directory ./output/opt_num_of_clust/
# If you want to change the range in which the metric graph is built, then add the -min and -max flags with the corresponding minimum and maximum number of clusters for the system
# Example: ${PYTHON} src/python/opt_num_of_clust.py -n ${NAME} -min 2 -max 10
