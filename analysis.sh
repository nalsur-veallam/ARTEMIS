NAME="v536e" # Project name
PYTHON="python3" # Your python launch codeword version >=3
NUM_OF_CLUST=5 # Desired number of clusters in the system
SOURCE_ACTIVE_SITE="test_system/v536e_active_site.json" # Path to a json file with a list of active site amino acid residues in the PARENT order 
                                                        # (can be understood from a json file obtained from a binary file using the "names" and "real _numbers" lists)
ACTIVE_SITE_NAME="active_site" # List name in json file
SOURCE_CUSTOM_MAP="test_system/v536e_map.json" # Path to additional matrix for comparison
MAT_NAME="map" # Matrix name in json file

mkdir output
mkdir output/analysis

${PYTHON} src/python/matrix_comparison.py -n ${NAME} -f ${SOURCE_CUSTOM_MAP} -matname ${MAT_NAME}  # Calculates the Frobenius norm of matrices and draws the difference matrix 
                                                                            # of the normalized matrices in the ./output/analysis/ directory
${PYTHON} src/python/allosteric_site_search.py -n ${NAME} -f ${SOURCE_ACTIVE_SITE} -asn ${ACTIVE_SITE_NAME} # Draws in the mao directory the intensity of association of amino acid residues with the active site

# If you do not need to execute any of the programs, then just comment out the corresponding line
