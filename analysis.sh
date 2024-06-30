NAME="v536e" # Project name
PYTHON="python3" # Your python launch codeword version >=3
SOURCE_ACTIVE_SITE="test_system/v536e_groups.json" # Path to a json file with a list of active site amino acid residues in the PARENT order
                                                        # (can be understood from a json file obtained from a binary file using the "names" and "real _numbers" lists)
SOURCE_ALLOSTERIC_SITE="test_system/v536e_groups.json" # Path to a json file with a list of allosteric site amino acid residues in the PARENT order
SOURCE_PDB="test_system/v536e.pdb" # Path to pdb file with protein
ALLOSTERIC_SITE_NAME="allosteric_site" # List name in json file
ACTIVE_SITE_NAME="active_site" # List name in json file
OTHER_NAME="..."

mkdir output &> /dev/null
mkdir output/${NAME} &> /dev/null
mkdir output/${NAME}/analysis &> /dev/null

${PYTHON} src/python/allostery_search.py -n ${NAME} -f ${SOURCE_ACTIVE_SITE} -asn ${ACTIVE_SITE_NAME} -table -top 10 # Draws the intensity of association
#                                                         of amino acid residues with the active site (Set the -filt or -sasa_filt flag to use a filtered matrix)
#                                                                      (Add the -noseq flag to ignore the interaction of adjacent residues in sequence)
#                                                 to show the intensity in the form of a table, use the -table flag; to draw intensity in zscore format use the -zscore flag
#                                                           To take into account only the top X percent of the intensity when drawing, use -top X

${PYTHON} src/python/allostery_paint.py -strc ${SOURCE_PDB} -n ${NAME} -f_act ${SOURCE_ACTIVE_SITE} -asn ${ACTIVE_SITE_NAME} -allsn ${ALLOSTERIC_SITE_NAME} -f_all ${SOURCE_ALLOSTERIC_SITE} -top 10
#                  Draws the intensity of association of amino acid residues with the active site in a pymol session (Set the -filt or -sasa_filt flag to use a filtered matrix)
#                                                                      (Add the -noseq flag to ignore the interaction of adjacent residues in sequence)
#                                                           To take into account only the top X percent of the intensity when drawing, use -top X

${PYTHON} src/python/allostery_analysis.py -n ${NAME} -f_act ${SOURCE_ACTIVE_SITE} -asn ${ACTIVE_SITE_NAME} -allsn ${ALLOSTERIC_SITE_NAME} -f_all ${SOURCE_ALLOSTERIC_SITE} -zscore
#                                                           To take into account only the top X percent of the intensity when drawing, use -top X
#                                                             to analyze intensity in zscore format use the -zscore flag

#${PYTHON} src/python/compare_systems.py -n1 ${NAME} -n2 ${OTHER_NAME} # Script to compare two systems

# All results are stored in the ./output/${NAME}/analysis/ directory
# If you do not need to execute any of the programs, then just comment out the corresponding line
