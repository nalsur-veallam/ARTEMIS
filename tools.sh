NAME="v536e" # Project name
PYTHON="python3" # Your python launch codeword version >=3

SOURCE_FIRST_MAP="test_system/v536e_map.json" # Path to first matrix for comparison
SOURCE_SECOND_MAP="test_system/dist_map.json" # Path to second matrix for comparison

SOURCE_FIRST_INTEN="test_system/v536e_intensity.json" # Path to first json file with intensities for comparison
SOURCE_SECOND_INTEN="test_system/v536e_entropy.json" # Path to second json file with intensities for comparison
FIRST_INTEN_NAME="intensity"
SECOND_INTEN_NAME="entropy"

SOURCE_AF="..." # Path to json file with alpha fold data for draw it and convert to our format
SOURCE_2D="test_system/2dproj.xvg" # Path to xvg gromacs file with with the projection of the trajectory onto two eigenvectors
SOURCE_XPM="test_system/dm.xpm" # Path to xmp gromacs matrix file for draw it and convert to our format
SOURCE_PDB="test_system/v536e.pdb" # Path to pdb file with protein and ligand
SOURCE_PSE="test_system/v536e.pse" # Path to pse file with selection

SOURCE_GROUP="test_system/v536e_groups.json" # Path to json with group to make ndx gromacs file
SOURCE_PDB_CLEAN="test_system/v536e.pdb" # Path to pdb file with protein
GROUP_NAME="active_site" # List name in json file

DIST_OUT="output/other/dist_map"
MAP_COMP_OUT="output/other/comp_map"
TD_PROJ_OUT="output/other/2dproj"
NDX_OUT="output/other/act_site"
INTEN_COMP_OUT="output/other/comp_inten"

mkdir output &> /dev/null
mkdir output/other/ &> /dev/null

${PYTHON} src/python/matrix_comparison.py -f1 ${SOURCE_FIRST_MAP} -f2 ${SOURCE_SECOND_MAP} -o ${MAP_COMP_OUT} -nodiag  # Calculates the Frobenius norm of matrices and draws the difference matrix 
#                                                                                                   of the normalized matrices (You can zero out the diagonal with the -nodiag flag)

${PYTHON} src/python/compare_intensities.py -f1 ${SOURCE_FIRST_INTEN} -f2 ${SOURCE_SECOND_INTEN} -o ${INTEN_COMP_OUT} # Calculates the Frobenius norm of intensities and draws the difference 
#                                                                                                   intensity

${PYTHON} src/python/jens_shann_dist.py -f1 ${SOURCE_FIRST_INTEN} -f2 ${SOURCE_SECOND_INTEN} -hist ${NAME} -arrn1 ${FIRST_INTEN_NAME} -arrn2 ${SECOND_INTEN_NAME} -prob
#                                                   Calculates Jensen-Shannon distance (metric) between two probability arrays and draws the histograms                                                                                   

#${PYTHON} src/python/alf2json.py -f ${SOURCE_AF} -o alfmap # Draw alpha fold data matrix and convert it to our json format

${PYTHON} src/python/2dplot.py -f ${SOURCE_2D} -o ${TD_PROJ_OUT} # Draw projection of the trajectory onto two eigenvectors

${PYTHON} src/python/frob_from_json.py -f ${SOURCE_FIRST_MAP} -matname map # Calculates the Frobenius norm of a matrix

${PYTHON} src/python/make_ndx.py -f ${SOURCE_GROUP} -strc ${SOURCE_PDB_CLEAN} -o ${NDX_OUT} -grn ${GROUP_NAME} # Create gromacs ndx file from custom group

${PYTHON} src/python/xpm2json.py -f ${SOURCE_XPM} -o ${DIST_OUT} # Draw xmp gromacs matrix and convert it to our json format

#${PYTHON} src/python/find_area.py -f ${SOURCE_GROUP} -sn ${GROUP_NAME} -chain A -ligname GLC -strc ${SOURCE_PDB} -cutoff 3.5 # Creates a list of residues that are in contact 
#            with the ligand (to use the protein-protein bond, instead of the -ligname flag, specify the -chain2 flag followed by the name of the second chain (ligand protein))

${PYTHON} src/python/create_group.py -f ${SOURCE_GROUP} -ps ${SOURCE_PSE} -sel ${GROUP_NAME} # Calculates group in JSON file from pymol session by selection.

# If you do not need to execute any of the programs, then just comment out the corresponding line
