PYTHON="python3" # Your python launch codeword version >=3

SOURCE_FIRST_MAP="output/glu1/map/glu1_map.json" # Path to first matrix for comparison
SOURCE_SECOND_MAP="output/glu2/map/glu2_map.json" # Path to second matrix for comparison

SOURCE_FIRST_INTEN="output/glu1/analysis/glu1_intensity.json" # Path to first json file with intensities for comparison
SOURCE_SECOND_INTEN="output/glu2/analysis/glu2_intensity.json" # Path to second json file with intensities for comparison

SOURCE_AF="../glucokinase/analysis/Glu_c370c_unrelaxed_rank_5_model_1_scores.json" # Path to json file with alpha fold data for draw it and convert to our format
SOURCE_2D="../glucokinase/analysis/2d_1.xvg" # Path to xvg gromacs file with with the projection of the trajectory onto two eigenvectors
SOURCE_XPM="../glucokinase/analysis/dm_1.xpm" # Path to xmp gromacs matrix file for draw it and convert to our format
SOURCE_PDB="test_system/1v4s.pdb" # Path to pdb file with protein and ligand

SOURCE_GROUP="test_system/glu_groups.json" # Path to json with group to make ndx gromacs file
SOURCE_PDB_CLEAN="test_system/1v4s_clean.pdb" # Path to pdb file with protein
GROUP_NAME="act_s" # List name in json file
SOURCE_ACTIVE_SITE="test_system/glu_groups.json"

${PYTHON} src/python/matrix_comparison.py -f1 ${SOURCE_FIRST_MAP} -f2 ${SOURCE_SECOND_MAP} -o output/glu1/map/glu1-2_comp -nodiag  # Calculates the Frobenius norm of matrices and draws the difference matrix 
#                                                                                                   of the normalized matrices (You can zero out the diagonal with the -nodiag flag)

${PYTHON} src/python/compare_intensities.py -f1 ${SOURCE_FIRST_INTEN} -f2 ${SOURCE_SECOND_INTEN} -o inten_comp # Calculates the Frobenius norm of intensities and draws the difference 
#                                                                                                   intensity

${PYTHON} src/python/jens_shann_dist.py -f1 ${SOURCE_FIRST_INTEN} -f2 ${SOURCE_SECOND_INTEN} -hist glu -arrn1 intensity -arrn2 intensity -prob
#                                                   Calculates Jensen-Shannon distance (metric) between two probability arrays and draws the histograms                                                                                   

${PYTHON} src/python/alf2json.py -f ${SOURCE_AF} -o alfmap # Draw alpha fold data matrix and convert it to our json format

${PYTHON} src/python/2dplot.py -f ${SOURCE_2D} -o 2dproj # Draw projection of the trajectory onto two eigenvectors

${PYTHON} src/python/frob_from_json.py -f ${SOURCE_FIRST_MAP} -matname map # Calculates the Frobenius norm of a matrix

${PYTHON} src/python/make_ndx.py -f ${SOURCE_GROUP} -strc ${SOURCE_PDB_CLEAN} -o glu_act_site -grn ${GROUP_NAME} # Create gromacs ndx file from custom group

${PYTHON} src/python/xpm2json.py -f ${SOURCE_XPM} -o xpmmap # Draw xmp gromacs matrix and convert it to our json format

${PYTHON} src/python/find_area.py -f ${SOURCE_ACTIVE_SITE} -sn act_s -chain A -ligname GLC -strc ${SOURCE_PDB} -cutoff 3.5 # Creates a list of residues that are in contact 
#            with the ligand (to use the protein-protein bond, instead of the -ligname flag, specify the -chain2 flag followed by the name of the second chain (ligand protein))

#${PYTHON} src/python/create_group.py -f output.json -ps source.pse -sel selection_name # Calculates group in JSON file from pymol session by selection.

# If you do not need to execute any of the programs, then just comment out the corresponding line
