import numpy as np
import json
from Bio import AlignIO
import time
from multiprocessing import Pool
from functools import partial
from tqdm import tqdm

# Convert the alignment to a numpy array
def al_to_mat(alignment, igc):
    mat = []
    length = alignment.get_alignment_length()
    i = 0
    if igc:
        while i < length:
            a = (np.array(list((str(alignment[:, i]).upper()))))
            mat.append(a)
            i += 1
        mat = np.transpose(mat)
    else:
        while i < length:
            a = (np.array(list((alignment[:, i]))))
            mat.append(a)
            i += 1
        mat = np.transpose(mat)
    return np.array(mat)


# Convert the array of characters to numerical msa
def mat_to_nmsa(matrix, igp, igg):
    place_holder_num = None
    unique_chars = list(np.unique(matrix))
    # Specify number of bins depending on unique characters and specified ignored symbols
    bins = len(unique_chars)
    # Convert to numercial
    char_to_int = {char: idx for idx, char in enumerate(unique_chars)}

    # Store place holder number and remove one bin
    if igp:
        place_holder_num = char_to_int[igp]
        bins -= 1
        unique_chars.remove(str(igp))

    # Store gap number and remove one bin
    gap = char_to_int['-']
    if igg:
        bins -= 1
        unique_chars.remove('-')
    print("These characters will be treated as unique in the MSA")
    print(unique_chars)
    int_matrix = np.vectorize(char_to_int.get)(matrix)
    return int_matrix, bins, gap, place_holder_num


def hobohm(matrix, gaps, cl):
    a = np.array(matrix)
    numseq, seqlen = a.shape

    gap_percentages = np.sum(a == gaps, axis=1)  # Count occurrences of ignore_value per sequence
    sort = np.argsort(gap_percentages)  # Get sorted indices by gap percentage
    a = a[sort]  # Sort the array by gap percentage

    # CALCULATE DISTANCE MATRIX, based on differences not including gaps
    d = np.zeros([numseq, numseq])  # Initialize distance matrix
    for i in tqdm(range(numseq), desc="Calculating distance matrix for Hobohm clustering"):
        valid_positions = a[i] != gaps  # Only consider non-gap positions
        d[i + 1:, i] = d[i, i + 1:] = 1 - np.sum(a[i] == a[i + 1:], axis=1) / float(np.sum(valid_positions))

    # APPOINT SEQUENCES TO CLUSTERS DEPENDING ON THEIR DISTANCES
    c = {}  # Dictionary for clusters
    s = np.ones(numseq)  # Array indicating unassigned sequences

    for i in range(numseq):
        if s[i]:  # If sequence is unassigned
            # Find matches which are unassigned and below the threshold
            m = np.nonzero(np.logical_and(d[i, i + 1:] < 1 - cl, s[i + 1:]))[0] + i + 1
            c[i] = list(m) if m.size else []  # Add matching sequences to cluster
            if m.size:
                s[m] = 0  # Mark the matched sequences as assigned

    # COMPUTING CLUSTER WEIGHTS
    cluster_weights = sorted([(x, 1.0 / len([k] + v)) for k, v in c.items() for x in [k] + v])
    sort, w = zip(*sorted(cluster_weights))  # Sort by the original order
    w = list(w)
    w = np.array(w)
    print(f'Identified {len(c)} clusters')
    return w

# MI between two columns
def calc_MI(col_pair, matrix, bins=20, ignore_values=None, weights=None):
    col1, col2 = col_pair
    data1, data2 = matrix[:, col1], matrix[:, col2]

    # Filter out values to be ignored for histogram calculations
    mask = np.isin(data1, ignore_values, invert=True) & np.isin(data2, ignore_values, invert=True)
    filtered_data1 = data1[mask]
    filtered_data2 = data2[mask]

    if weights is not None:
        filtered_weights = weights[mask]
    else:
        filtered_weights = np.ones_like(filtered_data1)

    # Compute histograms
    hist_2d, x_edges, y_edges = np.histogram2d(filtered_data1, filtered_data2, bins=bins, weights=filtered_weights)

    # Normalize the histogram to probabilities
    p_xy = hist_2d / np.sum(hist_2d)
    p_x = np.sum(p_xy, axis=1)  # marginal for x
    p_y = np.sum(p_xy, axis=0)  # marginal for y

    # Only keep non-zero entries to avoid log(0)
    nz_x, nz_y = np.nonzero(p_xy)
    p_xy_nonzero = p_xy[nz_x, nz_y]
    p_x_nonzero = p_x[nz_x]
    p_y_nonzero = p_y[nz_y]

    # Mutual information
    mi = np.sum(p_xy_nonzero * np.log(p_xy_nonzero / (p_x_nonzero * p_y_nonzero)))

    return col1, col2, mi

def comp_MI(matrix, bins=20, ignore_values=None, weights=None):
    if ignore_values is None:
        ignore_values = []
    num_cols = matrix.shape[1]
    col_pairs = [(i, j) for i in range(num_cols) for j in range(i + 1, num_cols)]

    # Use partial to pass the matrix, bins, and values to ignore to the worker function
    worker_func = partial(calc_MI, matrix=matrix, bins=bins, ignore_values=ignore_values, weights=weights)

    # Run in parallel
    with Pool() as pool:
        mutual_info_results = pool.map(worker_func, col_pairs)

    # Create an empty matrix for mutual information values
    mi_matrix = np.zeros((num_cols, num_cols))

    # Fill in the upper triangular part of the matrix
    for col1, col2, mi in mutual_info_results:
        mi_matrix[col1, col2] = mi
        # Symmetric matrix
        mi_matrix[col2, col1] = mi

    return mi_matrix


def remove_bg(mimat, iterations, ign, bins, ignored, int_matrix):
    # i = 0
    # num_cols = mimat.shape[1]
    # bg_matrix = np.zeros((num_cols, num_cols))
    # while i < iterations:
    #     shuffled = np.apply_along_axis(np.random.permutation, axis=1, arr=int_matrix)
    #     bg_matrix = np.add(bg_matrix, comp_MI(shuffled, bins, ignore_values=ignored))
    #     i += 1
    # bg_matrix = np.divide(bg_matrix, iterations)
    # corrected = np.subtract(mimat, bg_matrix)
    i = 0
    corrected = mimat
    while i < iterations:
        shuffled = np.apply_along_axis(np.random.permutation, axis=1, arr=int_matrix)
        bgmat = comp_MI(shuffled, bins, ignore_values=ignored)
        corrected = np.subtract(corrected, bgmat)
        i += 1
    if ign:
        corrected[corrected < 0] = 0
    return corrected


def calc_rcw(col_pair, matrix):
    col1, col2 = col_pair
    data1, data2 = matrix[:, col1], matrix[:, col2]
    imean = data1.sum() / (len(data1) - 1)
    jmean = data2.sum() / (len(data2) - 1)
    rcw = (imean + jmean) / 2
    return col1, col2, rcw


def apply_rcw(matrix):
    num_cols = matrix.shape[1]
    col_pairs = [(i, j) for i in range(num_cols) for j in range(i + 1, num_cols)]

    # Use partial to pass the matrix to the worker function
    worker_func = partial(calc_rcw, matrix=matrix)

    # Run in parallel
    with Pool() as pool:
        rcw_results = pool.map(worker_func, col_pairs)

    # Create an empty matrix for apc  values
    rcw_matrix = np.zeros((num_cols, num_cols))

    # Fill in the upper triangular part of the matrix
    for col1, col2, rcw in rcw_results:
        rcw_matrix[col1, col2] = rcw
        # Symmetric matrix
        rcw_matrix[col2, col1] = rcw
    np.seterr(invalid='ignore')
    corrected = np.divide(matrix, rcw_matrix)
    return corrected


def calc_apc(col_pair, matrix):
    col1, col2 = col_pair
    data1, data2 = matrix[:, col1], matrix[:, col2]
    mimean = np.sum(np.tril(matrix, -1)) * 2 / (len(data1) * (len(data1) - 1))
    imean = data1.sum() / (len(data1) - 1)
    jmean = data2.sum() / (len(data2) - 1)
    apc = (imean * jmean) / mimean
    return col1, col2, apc


def apply_apc(matrix, ign):
    num_cols = matrix.shape[1]
    col_pairs = [(i, j) for i in range(num_cols) for j in range(i + 1, num_cols)]

    # Use partial to pass the matrix to the worker function
    worker_func = partial(calc_apc, matrix=matrix)

    # Run in parallel
    with Pool() as pool:
        apc_results = pool.map(worker_func, col_pairs)

    # Create an empty matrix for apc  values
    apc_matrix = np.zeros((num_cols, num_cols))

    # Fill in the upper triangular part of the matrix
    for col1, col2, apc in apc_results:
        apc_matrix[col1, col2] = apc
        # Symmetric matrix
        apc_matrix[col2, col1] = apc
    corrected = np.subtract(matrix, apc_matrix)
    if ign:
        corrected[corrected < 0] = 0
    return corrected


def map_from_msa(name, out, format, apc, rcw, cl, igg, ss, igc, igp, ign):
    start = time.time()
    # Read alignment
    alignment = AlignIO.read(open(name), format)
    # Print information about alignment
    print("Alignment length %i" % alignment.get_alignment_length())
    print("Number of sequences %i" % len(alignment))
    ignored=[]
    matrix = al_to_mat(alignment, igc)
    int_matrix, bins, gap, place_holder_num = mat_to_nmsa(matrix, igp, igg)

    weights = None
    if cl:# or cl == float(0):
        weights = hobohm(int_matrix, gap, cl)

    if igg:
        ignored.append(gap)
    ignored.append(place_holder_num)

    mimat = comp_MI(int_matrix, bins=bins, ignore_values=ignored, weights=weights)
    res = mimat
    if ss:
        res = remove_bg(res, ss, ign, bins, ignored, int_matrix)
        print(f"Removed {ss} iterations of randomized alignment MI.")
    if rcw:
        res = apply_rcw(res)
        print(f"Applied RCW.")
        res = np.nan_to_num(res)
        res = (res - np.min(res)) / (np.max(res) - np.min(res))

    if apc:
        res = apply_apc(res, ign)
        print(f"Applied APC.")
    # Output
    res = res.tolist()
    data = {}
    data['NResidues'] = len(res)
    data['map'] = res
    data['names'] = list(alignment[0,:]) # list('R' * len(res))
    data['real_numbers'] = list(range(1, (len(res) + 1)))
    with open(f'{out + ".json"}', "w") as f:
        json.dump(data, f)
    end = time.time()
    print(f'Writing {out + ".json"}')
    print(f"Finished in {end - start} seconds.")
