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

    # Count occurrences of gaps per sequence
    gap_percentages = np.sum(a == gaps, axis=1)
    # Get sorted indices by gap percentage
    sort = np.argsort(gap_percentages)
    # Sort the array by gap percentage
    a = a[sort]

    # Calculation of distance matrix
    # Initialize distance matrix
    d = np.zeros([numseq, numseq])
    for i in tqdm(range(numseq), desc="Calculating distance matrix for Hobohm clustering"):
        # Only consider non-gap positions
        valid_positions = a[i] != gaps
        d[i + 1:, i] = d[i, i + 1:] = 1 - np.sum(a[i] == a[i + 1:], axis=1) / float(np.sum(valid_positions))

    # Initialize cluster dictionary and unassigned sequences array
    cls = {}
    seqs = np.ones(numseq)

    for i in range(numseq):
        # If sequence is unassigned
        if seqs[i]:
            # Find matches which are unassigned and below the threshold
            match = np.nonzero(np.logical_and(d[i, i + 1:] < 1 - cl, seqs[i + 1:]))[0] + i + 1
            # Add matching sequences to cluster
            cls[i] = list(match) if match.size else []
            if match.size:
                # Mark the matched sequences as assigned
                seqs[match] = 0

    # Compute cluster weights
    cluster_weights = sorted([(x, 1.0 / len([k] + v)) for k, v in cls.items() for x in [k] + v])
    # Sort by the original order
    sort, weights = zip(*sorted(cluster_weights))
    weights = list(weights)
    weights = np.array(weights)
    print(f'Identified {len(cls)} clusters')
    return weights

# MI between two columns
def calc_MI(col_pair, matrix, bins=20, ignore_values=None, weights=None):
    col1, col2 = col_pair
    data1, data2 = matrix[:, col1], matrix[:, col2]

    # Filter out values to be ignored for histogram calculations
    mask = np.isin(data1, ignore_values, invert=True) & np.isin(data2, ignore_values, invert=True)
    filtered_data1 = data1[mask]
    filtered_data2 = data2[mask]

    # Use weights for histogram calculations if present
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

# Function to compute a MI matrix column-wise using parallel processing (CPU)
def comp_MI(matrix, bins=20, ignore_values=None, weights=None):
    if ignore_values is None:
        ignore_values = []
    num_cols = matrix.shape[1]
    col_pairs = [(i, j) for i in range(num_cols) for j in range(i + 1, num_cols)]

    # Use partial to pass the matrix, bins, values to ignore and weights to the worker function
    worker_func = partial(calc_MI, matrix=matrix, bins=bins, ignore_values=ignore_values, weights=weights)

    # Run
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

# Function to calculate MI from scrambled alignment (randomized positions of AA)
def null_samples(int_matrix, iterations, ign, bins, ignored, rcw, apc, weights):

    i = 0
    null_list = []
    pbar = tqdm(total=iterations)
    while i < iterations:
        shuffled = int_matrix[:, np.random.permutation(int_matrix.shape[1])]
        sample = comp_MI(shuffled, bins, ignore_values=ignored, weights=weights)
        if rcw:
            sample = apply_rcw(sample)
            sample = (sample - np.min(sample)) / (np.max(sample) - np.min(sample))
        if apc:
            sample = apply_apc(sample, ign)
        null_list.append(sample)
        i += 1
    pbar.close()
    null_list = np.array(null_list)
    return null_list


def zscore_transform_matrices_against_population(sample,  null):

    z_mat = np.zeros((len(sample), len(sample)))

    # Iterate over every element in the matrix
    for i in range(len(sample)):
        for j in range(len(sample)):
            if i != j:  # Skip diagonal elements
                # Collect the corresponding elements from all matrices
                corresponding_elements = []
                for element in null:
                    corresponding_elements.append(element[i, j])
                #corresponding_elements = np.array(null[k][i, j] for k in range(num_matrices))

                # Calculate mean and std of corresponding elements
                mean = np.mean(corresponding_elements)
                std = np.std(corresponding_elements)

                # Compute Z-score for the current element
                if std != 0:
                    z_score = (sample[i, j] - mean) / std
                else:
                    z_score = 0  # Avoid division by zero

                # Update the matrix with the Z-score
                z_mat[i, j] = z_score

    return z_mat


#  Function to calculate column-wise row column weighting
def calc_rcw(col_pair, matrix):
    col1, col2 = col_pair
    data1, data2 = matrix[:, col1], matrix[:, col2]
    # Mean MI in column i
    imean = data1.sum() / (len(data1) - 1)
    # Mean MI in column j
    jmean = data2.sum() / (len(data2) - 1)
    rcw = (imean + jmean) / 2
    return col1, col2, rcw

#  Function to apply row column weighting to MI matrix
def apply_rcw(matrix):
    num_cols = matrix.shape[1]
    col_pairs = [(i, j) for i in range(num_cols) for j in range(i + 1, num_cols)]

    # Use partial to pass the matrix to the worker function
    worker_func = partial(calc_rcw, matrix=matrix)

    # Run
    with Pool() as pool:
        rcw_results = pool.map(worker_func, col_pairs)

    # Create an empty matrix for rcw  values
    rcw_matrix = np.zeros((num_cols, num_cols))

    # Fill in the upper triangular part of the matrix
    for col1, col2, rcw in rcw_results:
        rcw_matrix[col1, col2] = rcw
        # Symmetric matrix
        rcw_matrix[col2, col1] = rcw
    np.seterr(invalid='ignore')
    corrected = np.divide(matrix, rcw_matrix)
    return corrected

#  Function to calculate column-wise average product correction
def calc_apc(col_pair, matrix):
    col1, col2 = col_pair
    data1, data2 = matrix[:, col1], matrix[:, col2]
    # Mean MI of MI matrix
    mimean = np.sum(np.tril(matrix, -1)) * 2 / (len(data1) * (len(data1) - 1))
    # Mean MI in column i
    imean = data1.sum() / (len(data1) - 1)
    # Mean MI in column j
    jmean = data2.sum() / (len(data2) - 1)
    apc = (imean * jmean) / mimean
    return col1, col2, apc

#  Function to apply average product correction to MI matrix
def apply_apc(matrix, ign):
    num_cols = matrix.shape[1]
    col_pairs = [(i, j) for i in range(num_cols) for j in range(i + 1, num_cols)]

    # Use partial to pass the matrix to the worker function
    worker_func = partial(calc_apc, matrix=matrix)

    # Run
    with Pool() as pool:
        apc_results = pool.map(worker_func, col_pairs)

    # Create an empty matrix for apc values
    apc_matrix = np.zeros((num_cols, num_cols))

    # Fill in the upper triangular part of the matrix
    for col1, col2, apc in apc_results:
        apc_matrix[col1, col2] = apc
        # Symmetric matrix
        apc_matrix[col2, col1] = apc
    corrected = np.subtract(matrix, apc_matrix)
    # Set negative values to 0
    if ign:
        corrected[corrected < 0] = 0
    return corrected

# Main function
def map_from_msa(name, out, format, apc, rcw, cl, igg, zs, igc, igp, ign):
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
    if cl:
        weights = hobohm(int_matrix, gap, cl)
    if igg:
        ignored.append(gap)
    ignored.append(place_holder_num)

    mimat = comp_MI(int_matrix, bins=bins, ignore_values=ignored, weights=weights)
    res = mimat

    if rcw:
        res = apply_rcw(res)
        print(f"Applied RCW.")
        res = np.nan_to_num(res)
        res = (res - np.min(res)) / (np.max(res) - np.min(res))
    if apc:
        res = apply_apc(res, ign)
        print(f"Applied APC.")
    if zs:
        samples = null_samples(res, zs, ign, bins, ignored, rcw, apc, weights)
        z_mat = zscore_transform_matrices_against_population(res, samples)
        res = res*np.abs(z_mat)

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
