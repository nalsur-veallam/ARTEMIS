## ARTEMIS MAP

Module for constructing and processing a matrix of mutual information based on PARENT data. To view existing programs in a module, run:

```bash
artemis map -h
```

The `artemis map` module contains four programs. Two of them (`--gen` and `--denoise`) calculate the mutual information matrix from PARENT data. The `--draw` program draws the mutual information matrix as a heatmap. The `--contacts` program draws the Contacts enrichment matrix as a heatmap.

#### `--gen` program

To calculate the raw mutual information matrix (that is, the matrix that contains noise) from PARENT data, use the `artemis map --gen` program. The input is a binary file PARENT in `.par` format. The program outputs an internal ARTEMIS mutual information matrix format in the form of a `.json` file. So:

```bash
usage: artemis map --gen -o O file [.par-file]

positional arguments:
 file .par-file

options:
 -o O Outfile path.
```
And an example of using the program:

```bash
artemis map example.par --gen -o out_path.json
```

#### `--denoise` program

To obtain a mutual information matrix with noise removed, use the program `artemis map --denoise`. Two noise removal algorithms are implemented: *linear* and *approximation*. ***We do not recommend using linear if you do not understand what exactly you need. The approximation method produces the best values. Use the approximation method if you don't understand what both methods do.***

Two binary PARENT files in .par format are submitted for input. Both files must contain data of the same MD trajectory, but with a different number of frames (that is, the same fixed-length trajectory, but with a different frame recording time interval. For example, the trajectory is 1us with the record dt1=1ps and dt2=2ps - that is, a million and half a million frems, respectively). This is done in order to determine the noise intensity through their difference. The `-lin` flag is used to select the *linear* method (otherwise, the *approximation* method is used). The number of frames in the trajectories by which the files were counted .par is supplied via variables with the `-n1` and `-n2` flags.  You can also specify the number of frames of the approximate trajectory using the `-n0` flag (infinity by default. ***Don't use this if you don't understand how it works!***). At the output, the program outputs the internal format of the matrix of mutual information for ARTEMIS in the form of a `.json` file. So:

```bash
usage: artemis map --denoise -n0 N0 -n1 N1 -n2 N2 -o O -lin files [.par-files ...]

positional arguments:
 files .par-files

options:
 -n0 N0 Denoised map number of frames (default is inf).
 -n1 N1 First custom number of frames.
 -n2 N2 Second custom number of frames.
 -o O Outfile path.
```

And an example of using the program:

```bash
artemis map example1.par example2.par --denoise -n0 N0 -n1 N1 -n2 N2 -o out_path.json
```
#### `--draw` program

To draw the mutual information matrix, use the `artemis map --draw` program. As input, the program receives the ARTEMIS internal mutual information matrix format in the form of a `.json` file. The output is an image in the format specified by the user from those supported by python (default `.pdf`). Using the `-norm` flag, it is possible to normalize the matrix so that the maximum element of the matrix is ​​equal to 1. Also, the matrix is ​​drawn without diagonal values. If you want to draw them, use the `-diag` flag. So:

```bash
usage: artemis map --draw -norm -diag -o O file [.json-file]

positional arguments:
 files .json-file

options:
 -o O Outfile path.
 -norm       Normalize the MIE matrix before drawing.
 -diag       Draw the diagonal of the MIE matrix.
```

And an example of using the program:

```bash
artemis map example.json --draw -norm -diag -o out_path.pdf
```

#### `--contacts` program

To draw the Contacts enrichment matrix, use the `artemis map --contacts` program. As input, the program receives the ARTEMIS internal mutual information matrix format in the form of a `.json` file. The output is an image in the format specified by the user from those supported by python (default `.pdf`). Using the `-norm` flag, it is possible to zero out all pairs of residues for which enrichment is less than one. Also, the matrix is ​​drawn without diagonal values. If you want to draw them, use the `-diag` flag. To specify the maximum value of colorbar (this can be useful for comparing matrices of different systems), you can specify it using the `-vmax` flag. So:

```bash
usage: artemis map --contacts -norm -diag -vmax VMAX -o O file [.json-file]

positional arguments:
 files .json-file

options:
 -o O Outfile path.
 -norm       Normalize the MIE matrix before drawing.
 -diag       Draw the diagonal of the MIE matrix.
 -vmax VMAX  The maximum value of colorbar for Contacts map.
```

And an example of using the program:

```bash
artemis map example.json --contacts -norm -diag -vmax 2 -o out_path.pdf
```

***To see an example of use, look at the tutorial for using the program at the [link](https://nalsur-veallam.github.io/ARTEMIS/tutorial.html).***
