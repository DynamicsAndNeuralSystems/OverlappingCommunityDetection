# Empirical comparison of overlapping community detection algorithms

This package facilitates the computation of a wide range of overlapping community detection algorithms (OCDAs), and visualization of the results.

For the implementation used in our preprint ['Benchmarking overlapping community detection methods for applications in human connectomics'](https://www.biorxiv.org/content/10.1101/2025.03.19.643839v1.abstract), please refer to the designated repository: [https://github.com/DynamicsAndNeuralSystems/OverlappingCommunityDetection_HCP](https://github.com/DynamicsAndNeuralSystems/OverlappingCommunityDetection_HCP).

## Setup

Run `startup` (from any directory) to add the toolbox to the MATLAB path.

Methods implemented in MATLAB (SLPA, Clique, NNMF, WSBM) need nothing else.
The other methods call external programs, which are not included as binaries and must be compiled for your system:

1. **OSLOM**: in a terminal, navigate to `SourceCode/OSLOM` and run `./compile_all.sh`.
   Only `oslom_undir` is needed; errors compiling the other programs (pajek, louvain, infomap) can be ignored.
   See [the OSLOM website](http://www.oslom.org/software.htm) for more information.
2. **Infomap**: navigate to `SourceCode/Infomap` and run `make clean` and then `make`.
3. **SVINET** (Gopalan): requires the [GSL](https://www.gnu.org/software/gsl/) library.
   Navigate to `SourceCode/Gopalan` and run `./configure` and then `make` (see `SourceCode/Gopalan/README.md` for setting library paths on macOS).
4. **Link communities**: requires Python 3.
   MATLAB calls `python3` by default; to use a different interpreter, set the environment variable `OCDA_PYTHON` (e.g., `setenv('OCDA_PYTHON','/opt/homebrew/anaconda3/bin/python3')`).
5. **WSBM** (optional speed-up): run `InstallMEXFiles` in `SourceCode/Clauset` to compile its MEX files.

## Computation

### Usage

```matlab
[Final,computationTimes] = OLCD_Compute(networkAdj,methodList,isBenchmark,benchmarkFileName,sourceCodePath,outputPath);
```

Inputs (all but the first are optional):

- `networkAdj`: an adjacency matrix, or an edge list as a 3-column matrix `[node1, node2, weight]` with positive-integer node ids.
  An edge list is treated as undirected if `node1 <= node2` in every row, and as directed otherwise.
  A 3x3 input is ambiguous and is treated as an adjacency matrix (with a warning).
- `methodList`: cell array of method names (see below).
- `isBenchmark`: `true` if the network has known (benchmark) communities (default `false`).
- `benchmarkFileName`: text file defining the benchmark communities, one line per node: `node comm1 comm2 ...` (as produced by the LFR benchmark generator).
- `sourceCodePath`: path to the `SourceCode` directory (default: `SourceCode` in this toolbox).
- `outputPath`: directory in which results are saved (default: the current directory).

Outputs:

- `Final` (structure) contains the network, the benchmark communities (if provided), and one field per method output.
  Each method output is a structure with the method's `Name`, its parameters, and its `Result`.
- `computationTimes` (structure) gives the computation time (in seconds) of each method output.
- Results are also saved to `Computation_Result.mat` in `outputPath`, which can be visualized using code in `Visualization`.

### Methods

| Name in `methodList` (alias) | Method | Parameters swept by `OLCD_Compute` | Output field names |
| --- | --- | --- | --- |
| `'Infomap'` | Infomap (with `--overlapping`) | – | `Infomap` |
| `'OSLOM'` | Order Statistics Local Optimization Method | p-value tolerance 0.1, 0.2, …, 0.9 | `OSLOM_10`, …, `OSLOM_90` |
| `'NNMF'` | Non-negative matrix factorization | membership threshold 0.1–0.4 | `NNMF_10`, …, `NNMF_40` |
| `'SLPA'` (`'Jerry'`) | Speaker–listener label propagation | – | `Jerry` |
| `'Clique'` (`'Shen'`) | Clique-based method of Shen et al. (2009) | clique size 3, 4, 5, 6, 7, 9 | `Clique_3`, …, `Clique_9` |
| `'Link'` | Link communities (Ahn et al., 2010) | proportion of weakest links removed: 0, 0.01, 0.1 | `Link_prec_0`, … |
| `'SVINET'` (`'Gopalan'`) | Stochastic variational inference for the MMSB (Gopalan & Blei, 2013) | proportion of weakest links removed: 0, 0.1, 0.2, 0.3 | `Gopalan_prec_0`, … |
| `'WSBM'` (`'Clauset'`) | Weighted stochastic block model (Aicher et al., 2014), number of blocks chosen by maximum log-evidence over k = 2–8. Gives a *non-overlapping* partition. | – | `Clauset` |

### Algorithm output

Each method's `Result` is a matrix in which each row is a node and each column is a community.
Each node's memberships sum to 1, split equally across the communities it belongs to (or 0 if the method left the node unassigned).
Empty communities are removed.

### Structure

Each community detection method has its own directory in `Computation`, which contains three files:

- `call_[MethodName]`: called by `OLCD_Compute`; it runs the method and processes its output in one line.
- `run_[MethodName]`: runs the algorithm (for external programs, input and output files are written to a temporary directory).
- `process_[MethodName]`: converts the algorithm's output into the standard matrix format described above.

In addition:

- Source code for each algorithm is in a directory in `SourceCode`.
- Functions for converting between adjacency matrices and edge lists are in `Computation/Conversions`.
- Shared utilities (temporary directories, running external commands, reading community files, tidying the output matrix) are in `Computation/Helpers`.

### Adding a new method

1. Place the source code of the algorithm in its own directory inside `SourceCode`.
2. Write `run_[Method]`: this should take the network (in whatever format is needed) and parameters, and output a community result (in any format).
3. Write `process_[Method]`: this takes the output of the previous step and converts it to the standard matrix format (the helpers `readCommunityLines`, `communitiesToMatrix`, and `tidyCommunities` may be useful).
4. Write `call_[Method]`: this calls the previous two functions, so the method can be run in one line.
5. Add a case for the method in the `switch` statement of `OLCD_Compute` (looping over parameter values if needed), and add its name to the list of valid method names.

## Visualization

### Usage

```matlab
OLCD_Visualization(Methods,resultsFile,cmap)
```

- `Methods`: cell array of methods to plot, matched either by method name (e.g., `'NNMF'` plots all NNMF thresholds) or by output field name (e.g., `'NNMF_20'`).
  Include `'Benchmark'` to plot the benchmark communities.
  Default: all outputs in the results file.
- `resultsFile`: results file produced by `OLCD_Compute` (default `'Computation_Result.mat'`).
- `cmap`: colormap name, function handle, or matrix (default `'parula'`).

Nodes are ordered by the benchmark communities (if present), and otherwise by the first method plotted.

### Files

- `OLCD_Visualization`: plots the network adjacency matrix alongside the communities of each method, for comparison.
- `Node_Reorder`: orders each method's communities from largest to smallest (colors depend on the community labels).
- `Node_Sorter`: orders the nodes according to the benchmark communities, so that communities are easier to see.

### Comparing communities

- `ENMI_calc`: the extended normalized mutual information for overlapping communities (Lancichinetti, Fortunato & Kertész, 2009), between two community matrices (uses `extendNMI_calcs`).
- `NMI_calc`: normalized mutual information between two non-overlapping community matrices (errors on overlapping input).
- `NMI_calc_vecs`: normalized mutual information between two vectors of community labels.
- `F1_overlapcalc`: extracts the set of overlapping nodes from two community matrices and computes the F1 score for recovering them (using `F1_calc`).
