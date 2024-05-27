# Empirical comparison of overlapping community detection algorithms

This package facilitates the computation of a wide range of overlapping community detection algorithms, and visualization of  the results.

## Compilation

In the `SourceCode` directory, some algorithms need to be compiled from the terminal for your system.

1. In a terminal, navigate to `SourceCode/OSLOM` and run `./compile_all.sh', and then check that the result works. See [here](http://www.oslom.org/software.htm) for more information.
2. Now, navigate to `Source/Infomap`
   1. Run `make clean` and execute and terminal suggested code.
   2. Now run `make` 

## Computation

### Usage/Syntax:

`Computation` computes a range of community-detection algorithms on a given network, as:

```matlab
Final = OLCD_Compute(networkAdj,Methods,isBenchmark,benchFileName);
```

INPUTS:

- `networkAdj` -	Can be an adjacency matrix, (matrix) list of undirected links, or (matrix) list of directed links.
- `Methods` -	Cells of the names of methods to be computed.
- `isBenchmark` - Logical: `false` for no benchmark (default), or `true` for benchmark.
- `benchFileName` - String name of the text file within which the benchmark communities are defined

OUTPUTS:

- `Final` (structure) contains all data from the algorithms (Method name, community result).
- Results are also saved to `Computation_Result.mat`, which can be visualized using code in the Visualization module.

### Structure

Each community detection method has its own directory within the `Computation_Module` directory.
Each directory, `MethodName`, contains three files:

- `call_[MethodName]` - This function is called by the overarching `computation.m` function. It simplifies the entire process down to one line.
- `run_[MethodName]` - This function is called to do the actual calculations; running the algorithm and spitting out its output (Note: If calculations need to be hard coded to specific directories, this is done within this `.m` file.
  For examples, view `run_Gopalan`).
- `process_[MethodName]` - Called after running the algorithm, as this function will take the algorithm output and convert it to the universal matrix format, explained below.

In addition, source code for each algorithm is in a directory in the `SourceCode` directory.
And scripts for converting between the three types of inputs are in `Conversions`.

### Algorithm Output

After the `process_[Method]` function is run, the output should be in a matrix format, where each row is a node, and each column is a community.

### Adding a new method

If you wish to add a module, here is how to do so:

1. Place source code of the algorithm within its own named directory inside  the `Computation_Module/SourceCode` directory.
2. Code `run_[Method]` - This should take the network (in whatever format that is needed) and parameters, and output a final community result.
   This output does not have to be in the matrix format, the processing function should do that.
3. Code `process_[Method]` - This will take the output of the previous step as an input, and then output the matrix format for the visualizer.
4. Code `call_[Method]` - This is essentially a function that calls both the previous two functions, and therefore can be called in one line.
5. Add the function created in Step 4 into `Computation`, in the switch case at the later half of the code. Look to the other methods coded for structure.
   This can also be altered to loop multiple times for different variables.
6. Congratulations, if it gives an output then your method has been incorporated successfully!


## Visualization

### Usage/Syntax:

```matlab
OLCD_Visualization(Methods)
```

`Methods` is a cell array of methods to be visualized.
If benchmark exists, and you want to visualize it as well, then include `'Benchmark'` in the list.

### File Structure

- `Visualization` - This function is the overarching code that plots the network, as well as communities to compare their results.
- `Node_Reorder` - This moves the labeling of communities to the largest being number 1.
  This allows for a somewhat easier view of the communities, as the colours within the visualiser depend on the labels.
- `Node_Sorter` - This reorders the nodes according to the benchmark (so that communities are easier to see).

Unused code (so far)

- `NMI_calc` - Given two community results, will spit out a number between 0 and 1 (NMI calculation).
  Does not work with overlapping communities.
* `ENMI_calc` & `extendNMI_calcs` are functions for calculating the NMI for overlapping communities.
  To use, just input the two community matrices into `ENMI_calc`
* `F1_overlapcalc` extracts the set of overlapping nodes from two overlapping community structures and computes their F1 score using F1_calc.
  To use, input the two community structures into F1_overlapcalc.
