# Empirical comparison of overlapping community detection algorithms

These modules are separated into the computation and visualisation modules, where computation does all the calculations, and visualisation allows you to view the results and compare them.

## Computation Module Information

### Usage/Syntax:

```
Final = Computation(Input, Methods, isBenchmark, benchFileName)
```

Where:

* `Final` -	The outputted structure, which contains all data from the algorithms (Method name, community result)
* `Input` -	Can be an adjacency matrix, matrix list of undirected links, or matrix list of directed links
* `Methods` -	Cells of the names of methods to be computed. *Default is Jerry and Shen methods*
* `isBenchmark` -	0 for no benchmark, 1 for benchmark. *Default is 0*
* `benchFileName` - String name of the text file within which the benchmark communities are defined

This function also saves the Final structure as a file called `Computation_Result.mat`. This file will be accessed by the Visualisation code.

### File structure

These files are organised as such:

Within the `Computation_Module`, there are many directories.
Most of these are names of methods currently implemented.
Within these directories, there are 3 files -

* `call_[MethodName]` - This function is called by the overarching `computation.m` function. It simplifies the entire process down to one line.
* `run_[MethodName]` - This function is called to do the actual calculations; running the algorithm and spitting out its output (Note: If calculations need to be hard coded to specific directories, this is done within this `.m` file. For examples, view `run_Gopalan`).
* `process_[MethodName]` - Called after running the algorithm, as this function will take the algorithm output and convert it to the universal matrix format, explained below.

Apart from these directories, there are two others called `Conversions` and `SourceCode`.

* `Conversions` is a directory that holds the `.m` files that convert any of the 3 types of inputs into other types. This directory should not be changed, unless there is a different type of input that is not adjacency matrix, undirected list or directed list.
* `SourceCode` is the directory where all the algorithms are stored within their own folder.

Finally, the main overarching function, `Computation`.
This function is the one to call (as seen above) when the computation of algorithms is needed.

### Algorithm Output ###

After the `process_[Method]` function is run, the output should be in a matrix format, where each row is a node, and each column is a community.

![](http://i.imgur.com/LP5r46R.png?1)

In the above example:

* Node 1 is in community 3
* Node 2 is in community 1
* Node 3 is in both communities 1 and 3
* Node 4 is in all three communities
* Node 5 is homeless (not linked to any communities)

Nodes can also be more strongly linked to a community, but the sum of numbers for each node must be 1.

### Adding a new method

If you wish to add a module, here is how to do so:

1. Place source code of the algorithm within its own named directory inside  the `Computation_Module/SourceCode` directory.
2. Code `run_[Method]` - This should take the network (in whatever format that is needed) and parameters, and output a final community result.
   This output does not have to be in the matrix format, the processing function should do that.
3. Code `process_[Method]` - This will take the output of the previous step as an input, and then output the matrix format for the visualiser.
4. Code `call_[Method]` - This is essentially a function that calls both the previous two functions, and therefore can be called in one line.
5. Add the function created in step 4 into `Computation`, in the switch case at the later half of the code. Look to the other methods coded for structure.
   This can also be altered to loop multiple times for different variables.
6. Congratulations, assuming it gives an output, your method is now implemented!


## Visualisation Module Information ##

### Usage/Syntax: ###

```
Visualisation(Methods)
```

Where:

* Methods - Cell list of methods to be visualised.
  If benchmark exists, and you want to visualise it as well, the name 'Benchmark' needs to be within this cell list.
  *If no Methods are supplied, then the default is everything within the results of the Computation module*

### File Structure ###

Within the `Visualisation_Module` directory, there are only `.m` files.

* `Visualisation` - This function is the overarching code that plots the network, as well as communities to compare their results.
* `Node_Reorder` - This moves the labeling of communities to the largest being number 1.
  This allows for a somewhat easier view of the communities, as the colours within the visualiser depend on the labels.
* `Node_Sorter` - This reorders the nodes according to the benchmark (so that communities are easier to see).

Unused code (So far)

* `NMI_calc` - Given two community results, will spit out a number between 0 and 1 (NMI calculation).
  Does not work with overlapping communities.
* `ENMI_calc` & `extendNMI_calcs` are functions for calculating the NMI for overlapping communities.
  To use, just input the two community matrices into `ENMI_calc`

### Adding a new method

Because of the fact that the visualiser runs off the output of the computation module, nothing new needs to be added into this code!
