# Community detection & visualisation software
Brandon Lam, 2015

These modules are separated into the computation and visualisation modules, where computation does all the calculations, and visualisation allows you to view the results and compare.

[File structure]()

[File](#fs)

##Computation Module Information##
Usage/Syntax:

```
Final = Computation(Input, Methods, isBenchmark, benchfilename)
```

Where:

* Final -	The outputted structure, which contains all data from the algorithms (Method name, community result)
* Input -	Can be an adjacency matrix, matrix list of undirected links, or matrix list of directed links
* Methods -	Cells of the names of methods to be computed
* isBenchmark -	0 for no benchmark, 1 for benchmark
* benchfilename - String name of the text file within which the benchmark communities are defined

<a name="FS"/></a>
### File structure ###


These files are organised as such:

Within the `/Computation_Module`, there are many directories. Most of these are names of methods currently implemented. Within these directories, there are 3 files - 

* `call_[MethodName].m` - This function is called by the overarching `computation.m` function. It simplifies the entire process down to one line.
* `run_[MethodName].m` - This function is called to do the actual calculations; running the algorithm and spitting out its output (Note: If calculations need to be hard coded to specific directories, this is done within this m file. For examples, view `run_Gopalan.m`)
* `process_[MethodName].m` - Called after running the algorithm, as this function will take the algorithm output and convert it to the universal matrix format, explained below.



### Adding a new method

If you wish to add a module, you ought to add the following files in the following formats:

### Summary of file structure

The repository is organized as follows:


*Comp module*:
* `process_Gopalan.m` this file does...


*Vis module*:
* `process_Gopalan.m` this file ...
