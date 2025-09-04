# Generalized Bipartite Modularity Density 
Maximization of Generalized Bipartite Modularity Density ($Q_bg$) using bipartite version of RenEEL algorithm


This is an implementation of the Reduced Network Extremal Ensemble Learning (RenEEL) scheme for community detection in bipartite complex networks. The example network used here for illustration is the Southern Women–Event Bipartite Network (Davis, A., Gardner, B., & Gardner, M. (1941). Deep South: A social anthropological study of caste and class. University of Chicago Press).

For comments/questions and reporting any bugs that you encounter in the program please contact Kevin E. Bassler (bassler@uh.edu)

Usage: 
To use the code, follow the steps below:

1. Prepare data

1.1 Use an edgelist file including 3 columns (separated by space or tab with no header) with the last one representing the weight. (If the network is unweighted, just put 1 to the whole column. As the network is bipartite, nodes of the first type should ne indexed from 0 to n and the nodes of the second type are indexed from n+1 onward.
See example in Women_Event.txt)  
1.2 Use bash script (work.sh) to generate the three files required by the program.

Example:
	sh work.sh Women_Event.txt 


2. Compile code

compile main.c, help.c and rg.c with required libraries (math).

Example:
	gcc-9 main.c help.c rg.c  -lm

This will generate file a.out.

3. Run a.out with 3 arguments.

argument 1: Positive Integer, parameter for Randomized Greedy  (usually 2)  
argument 2: Positive Integer, maximum and initial ensemble size of partitions used in RenEEL  
argument 3: Positive Integer, ensemble size of partitions of the reduced network for iteration part in RenEEL  
(seed of random number will be generated using system time)  
argument 4: Real number, for parameter $\chi$ in Generalized Bipartite Modularity Density formula ($\chi$ = 0 for Modularity) (We mainly use $\chi>=0$ to obtain finer structure beyond Modulartiy. But if needed, negative $\chi$ might also be used to obtain coarser sturture than Modularity.)
argument 5: Input file name (edgelist), for example `Women_Event.txt`. The file should have already been processed using `work.sh` as described in 1
  
Example:
	./a.out 2 10 5 1.0 Women_Event.txt

4. Collect results

file 1: `partition_Women_Event.txt` (or partition_[inputfile])
file 2: `results_Women_Event.txt` (or results_[inputfile]), a copy will also be printed to stdout
