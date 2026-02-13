# BnB-Bipart
A balanced bi-paritioning implementation for a netlist represented as a hypergraph using Branch-and-Bound. Takes in a set of blocks connected by a multi-pin net, and defined community contraints (pairs of blocks that prefer to stay together), and splits the block into left and right side with equal size while minimizing a combined cost.
Input netlists or circuits are represented as blocks connected and nets where each block is the first column and is connected to multiple nets. The cost function or the lower bound function is defined as: <br>
Cost = # crossing nets + # community mismatches <br>
The program gurantees global optimal meaning the partition will always have the least cost possible.   
Runtime is reduced by initial placement heuristics, greedy initialization assignment to create an initial upper bound, and a stronger lower bound function by predicting future cost to prune the tree aggressively. 


## Usage
Run the cpp file to create an executable using g++
```bash
g++ bipartition.cpp -o bipartition
```
The executable takes the following arguments
```bash
./bipartition -file_name <path> -init_place <base|fanout>
```
Where init_place creates build order based on ascending block number (base) or descending fanout (number of nets a block is connected to) order. Use the test_circuit* files as input files or files with similar format. 
