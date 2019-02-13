# PC with p-values (PC-p)

This is a MATLAB implementation of the PC with p-values (PC-p) algorithm. The algorithm returns edge-specific p-values for all edges in the output of PC (i.e., the PDAG) by:

1. Deriving complex edge-specific hypothesis tests by conjoining and disjoining the results of primitive conditional independence tests, and
2. Using modified (a) skeleton discovery, (b) v-structure orientation, and (c) orientation rule application procedures to help ensure that the p-value estimates are valid.

PC-p subsequently estimates and/or controls the false discovery rate (FDR) across every edge, so you don't have to worry about multiple comparison problems.

The associated manuscript is currently under submission (arXiv: http://arxiv.org/abs/1607.03975).

# Installation

Download all of the files and add all subfolders to your MATLAB path.

# Run PC-p
There are basically four main functions: 

> [pdag,p_val,IDs] = PC_with_pval(@gaussCItest, [], [], size(data,2),size(data,1),data);

See `pcp_demo.m` for a demo of the functions.

# Reproduce Experimental Results

Run `run_algorithms.m` to collect algorithm outputs

Run `get_results.m` to collect mean FDR, control bias and estimation bias results
