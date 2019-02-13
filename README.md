# PC with p-values (PC-p)

This is a MATLAB implementation of the PC with p-values (PC-p) algorithm. The algorithm returns edge-specific p-values for all edges in the output of PC (i.e., the PDAG) by:

1. Deriving complex edge-specific hypothesis tests by conjoining and disjoining the results of primitive conditional independence tests, and
2. Using modified (a) skeleton discovery, (b) v-structure orientation, and (c) orientation rule application procedures to help ensure that the p-value estimates are valid.

PC-p subsequently estimates and/or controls the false discovery rate (FDR) across every edge, so you don't have to worry about multiple comparison problems.

The associated manuscript is currently under submission (arXiv: http://arxiv.org/abs/1607.03975).

# The Files
There are basically four main functions: 

1. `PC_with_pval.m` -- this is PC-p.
2. `stable_skeleton_discovery.m` -- PC-stable's skeleton discovery procedure with p-values; suitable, for example, when you have time information to automatically orient the edges.
3. `control_FDR.m` -- controls the FDR at a given FDR level q and outputs an FDR corrected graph
4. `estimate_FDR.m` -- estimates the FDR for a given graph

See `pcp_demo.m` for a demo of the functions.

# Reproduce Experimental Results

Run `run_algorithms.m` to collect algorithm outputs

Run `get_results.m` to collect mean FDR, control bias and estimation bias results
