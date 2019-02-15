# PC with p-values (PC-p)

This is a MATLAB implementation of the PC with p-values (PC-p) algorithm. The algorithm returns edge-specific p-values for all edges in the output of PC. The algorithm then estimates and/or controls the false discovery rate (FDR) across every edge, so you don't have to worry about multiple comparisons.

The associated manuscript is currently under submission (arXiv: http://arxiv.org/abs/1607.03975).

# Installation

Download all of the files and add all subfolders to your MATLAB path.

# Run PC-p

See `pcp_demo.m` for a demo of the main functions.

# Reproduce Experimental Results

Run `run_algorithms.m` to collect algorithm outputs.

Then run `get_results.m` to collect mean FDR, control bias and estimation bias results.
