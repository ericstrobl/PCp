# PC with p-values (PC-p)

This is a MATLAB implementation of the PC with p-values (PC-p) algorithm. The algorithm returns edge-specific p-values for all edges in the output of PC (i.e., the PDAG) by:

1. Deriving complex edge-specific hypothesis tests by conjoining and disjoining the results of primitive conditional independence tests, and
2. Using modified (a) skeleton discovery, (b) v-structure orientation, and (c) orientation rule application procedures to help ensure that the p-value estimates are valid.

PC-p subsequently controls the false discovery rate (FDR) across every edge, so you don't have to worry about multiple comparison problems.

The associated manuscript is currently under submission (arXiv: http://arxiv.org/abs/1607.03975).

# The Files
There are basically four main functions: 

1. PC_with_pval -- this is PC-p.
2. stable_skeleton_discovery -- PC-stable's skeleton discovery procedure with p-values; suitable, for example, when you have time information to automatically orient the edges.
3. control_FDR -- controls the FDR at a given FDR level q and outputs an FDR corrected graph
4. estimate_FDR -- estimates the FDR for a given graph

See pcp_demo.m for a demo of the functions.

Please let me know if you find any bugs or have any suggestions by emailing me at ericvonstrobl at google's email dot com.

Coded using scripts from the Bayes Net Toolbox as a base (https://github.com/bayesnet/bnt). Tested on MATLAB 2015a.
