# PC with p-values (PC-p)

This is a MATLAB implementation of the PC with p-values (PC-p) algorithm. The algorithm returns edge-specific p-values for all edges in the output of PC (i.e., the PDAG) by:

1. Deriving complex edge-specific hypothesis tests by conjoining and disjoining the results of primitive conditional independence tests, and
2. Using modified (a) skeleton discovery, (b) v-structure orientation, and (c) orientation rule application procedures to help ensure that the p-value estimates are valid.

PC-p subsequently controls the false discovery rate across every edge, so you don't have to worry about multiple comparison problems.

The associated manuscript is currently under submission.

# The Files
There are basically three functions: 

1. orig_PC_with_pval -- original PC but with p-value computations added,
2. PC_with_pval -- a conservative algorithm that takes extra measures in order to ensure that the p-value bounds are correct with finite sample sizes; best for publication quality results,
3. stable_skeleton_discovery -- just skeleton discovery (no edge orientations) with p-values; suitable for example when you have time information to automatically orient the edges.

See pcp_demo.m for a demo of all three functions. The demo also demonstrates how to control or estimate the FDR for each of the above functions.

Please let me know if you find any bugs or have any suggestions by emailing me at ericvonstrobl at google's email dot com.

Coded using scripts from the Bayes Net Toolbox as a base (https://github.com/bayesnet/bnt).
