# PC with p-values (PC-p)

This is a MATLAB implementation of the PC with p-values (PC-p) algorithm. The algorithm returns edge-specific p-values for all edges in the output of PC (i.e., the PDAG) by:

1. Deriving complex edge-specific hypothesis tests by conjoining and disjoining the results of primitive conditional independence tests, and
2. Using modified (a) skeleton discovery, (b) v-structure orientation, and (c) orientation rule application procedures to help ensure that the p-value estimates are valid.

PC-p subsequently controls the false discovery rate (FDR) across every edge, so you don't have to worry about multiple comparison problems.

The associated manuscript is currently under submission.

# The Files
There are basically two main functions: 

1. PC_with_pval -- this is PC-p,
3. stable_skeleton_discovery -- PC-stable's skeleton discovery procedure (no edge orientations) with p-values; suitable for example when you have time information to automatically orient the edges.

See pcp_demo.m for a demo of both functions. The demo also demonstrates how to control or estimate the FDR.

Please let me know if you find any bugs or have any suggestions by emailing me at ericvonstrobl at google's email dot com.

Coded using scripts from the Bayes Net Toolbox as a base (https://github.com/bayesnet/bnt). Tested on MATLAB 2015a.
