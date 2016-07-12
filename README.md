# PC with p-values (PC-p)

This is a MATLAB implementation of the PC with p-values (PC-p) algorithm. The algorithm returns edge-specific p-values for all edges in the output of PC (i.e., the PDAG) by:

1. Deriving complex edge-specific hypothesis tests by conjoining and disjoining the results of primitive conditional independence tests, and
2. Using modified (a) skeleton discovery, (b) v-structure orientation, and (c) orientation rule application procedures to help ensure that the p-value estimates are valid.

PC-p subsequently controls the false discovery rate across every edge, so you don't have to worry about multiple comparison problems.

# The Files

- PC_with_pval.m is the ``main'' main function.

- See pcp_demo.m for a demo. Also demonstrates how to return p-values for the skeleton, if you don't want to infer causal directions.

Please let me know if you find any bugs by emailing me at ericvonstrobl at google's email dot com.

Coded using scripts from the Bayes Net Toolbox as a base (https://github.com/bayesnet/bnt).
