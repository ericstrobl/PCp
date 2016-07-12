# PCp

This is a MATLAB implementation of PC with p-values (PC-p). The algorithm returns edge-specific p-values for all edges in the output of PC (the PDAG) by:

1. Deriving complex edge-specific hypothesis tests by conjoining and disjoining the results of primitive conditional independence tests, and
2. Using modified (a) skeleton discovery, (b) v-structure orientation, and (c) orientation rule application procedures. 

PC-p subsequently controls the false discovery rate across every edge. 

This work was inspired by the idea of getting at a "causal p-value."

- PC_with_pval.m is the main file.

- See pcp_demo.m for a demo.

Coded using scripts from the Bayes Net Toolbox as a base (https://github.com/bayesnet/bnt).
