# PCp

This is a MATLAB implementation of PC with p-values (PC-p). The algorithm returns edge-specific p-values for all edges in the CPDAG by:
1. Decomposing complex edge-specific hypothesis tests into primitive conditional independence tests,
2. Using order-independent (a) skeleton discovery (as in PC-stable), (b) v-structure orientation, and (c) orientation rule application. This work was inspired by the idea of getting at a "causal p-value."

- PC_with_pval.m is the main file

- See pcp_demo.m for a demo

Coded using the Bayes Net Toolbox as a base (https://github.com/bayesnet/bnt).
