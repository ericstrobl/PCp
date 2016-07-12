# PCp

This is a MATLAB implementation of PC with p-values (PC-p). The algorithm returns edge-specific p-values for all edges in the output of PC (the PDAG) by:

1. Deriving complex edge-specific hypothesis tests by conjoining and disjoining the results of primitive conditional independence tests, and
2. Using modified (a) skeleton discovery, (b) v-structure orientation, and (c) orientation rule application procedures. 

PC-p subsequently controls the false discovery rate across every edge. I completed this project because I really wanted to get at a "causal p-value"; I hope that you also believe the algorithm helps achieves this goal.

# The Files

- PC_with_pval.m is the ``main'' main file.

- See pcp_demo.m for a demo.

- I have included 5 example DAGs that you can play with.

Please let me know if you find any bugs by emailing me at ericvonstrobl at google's email dot com.

Coded using scripts from the Bayes Net Toolbox as a base (https://github.com/bayesnet/bnt).
