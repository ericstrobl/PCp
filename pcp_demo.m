clear all
clc

%% create a DAG, generate a Gaussian dataset

n=1000; %sample size
d=20; %number of vertices
en=2; %expected neighborhood size
B=create_dag(en,d);
data=create_dataset_dag(B>0,n,B);

%% run PC-p

%get pdag with raw p-values using Fizher's z test
C=corr(data);
[pdag,p_val,IDs] = PC_with_pval(@gaussCItest, [], size(data,2),n,C,n);

%control the FDR at a given FDR level q
q=0.1;
[pdag_adj,alpha_star] = control_FDR(pdag,p_val, IDs, q)

%estimate the FDR at a given alpha threshold
alpha=0.05;
FDR_est = estimate_FDR(p_val,IDs,alpha)



