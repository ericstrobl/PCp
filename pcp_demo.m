clear all
clc

%load a DAG, generate a Gaussian dataset
load DAG_graph1.mat
% load DAG_graph2.mat
% load DAG_graph3.mat
% load DAG_graph4.mat
% load DAG_graph5.mat

d = size(graph,2);
data=create_dataset_dag(graph,10000,randn(d,d));

%%
%get pdag with raw p-values using Fizher's z test
[pdag,p_val,IDs] = PC_with_pval(@lin_test_PC, [], [], size(data,2),size(data,1),data);

%control the FDR at a given FDR level q
q=0.1;
[pdag_adj,alpha_star] = control_FDR(pdag,p_val, IDs, q)

%estimate the FDR at a given alpha threshold
alpha=0.05;
FDR_est=estimate_FDR(p_val,IDs,alpha)



%%
%get skeleton with raw p-values
[pdag_sk, ~, p_val_sk, IDs_sk] = stable_skeleton_discovery(@lin_test_PC, [], [], size(data,2), size(data,1),data);

%control the FDR at a given FDR level q
q=0.1;
[pdag_adj,alpha_star] = control_FDR(pdag_sk, p_val_sk, IDs_sk, q)

%estimate the FDR at a given alpha threshold
alpha=0.05;
FDR_est=estimate_FDR(p_val_sk, IDs_sk, alpha)



