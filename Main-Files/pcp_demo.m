clear all
clc

%load a DAG, generate a Gaussian dataset
load DAG_graph1.mat

d = size(graph,2);
data=create_dataset_dag(graph,1000,randn(d,d));

%%
%get pdag with raw p-values using Spearman's rho
[pdag,p_val,IDs] = PC_with_pval(@rho_test_PC, [], [], size(data,2),data);

%get alpha_star using a given FDR level q
q=0.10;
alpha_star = get_alpha_star(p_val, IDs, q)

%approximate FDR for a given graph using a given alpha_threshold
alpha=0.05;
FDR_est=get_FDR(p_val,IDs,alpha)



%%
%get skeleton with raw p-values
[pdag_sk, sep_sk, p_val_sk, IDs_sk] = stable_skeleton_discovery(@rho_test_PC, [], [], size(data,2), data);

%get alpha_star using a given FDR level q
q=0.10;
alpha_star = get_alpha_star(p_val_sk, IDs_sk, q)

%approximate FDR for a given graph using a given alpha_threshold
alpha=0.05;
FDR_est=get_FDR(p_val_sk,IDs_sk,alpha)


