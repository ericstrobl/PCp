clear all
clc

%load a DAG, generate a Gaussian dataset
load DAG_graph1.mat

d = size(graph,2);
data=create_dataset_dag(graph,1000,randn(d,d));

%get pdag with raw p-values using Spearman's rho
[pdag,p_val,IDs] = PC_with_pval(@rho_test_PC, [], [], size(data,2),data);

%get alpha_star using a given FDR level q
alpha_star = get_alpha_star(p_val, IDs, 0.10);

%approximate FDR for a given graph using a given alpha_threshold
FDR_est=get_FDR(p_val,IDs,0.05);

%get FDR and FWER corrected p-values
[FDR_p_val,FWER_p_val]=correct_p_values(p_val,IDs);





