function [pdag,p_val,IDs] = PC_with_pval(cond_indep, alpha, k, d, varargin) 
%
% Runs PC with p-values, outputs *raw* p-values
%
% Inputs:
% 1) cond_indep = conditional independence test
% 2) alpha = alpha threshold. If empty, uses default value of 0.20.
% 3) k = conditioning set size. If empty, uses default value of 3.
% 4) d = number of variables
% 5) varargin = stuff necessary for the conditional independence test
% specified with cond_indep (e.g., the data, hyperparameters)
%
% Outputs:
% 1) pdag = the partially directed acyclic graph (PDAG). Edge is
% undirected when pdag(i,j)=pdag(i,j)=1. Edge is directed from i to j when
% pdag(i,j)=-1
% 2) p_val = edge-specific p-values
% 3) IDs = unique identifier per edge-specific hypothesis test
%
% Example call: 
% [pdag,p_val,IDs] = PC_with_pval(@rho_test_PC, 0.20, 3, data);
%
%

if isempty(alpha), alpha=0.20; end
if isempty(k), k=3; end

[G, sep, cell_p, num_struc] = stable_skeleton_discovery(@rho_test_PC, d, k, alpha, varargin{:}); %stable skeleton discovery

[pdag, cell_p2, num_struc2] = get_v_structures2(G, sep, cell_p, num_struc, k, cond_indep, varargin{:}); % get v-structures, unconstrained edge directions

[pdag,G,cell_p2,num_struc2]=clamp_edges(pdag,G,cell_p2,cell_p,num_struc2); %clamp all bidirected edges (refuse to propagate directions with clamped
% edges), revert back to undirected p-values

[pdag, cell_p2, num_struc2] = orientation_rules(G, pdag, cell_p2, num_struc2); %unconstrianed edge propagation via orientation rules

idx=find(pdag==2); %revert back to undirected edges if conflicting rules
for t=1:length(idx),
    [i,j]= ind2sub(size(G),idx(t));
    G(i,j)=1;
    G(j,i)=1;
    pdag(i,j)=1;
    pdag(j,i)=1;
    
    cell_p2{i,j}=cell_p{i,j};
    cell_p2{j,i}=cell_p{j,i};
    num_struc2(i,j)=num_struc(i,j);
    num_struc2(j,i)=num_struc(j,i);
end

p_val=cell_p2;
IDs = num_struc2;