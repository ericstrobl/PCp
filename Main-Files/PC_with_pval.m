function [pdag,p_val,IDs] = PC_with_pval(cond_indep, alpha, k, d, n, varargin) 
%
% Runs PC with p-values, outputs *raw* p-values. This is a conservative
% algorithm. Use orig_PC_with_pval for a less conservative algorithm 
% (it is still theoretically correct).
%
% Inputs:
% 1) cond_indep = conditional independence test function handle. The test
% should be of the form CItest(x,y,Z,varargin).
% 2) alpha = alpha threshold. If empty, uses heuristic default value that 
% starts at 0.20 and then decreases at rate root n, where n denotes sample 
% size; this is again just a heuristic but decreasing the alpha value
% ensures asymptotic consistency. If you want to be extra conservative, set 
% alpha=0.20 no matter what sample size. The theory says that the alpha
% threshold should be high *for the problem at hand*, so the Type II error 
% rate is low.
% 3) k = maximum conditioning set size. If empty, algorithm dynamically 
% adjusts the conditioning set size for each node according to whatever 
% size is necessary. Set this if you know the maximum in-degree.
% 4) d = number of variables
% 5) n = number of samples
% 6) varargin = stuff necessary for the conditional independence test
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
% [pdag,p_val,IDs] = PC_with_pval(@lin_test_PC, [], [], size(data,2),size(data,1),data);
%
% Code written by Eric V. Strobl using the Bayes Net Toolbox
% (https://github.com/bayesnet/bnt) as a base.
%

if isempty(alpha),
%     alpha=0.20;
    if n<=500, alpha=0.20;
    else alpha=0.20/sqrt(n/500);
    end
end

[G, sep, cell_p, num_struc] = stable_skeleton_discovery(cond_indep, alpha, k, d, n, varargin{:}); %stable skeleton discovery

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