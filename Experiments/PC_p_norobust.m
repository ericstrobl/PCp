function [pdag,cell_p2,num_struc2,G]=PC_p_norobust(G, sep, cell_p, num_struc, k, cond_indep, varargin)

[pdag, cell_p2, num_struc2] = get_v_structures2_norobust(G, sep, cell_p, num_struc, k, cond_indep, varargin{:}); % get v-structures, unconstrained edge directions
%clamp all bidirected edges (refuse to propagate directions with clamped
%edges), revert back to undirected p-values
[pdag,G,cell_p2,num_struc2]=clamp_edges(pdag,G,cell_p2,cell_p,num_struc2);

[pdag, cell_p2, num_struc2] = orientation_rules_norobust(G, pdag, cell_p2, num_struc2); %unconstrianed edge propagation via orientation rules

idx=find(pdag==2);
for t=1:length(idx),
    [i,j]= ind2sub(size(G),idx(t));
    G(i,j)=1;
    G(j,i)=1;
    pdag(i,j)=1;
    pdag(j,i)=1;
    
    cell_p2{i,j}=max(cell_p{i,j});
    cell_p2{j,i}=max(cell_p{j,i});
    num_struc2(i,j)=num_struc(i,j);
    num_struc2(j,i)=num_struc(j,i);
end