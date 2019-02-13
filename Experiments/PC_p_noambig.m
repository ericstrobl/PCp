function [pdag,cell_p2,num_struc2,G]=PC_p_noambig(G, sep, cell_p, num_struc, k, cond_indep, varargin)

[pdag, cell_p2, num_struc2] = get_v_structures2(G, sep, cell_p, num_struc, k, cond_indep, varargin{:}); % get v-structures, unconstrained edge directions
%clamp all bidirected edges (refuse to propagate directions with clamped
%edges), revert back to undirected p-values
% [pdag,G,cell_p,num_struc]=clamp_edges(pdag,G,cell_p2,cell_p,num_struc);

[pdag, cell_p2, num_struc2] = orientation_rules_nolock(G, pdag, cell_p2, num_struc2); %unconstrianed edge propagation via orientation rules

idx=find(pdag==-1);
for t=1:length(idx),
    [r,c]=ind2sub(size(pdag),idx(t));
    if pdag(c,r)==-1;
        pdag(c,r)=1;
        pdag(r,c)=1;
        
        cell_p2{r,c}=max(cell_p{r,c});
        cell_p2{c,r}=max(cell_p{c,r});
        
        num_struc2(r,c)=num_struc(r,c);
        num_struc2(c,r)=num_struc(c,r);
    else
        pdag(c,r)=0;
        cell_p2{c,r}=[];
        num_struc2(c,r)=0;
    end
        
    end
end