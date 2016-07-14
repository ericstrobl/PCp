function [pdag_adj,alpha_star] = control_FDR(pdag, p_val, IDs, q)
%
% Gets alpha* with a pre-specified FDR level q and then controls the FDR of 
% the pdag at level q
%
% Inputs:
% 1) pdag
% 2) p_val = raw edge-specific p-values
% 3) IDs = unique identifier per edge-specific hypothesis test
% 4) q = FDR level

% Outputs:
% 1) pdag_adj = FDR controlled pdag with edges associated with p-values 
% above alpha_star removed
% 2) alpha_star
%
% Example call: 
% [pdag_adj,alpha_star] = control_FDR(pdag,p_val, IDs, 0.05);
%


[~,idx]=unique(IDs);
idx(1)=[]; %remove 0

p=[];
for t=idx,
    p=[p;p_val{t}];
end
% p=unique(p);
p(p>1)=1;

alpha_star=binary_search(1E-100,1,p,q);

idx2 = find(p>alpha_star);

pdag_adj=pdag;

for t=idx2,
    [r_s,c_s] = find(IDs==IDs(idx(t)));
    pdag_adj(r_s,c_s)=0;
end
