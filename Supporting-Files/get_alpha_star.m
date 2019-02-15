function alpha_star = get_alpha_star(p_val, IDs, q)
%
% Gets alpha* with a pre-specified FDR level q
%
% Inputs:
% 1) p_val = edge-specific p-values
% 2) IDs = unique identifier per edge-specific hypothesis test
% 3) q = FDR level

% Outputs:
% 1) alpha_star
%
% Example call: 
% alpha_star = get_alpha_star(0.10, p_val, IDs);
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
