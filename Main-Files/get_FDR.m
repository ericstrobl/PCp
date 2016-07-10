function FDR_est=get_FDR(p_val,IDs,alpha_thres)
%
% Approximates FDR with Benjamini-Yekutielli FDR estimate
%
% Inputs:
% 1) p_val = edge-specific p-values
% 2) alpha_thres = some alpha threshold (e.g., 0.05)
%
% Outputs:
% 1) FDR_est = FDR estimate
%
% Example call: 
% FDR=get_FDR(p_val,0.05);
%

[~,idx]=unique(IDs);
idx(1)=[]; %remove 0

p=[];
for t=idx,
    p=[p;p_val{t}];
end
% p=unique(p);
p(p>1)=1;

nt=length(p);
den=sum(p<=alpha_thres);
if den==0, den=1; end

FDR_est=(nt*alpha_thres*sum(1./(1:nt)))/den;
