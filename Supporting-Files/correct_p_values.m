function p=correct_p_values(p_val,IDs)

[~,idx]=unique(IDs);
%idx(1)=[]; %remove 0

p=[];
for t=idx,
    p=[p;p_val{t}];
end
% p=unique(p);
p(p>1)=1;
% 
% [~, c_p,FDRp]=fdr_bh(p,0.05,'dep');
% 
% FWERp=bonf_holm(p,0.05); %FWER (Holm, 1979)
% 
% FDR_p_val=p_val;
% FWER_p_val=p_val;
% 
% for t=1:length(p),
%    FDR_p_val{idx(t)}=FDRp(t);
%    FWER_p_val{idx(t)}=FWERp(t);
% end
