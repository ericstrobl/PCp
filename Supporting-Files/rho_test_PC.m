function [p,r]=rho_test_PC(x,y,z,data)

if isempty(z),
    [r,p]=corr(data(:,x),data(:,y),'Type','Spearman'); %%%
else
    [r,p]=partialcorr(data(:,x),data(:,y),data(:,z),'Type','Spearman'); %%%
end