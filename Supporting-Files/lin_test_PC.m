function [p,r]=lin_test_PC(x,y,z,data)

if isempty(z),
    [r,p]=corr(data(:,x),data(:,y)); %%%
else
    [r,p]=partialcorr(data(:,x),data(:,y),data(:,z)); %%%
end